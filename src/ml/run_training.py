import argparse
import csv
import datetime
import os
import time
from contextlib import nullcontext
from typing import List

import matplotlib
import torch
import torch.distributed as dist
import torch.multiprocessing as mp
from torch.nn.parallel import DistributedDataParallel as DDP

from ml.train import (
    get_training_objects,
    load_config,
    run_batch_diffusion,
    set_seed,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def save_checkpoint(ckpt_dir, step, model, optimizer, cfg):
    os.makedirs(ckpt_dir, exist_ok=True)
    ckpt_name = f"diffusion_step_{step:08d}.pt"
    ckpt_path = os.path.join(ckpt_dir, ckpt_name)
    model_to_save = model.module if hasattr(model, "module") else model
    torch.save(
        {
            "step": step,
            "model_state_dict": model_to_save.state_dict(),
            "optimizer_state_dict": optimizer.state_dict(),
            "config": cfg,
        },
        ckpt_path,
    )
    print(f"saved checkpoint: {ckpt_path}")
    return ckpt_path


def load_checkpoint(ckpt_path, model, optimizer, device, rank=0):
    if not ckpt_path:
        return 0
    if not os.path.exists(ckpt_path):
        raise FileNotFoundError(f"Checkpoint not found: {ckpt_path}")

    ckpt = torch.load(ckpt_path, map_location=device, weights_only=False)
    raw_model_state = ckpt.get("model_state_dict")
    if not isinstance(raw_model_state, dict):
        raise RuntimeError("Checkpoint is missing a valid model_state_dict.")

    model_state = model.state_dict()
    if not any(k in model_state for k in raw_model_state.keys()):
        if any(k.startswith("module.") for k in raw_model_state.keys()):
            raw_model_state = {
                (k[7:] if k.startswith("module.") else k): v
                for k, v in raw_model_state.items()
            }
            if rank == 0:
                print("checkpoint load: stripped leading 'module.' prefix")

    loadable_state = {}
    unexpected_keys = []
    shape_mismatches = []
    for k, v in raw_model_state.items():
        if k not in model_state:
            unexpected_keys.append(k)
            continue
        if model_state[k].shape != v.shape:
            shape_mismatches.append((k, tuple(v.shape), tuple(model_state[k].shape)))
            continue
        loadable_state[k] = v

    incompatible = model.load_state_dict(loadable_state, strict=False)
    missing_keys = list(incompatible.missing_keys)
    unexpected_keys.extend(list(incompatible.unexpected_keys))

    if rank == 0:
        print(
            "checkpoint load summary: "
            f"loaded={len(loadable_state)} "
            f"missing={len(missing_keys)} "
            f"unexpected={len(unexpected_keys)} "
            f"shape_mismatch={len(shape_mismatches)}"
        )
        for k, ckpt_shape, model_shape in shape_mismatches:
            print(f"[ckpt][shape_mismatch] {k}: ckpt={ckpt_shape} model={model_shape}")

    can_restore_optimizer = (
        len(missing_keys) == 0
        and len(unexpected_keys) == 0
        and len(shape_mismatches) == 0
    )
    if "optimizer_state_dict" in ckpt and ckpt["optimizer_state_dict"] is not None:
        if can_restore_optimizer:
            optimizer.load_state_dict(ckpt["optimizer_state_dict"])
            if rank == 0:
                print("checkpoint load: optimizer state loaded")
        elif rank == 0:
            print("checkpoint load: optimizer state skipped after partial model restore")

    step = int(ckpt.get("step", 0))
    if rank == 0:
        print(f"resumed from checkpoint: {ckpt_path} (step={step})")
    return step


def parse_devices(devices_arg: str) -> List[int]:
    if not devices_arg:
        return []
    return [int(x.strip()) for x in devices_arg.split(",") if x.strip()]


def append_metric(metrics_path, row):
    metrics_dir = os.path.dirname(metrics_path)
    if metrics_dir:
        os.makedirs(metrics_dir, exist_ok=True)
    exists = os.path.exists(metrics_path)
    fields = ["step", "split", "loss", "lr", "steps_per_sec"]
    with open(metrics_path, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        if not exists:
            writer.writeheader()
        writer.writerow({field: row.get(field, "") for field in fields})


def save_probability_figure(fig, vis_dir, step):
    os.makedirs(vis_dir, exist_ok=True)
    path = os.path.join(vis_dir, f"all_masked_probability_step_{step:08d}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def run_training_loop(cfg, rank=0, world_size=1, devices=None, distributed=False, resume_ckpt=None):
    seed = int(cfg.get("seed", 123))
    set_seed(seed + rank)
    if torch.cuda.is_available():
        torch.backends.cudnn.benchmark = True
        if hasattr(torch.backends, "cuda") and hasattr(torch.backends.cuda, "matmul"):
            torch.backends.cuda.matmul.allow_tf32 = True
        if hasattr(torch.backends, "cudnn") and hasattr(torch.backends.cudnn, "allow_tf32"):
            torch.backends.cudnn.allow_tf32 = True
        if hasattr(torch, "set_float32_matmul_precision"):
            torch.set_float32_matmul_precision("high")

    devices = devices or []
    data_cfg = cfg["data"]
    model_cfg = cfg["model"]
    train_cfg = cfg["train"]

    if distributed:
        if not torch.cuda.is_available():
            raise RuntimeError("DDP requested but CUDA is not available.")
        device_id = int(devices[rank])
        torch.cuda.set_device(device_id)
        device = torch.device(f"cuda:{device_id}")
    else:
        if len(devices) == 1:
            device = torch.device(f"cuda:{int(devices[0])}")
        else:
            device = torch.device(train_cfg.get("device", "cpu"))

    train_loader, val_loader, model, optimizer = get_training_objects(
        data_cfg,
        model_cfg,
        train_cfg,
        device,
        distributed=distributed,
        rank=rank,
        world_size=world_size,
    )

    step = load_checkpoint(resume_ckpt, model, optimizer, device, rank=rank) if resume_ckpt else 0

    if distributed:
        model = DDP(model, device_ids=[device.index], output_device=device.index)

    max_steps = int(train_cfg["max_steps"])
    val_every_steps = int(train_cfg["val_every_steps"])
    max_val_batches = int(train_cfg.get("max_val_batches", 0))
    val_visual_every = int(train_cfg.get("val_visual_every", 1))
    log_every = int(train_cfg["log_every"])
    ckpt_every_n_steps = int(train_cfg.get("ckpt_every_n_steps", 0))
    grad_accum_steps = max(1, int(train_cfg.get("grad_accum_steps", 1)))
    max_grad_norm = float(train_cfg.get("max_grad_norm", 0.0) or 0.0)
    target_lr = float(train_cfg.get("lr", 1e-3))
    warmup_steps = max(0, int(train_cfg.get("warmup_steps", 0)))
    warmup_start_lr = float(train_cfg.get("warmup_start_lr", target_lr))

    run_dir = train_cfg.get("run_dir", "runs/diffusion")
    ckpt_dir = train_cfg.get("ckpt_dir") or os.path.join(run_dir, "checkpoints")
    vis_dir = train_cfg.get("vis_dir") or os.path.join(run_dir, "probability_maps")
    metrics_path = train_cfg.get("metrics_path") or os.path.join(run_dir, "metrics.csv")
    is_master = rank == 0

    last_speed_step = int(step)
    last_speed_time = time.perf_counter()
    last_checkpoint_step = -1
    last_validation_step = -1

    train_batches_per_epoch = max(1, len(train_loader))
    consumed_train_batches = int(step) * grad_accum_steps
    sampler_epoch = consumed_train_batches // train_batches_per_epoch
    if consumed_train_batches % train_batches_per_epoch > 0:
        sampler_epoch += 1
    if distributed and hasattr(train_loader, "sampler") and hasattr(train_loader.sampler, "set_epoch"):
        train_loader.sampler.set_epoch(sampler_epoch)

    train_iter = iter(train_loader)
    model.train()
    accum_step_idx = 0

    while step < max_steps:
        if warmup_steps > 0 and step < warmup_steps:
            progress = float(step) / float(max(1, warmup_steps))
            current_lr = warmup_start_lr + progress * (target_lr - warmup_start_lr)
        else:
            current_lr = target_lr
        for param_group in optimizer.param_groups:
            param_group["lr"] = current_lr

        try:
            batch = next(train_iter)
        except StopIteration:
            sampler_epoch += 1
            if distributed and hasattr(train_loader, "sampler") and hasattr(train_loader.sampler, "set_epoch"):
                train_loader.sampler.set_epoch(sampler_epoch)
            train_iter = iter(train_loader)
            batch = next(train_iter)

        accum_end = ((accum_step_idx + 1) % grad_accum_steps) == 0
        sync_ctx = (
            model.no_sync()
            if distributed and grad_accum_steps > 1 and not accum_end and hasattr(model, "no_sync")
            else nullcontext()
        )
        with sync_ctx:
            step, metrics = run_batch_diffusion(
                step,
                batch,
                model,
                optimizer,
                device,
                training=True,
                grad_accum_steps=grad_accum_steps,
                accum_step_idx=accum_step_idx,
                max_grad_norm=max_grad_norm,
            )
        accum_step_idx += 1

        if is_master and metrics["optimizer_step"]:
            steps_per_sec = ""
            if log_every > 0 and step % log_every == 0:
                now = time.perf_counter()
                delta_steps = max(step - last_speed_step, 1)
                delta_time = max(now - last_speed_time, 1e-9)
                steps_per_sec = delta_steps / delta_time
                last_speed_step = step
                last_speed_time = now
                print(
                    f"step {step} train_loss {metrics['loss']:.6f} "
                    f"lr {current_lr:.6g} steps_per_sec {steps_per_sec:.3f}"
                )
            append_metric(
                metrics_path,
                {
                    "step": step,
                    "split": "train",
                    "loss": metrics["loss"],
                    "lr": current_lr,
                    "steps_per_sec": steps_per_sec,
                },
            )

        if (
            is_master
            and ckpt_every_n_steps > 0
            and step > 0
            and step % ckpt_every_n_steps == 0
            and step != last_checkpoint_step
        ):
            save_checkpoint(ckpt_dir, step, model, optimizer, cfg)
            last_checkpoint_step = step

        if step > 0 and step % val_every_steps == 0 and step != last_validation_step:
            if distributed:
                dist.barrier()
            if is_master:
                eval_model = model.module if hasattr(model, "module") else model
                eval_model.eval()
                val_round = step // val_every_steps
                do_visual = (val_round % max(1, val_visual_every)) == 0
                val_losses = []
                with torch.no_grad():
                    for i, val_batch in enumerate(val_loader):
                        _step, val_metrics = run_batch_diffusion(
                            step,
                            val_batch,
                            eval_model,
                            optimizer,
                            device,
                            training=False,
                            visualize=(do_visual and i == 0),
                        )
                        val_losses.append(val_metrics["loss"])
                        fig = val_metrics.get("figure")
                        if fig is not None:
                            out_path = save_probability_figure(fig, vis_dir, step)
                            print(f"saved probability map: {out_path}")
                        if max_val_batches > 0 and (i + 1) >= max_val_batches:
                            break
                if val_losses:
                    val_loss = float(sum(val_losses) / len(val_losses))
                    print(f"step {step} val_loss {val_loss:.6f}")
                    append_metric(
                        metrics_path,
                        {
                            "step": step,
                            "split": "val",
                            "loss": val_loss,
                            "lr": current_lr,
                            "steps_per_sec": "",
                        },
                    )
                model.train()
                last_validation_step = step
            if distributed:
                dist.barrier()

    if is_master:
        save_checkpoint(ckpt_dir, step, model, optimizer, cfg)
        print("training complete")


def ddp_worker(rank, world_size, devices, cfg, resume_ckpt, ddp_timeout_minutes):
    timeout = datetime.timedelta(minutes=int(ddp_timeout_minutes))
    dist.init_process_group(
        backend="nccl",
        rank=rank,
        world_size=world_size,
        timeout=timeout,
    )
    try:
        run_training_loop(
            cfg=cfg,
            rank=rank,
            world_size=world_size,
            devices=devices,
            distributed=True,
            resume_ckpt=resume_ckpt,
        )
    finally:
        dist.destroy_process_group()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=str,
        default=os.path.join(os.path.dirname(__file__), "configs_diff.yaml"),
    )
    parser.add_argument(
        "--devices",
        type=str,
        default="",
        help="Comma-separated CUDA device ids, for example 0 or 0,1.",
    )
    parser.add_argument("--master_port", type=int, default=29500)
    parser.add_argument("--batch_size", type=int, default=None)
    parser.add_argument("--base_path", type=str, default="")
    parser.add_argument("--run_dir", type=str, default="")
    parser.add_argument("--resume_ckpt", type=str, default="")
    parser.add_argument("--ddp_timeout_minutes", type=int, default=60)
    args = parser.parse_args()

    cfg = load_config(args.config)
    if args.batch_size is not None:
        cfg.setdefault("data", {})["batch_size"] = int(args.batch_size)
    if args.base_path:
        cfg.setdefault("data", {})["base_path"] = args.base_path
    if args.run_dir:
        cfg.setdefault("train", {})["run_dir"] = args.run_dir
    devices = parse_devices(args.devices)

    if len(devices) <= 1:
        run_training_loop(
            cfg=cfg,
            rank=0,
            world_size=1,
            devices=devices,
            distributed=False,
            resume_ckpt=args.resume_ckpt or None,
        )
        return

    if not torch.cuda.is_available():
        raise RuntimeError("Requested multi-GPU DDP but CUDA is not available.")
    visible = torch.cuda.device_count()
    for d in devices:
        if d < 0 or d >= visible:
            raise RuntimeError(
                f"Requested device {d} not in visible CUDA range [0, {visible - 1}]."
            )

    os.environ.setdefault("MASTER_ADDR", "127.0.0.1")
    os.environ["MASTER_PORT"] = str(args.master_port)

    world_size = len(devices)
    mp.spawn(
        ddp_worker,
        args=(world_size, devices, cfg, args.resume_ckpt or None, args.ddp_timeout_minutes),
        nprocs=world_size,
        join=True,
    )


if __name__ == "__main__":
    main()
