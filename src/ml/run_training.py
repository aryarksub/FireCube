import argparse
import os
from typing import List

import torch
import torch.distributed as dist
import torch.multiprocessing as mp
from torch.nn.parallel import DistributedDataParallel as DDP

import wandb

from ml.train import (
    load_config,
    set_seed,
    get_training_objects,
    run_batch_ce,
    run_batch_diffusion,
)


def save_checkpoint(ckpt_dir, step, model, optimizer, cfg, loss_type, wandb_id):
    os.makedirs(ckpt_dir, exist_ok=True)
    run_tag = wandb_id if wandb_id else "no_wandb"
    ckpt_name = f"{loss_type}_{run_tag}_step_{step:08d}.pt"
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


def parse_devices(devices_arg: str) -> List[int]:
    if not devices_arg:
        return []
    return [int(x.strip()) for x in devices_arg.split(",") if x.strip()]


def run_training_loop(cfg, rank=0, world_size=1, devices=None, distributed=False):
    seed = int(cfg.get("seed", 123))
    set_seed(seed + rank)

    devices = devices or []
    data_cfg = cfg["data"]
    model_cfg = cfg["model"]
    train_cfg = cfg["train"]
    wandb_cfg = cfg.get("wandb", {})

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

    if distributed:
        model = DDP(model, device_ids=[device.index], output_device=device.index)

    max_steps = int(train_cfg["max_steps"])
    val_every_steps = int(train_cfg["val_every_steps"])
    max_val_batches = int(train_cfg.get("max_val_batches", 0))  # 0 -> full val loader
    val_visual_every = int(train_cfg.get("val_visual_every", 1))
    log_every = int(train_cfg["log_every"])
    log_images_every = int(train_cfg["log_images_every"])
    sample_steps = int(train_cfg["sample_steps"])
    loss_type = train_cfg.get("loss_type", "diffusion")
    show_pred_vs_gt = bool(train_cfg.get("show_pred_vs_gt", loss_type != "diffusion"))
    ckpt_every_n_steps = int(train_cfg.get("ckpt_every_n_steps", 0))
    ckpt_dir = train_cfg.get("ckpt_dir", "checkpoints")
    is_master = rank == 0
    local_log_every = log_every if is_master else (max_steps + 1)
    local_log_images_every = log_images_every if is_master else 0

    use_wandb = bool(wandb_cfg.get("enabled", False))
    wandb_run = None
    if use_wandb and is_master:
        wandb_run = wandb.init(
            project=wandb_cfg.get("project", "firecube"),
            name=wandb_cfg.get("run_name"),
            mode=wandb_cfg.get("mode", "online"),
            config=cfg,
        )
    wandb_id = wandb_run.id if wandb_run is not None else None

    run_batch = run_batch_diffusion if loss_type == "diffusion" else run_batch_ce

    step = 0
    sampler_epoch = 0
    if distributed and hasattr(train_loader, "sampler") and hasattr(train_loader.sampler, "set_epoch"):
        train_loader.sampler.set_epoch(sampler_epoch)

    train_iter = iter(train_loader)
    model.train()
    while step < max_steps:
        try:
            batch = next(train_iter)
        except StopIteration:
            sampler_epoch += 1
            if distributed and hasattr(train_loader, "sampler") and hasattr(train_loader.sampler, "set_epoch"):
                train_loader.sampler.set_epoch(sampler_epoch)
            train_iter = iter(train_loader)
            batch = next(train_iter)

        step = run_batch(
            step,
            batch,
            model,
            optimizer,
            wandb_run,
            device,
            local_log_every,
            local_log_images_every,
            training=True,
            show_pred_vs_gt=show_pred_vs_gt,
            is_master=is_master,
        )

        if is_master and ckpt_every_n_steps > 0 and step > 0 and (step % ckpt_every_n_steps) == 0:
            ckpt_path = save_checkpoint(
                ckpt_dir,
                step,
                model,
                optimizer,
                cfg,
                loss_type=loss_type,
                wandb_id=wandb_id,
            )
            if wandb_run is not None:
                wandb_run.log({"checkpoint/step": step, "checkpoint/path": ckpt_path})

        if step > 0 and step % val_every_steps == 0:
            if distributed:
                dist.barrier()
            if is_master:
                # Avoid rank-0-only DDP forward during validation/sampling:
                # DDP can trigger collectives (e.g., buffer sync) and deadlock
                # while other ranks are waiting at barriers.
                eval_model = model.module if hasattr(model, "module") else model
                eval_model.eval()
                val_round = step // val_every_steps
                do_val_visual = (val_round % max(1, val_visual_every)) == 0
                with torch.no_grad():
                    for i, val_batch in enumerate(val_loader):
                        run_batch(
                            step,
                            val_batch,
                            eval_model,
                            optimizer,
                            wandb_run,
                            device,
                            local_log_every,
                            local_log_images_every,
                            training=False,
                            visualize=(loss_type == "diffusion" and do_val_visual and i == 0),
                            sample_steps=sample_steps,
                            show_pred_vs_gt=show_pred_vs_gt,
                            is_master=is_master,
                        )
                        if max_val_batches > 0 and (i + 1) >= max_val_batches:
                            break
                model.train()
            if distributed:
                dist.barrier()

    if is_master:
        final_ckpt_path = save_checkpoint(
            ckpt_dir,
            step,
            model,
            optimizer,
            cfg,
            loss_type=loss_type,
            wandb_id=wandb_id,
        )
        if wandb_run is not None:
            wandb_run.log({"checkpoint/final_path": final_ckpt_path, "checkpoint/final_step": step})

    if wandb_run is not None:
        wandb_run.finish()
    if is_master:
        print("training complete")


def ddp_worker(rank, world_size, devices, cfg):
    dist.init_process_group(backend="nccl", rank=rank, world_size=world_size)
    try:
        run_training_loop(
            cfg=cfg,
            rank=rank,
            world_size=world_size,
            devices=devices,
            distributed=True,
        )
    finally:
        dist.destroy_process_group()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=str,
        default=os.path.join(os.path.dirname(__file__), "configs.yaml"),
    )
    parser.add_argument(
        "--devices",
        type=str,
        default="",
        help="Comma-separated CUDA device ids (e.g. 0 or 0,1).",
    )
    parser.add_argument("--master_port", type=int, default=29500)
    args = parser.parse_args()

    cfg = load_config(args.config)
    devices = parse_devices(args.devices)

    if len(devices) <= 1:
        run_training_loop(cfg=cfg, rank=0, world_size=1, devices=devices, distributed=False)
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
    mp.spawn(ddp_worker, args=(world_size, devices, cfg), nprocs=world_size, join=True)


if __name__ == "__main__":
    main()
