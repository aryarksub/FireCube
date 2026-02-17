import os
import sys
import time
import random

import numpy as np
import torch
from torch.utils.data import Subset, DataLoader
from torch.utils.data.distributed import DistributedSampler

from ml.data import build_onestep_loader
from ml.model import ViT

import yaml
import wandb
import matplotlib.pyplot as plt
from ml.sample import get_diffusion_samples

# Caldor fire forced to be out of the training set
FORCED_TEST_EVENT_IDS = {"CA3858612053820210815"}


def load_config(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def set_seed(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def plot_gt_vs_probs(prev_t, y, probs):
    i = random.randint(0, y.shape[0] - 1)
    prev = prev_t[i, 0].detach().cpu().numpy()
    gt = y[i, 0].detach().cpu().numpy()
    pr = probs[i, 0].detach().cpu().numpy()
    fig, axes = plt.subplots(1, 3, figsize=(9, 3))
    axes[0].imshow(prev, cmap="gray", vmin=0.0, vmax=1.0)
    axes[0].set_title("prev timestep")
    axes[0].axis("off")
    axes[1].imshow(gt, cmap="gray", vmin=0.0, vmax=1.0)
    axes[1].set_title("y (gt)")
    axes[1].axis("off")
    axes[2].imshow(pr, cmap="viridis", vmin=0.0, vmax=1.0)
    axes[2].set_title("p(y=1)")
    axes[2].axis("off")
    return fig


def plot_diffusion_visualization(prev_t, y, y_masked, probs, sample):
    i = random.randint(0, y.shape[0] - 1)
    prev = prev_t[i, 0].detach().cpu().numpy()
    gt = y[i, 0].detach().cpu().numpy()
    pr = probs[i, 0].detach().cpu().numpy()
    smp = sample[i, 0].detach().cpu().numpy()

    masked = y_masked[i, 0].detach().clone()
    masked = torch.where(masked > 1.5, torch.full_like(masked, 0.5), masked)
    masked = masked.cpu().numpy()

    fig, axes = plt.subplots(1, 5, figsize=(15, 3))
    axes[0].imshow(prev, cmap="gray", vmin=0.0, vmax=1.0)
    axes[0].set_title("prev timestep")
    axes[0].axis("off")
    axes[1].imshow(gt, cmap="gray", vmin=0.0, vmax=1.0)
    axes[1].set_title("y (gt)")
    axes[1].axis("off")
    axes[2].imshow(masked, cmap="gray", vmin=0.0, vmax=1.0)
    axes[2].set_title("masked input")
    axes[2].axis("off")
    axes[3].imshow(pr, cmap="viridis", vmin=0.0, vmax=1.0)
    axes[3].set_title("p(y=1)")
    axes[3].axis("off")
    axes[4].imshow(smp, cmap="gray", vmin=0.0, vmax=1.0)
    axes[4].set_title("sample (scratch)")
    axes[4].axis("off")
    return fig


def run_batch_ce(
    step,
    batch,
    model,
    optimizer,
    wandb_run,
    device,
    log_every,
    log_images_every,
    training=True,
    visualize=False,
    sample_steps=32,
    show_pred_vs_gt=True,
    is_master=True,
):
    x = batch["x"].to(device, non_blocking=True)
    y = batch["y"].to(device, non_blocking=True)
    mask = batch.get("mask")
    if mask is not None:
        mask = mask.to(device, non_blocking=True).float()

    logits = model(x, mask=mask) + x[:, 1].clamp(min=1e-3, max=1 - 1e-3).log().unsqueeze(1)
    probs = torch.sigmoid(logits)
    target = (y > 0.5).float()
    per_pixel = torch.nn.functional.binary_cross_entropy(
        probs, target, reduction="none"
    )
    if mask is not None:
        loss = (per_pixel * mask).sum() / mask.sum().clamp_min(1.0)
    else:
        loss = per_pixel.mean()

    split = "train" if training else "val"
    if training:
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()

        if is_master and step % log_every == 0:
            print(f"step {step} {split}_loss {loss.item():.6f}")

    if is_master and wandb_run is not None:
        wandb_run.log({f"{split}/loss": loss.item(), "step": step})
        if show_pred_vs_gt and training and log_images_every > 0 and step % log_images_every == 0:
            fig = plot_gt_vs_probs(x[:, 1:2], y, probs)
            wandb_run.log({f"{split}/pred_vs_gt": wandb.Image(fig), "step": step})
            plt.close(fig)

    if training:
        step += 1

    return step

def apply_diffusion_mask(t, y, mask):
    t = t.to(y.device).float()
    pr = 1 - t[:, None, None, None]  # probability of being masked (linear schedule)
    samples_mask = torch.rand_like(y) < pr
    valid_mask = mask > 0
    diffusion_mask = (samples_mask & valid_mask).float()
    y_masked = y.clone()
    y_masked[diffusion_mask.bool()] = 2.0
    weight = 1.0 / torch.clamp_min(1 - pr, 1e-6)
    return y_masked, diffusion_mask, weight

def run_batch_diffusion(
    step,
    batch,
    model,
    optimizer,
    wandb_run,
    device,
    log_every,
    log_images_every,
    training=True,
    visualize=False,
    sample_steps=512,
    show_pred_vs_gt=True,
    is_master=True,
):
    x = batch["x"].to(device, non_blocking=True)
    y = batch["y"].to(device, non_blocking=True)
    mask = batch["mask"].to(device, non_blocking=True).float()

    t = torch.rand(y.shape[0])
    y_masked, diffusion_mask, weight = apply_diffusion_mask(t, y, mask)

    # Explicitly provide masked-token indicator as an additional input channel.
    masked_indicator = (y_masked > 1.5).float()
    x_in = torch.cat([y_masked, masked_indicator, x], dim=1)
    logits = model(x_in, mask=mask) + x[:, 1].clamp(min=1e-3, max=1 - 1e-3).log().unsqueeze(1)
    probs = torch.sigmoid(logits)

    target = (y > 0.5).float()
    per_pixel = torch.nn.functional.binary_cross_entropy(
        probs, target, reduction="none"
    )

    # diffusion loss, only computed on pixel with diffusion_mask
    # per_pixel = weight * per_pixel * diffusion_mask
    per_pixel = per_pixel * diffusion_mask
    if mask is not None:
        loss = (per_pixel * mask).sum() / mask.sum().clamp_min(1.0)
    else:
        loss = per_pixel.mean()

    split = "train" if training else "val"
    if training:
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()

        if is_master and step % log_every == 0:
            print(f"step {step} {split}_loss {loss.item():.6f}")

    if is_master and wandb_run is not None:
        wandb_run.log({f"{split}/loss": loss.item(), "step": step})
        if visualize:
            sample = get_diffusion_samples(
                n_samples=1,
                valid_mask=mask,
                x=x,
                model=model,
                n_steps=sample_steps,
                device=device,
            )[:, 0]
            fig = plot_diffusion_visualization(x[:, 1:2], y, y_masked, probs, sample)
            wandb_run.log({f"{split}/diffusion_vis": wandb.Image(fig), "step": step})
            plt.close(fig)

    if training:
        step += 1

    return step

def _dataset_event_name(dataset, idx):
    # Frame-level dataset: index maps into dataset.sample_index -> event index.
    if hasattr(dataset, "sample_index") and hasattr(dataset, "base") and hasattr(dataset.base, "fire_events"):
        event_idx, _t = dataset.sample_index[idx]
        return dataset.base.fire_events[event_idx][0]
    # OneStepDatasetSimple wraps GeoTiffDatasetStructured in `base`.
    if hasattr(dataset, "base") and hasattr(dataset.base, "fire_events"):
        return dataset.base.fire_events[idx][0]
    # GeoTiffDatasetStructured exposes fire_events directly.
    if hasattr(dataset, "fire_events"):
        return dataset.fire_events[idx][0]
    # Fallback (slower): load sample.
    return dataset[idx]["event_name"]


def hardcoded_split(dataset, forced_test_ids=None):
    # Train/val/test split (hardcoded seed + optional forced test event IDs)
    split_seed = 123
    val_frac = 0.1
    if forced_test_ids is None:
        forced = set(FORCED_TEST_EVENT_IDS)
    elif isinstance(forced_test_ids, str):
        forced = {forced_test_ids}
    else:
        forced = set(forced_test_ids)

    n = len(dataset)
    indices = list(range(n))
    event_by_idx = {i: _dataset_event_name(dataset, i) for i in indices}

    test_idx = [i for i in indices if event_by_idx[i] in forced]
    remaining = [i for i in indices if event_by_idx[i] not in forced]

    remaining_events = sorted({event_by_idx[i] for i in remaining})

    rng = random.Random(split_seed)
    rng.shuffle(remaining_events)

    if len(remaining_events) == 0:
        return [], [], test_idx

    n_val_events = max(1, int(len(remaining_events) * val_frac))
    val_events = set(remaining_events[:n_val_events])
    val_idx = [i for i in remaining if event_by_idx[i] in val_events]
    train_idx = [i for i in remaining if event_by_idx[i] not in val_events]
    return train_idx, val_idx, test_idx


def get_training_objects(
    data_cfg,
    model_cfg,
    train_cfg,
    device,
    distributed: bool = False,
    rank: int = 0,
    world_size: int = 1,
):
    dataset, loader = build_onestep_loader(
        base_path=data_cfg["base_path"],
        batch_size=data_cfg.get("batch_size", 8),
        shuffle=bool(data_cfg.get("shuffle", True)),
        num_workers=int(data_cfg.get("num_workers", 0)),
        out_type="flattened", # default for train and val
        required_vars=data_cfg.get("required_vars"),
        target_var=data_cfg.get("target_var", "fire_spread/nfp"),
        step_hours=int(data_cfg.get("step_hours", 12)),
        horizon_hours=int(data_cfg.get("horizon_hours", 12)),
        hourly_agg=data_cfg.get("hourly_agg", "concat"),
        missing_value=float(data_cfg.get("missing_value", -1.0)),
        compute_stats=bool(data_cfg.get("compute_stats", True)),
        stats_path=data_cfg.get("stats_path"),
        stats_sample_limit=data_cfg.get("stats_sample_limit"),
        distributed=distributed,
        rank=rank,
        world_size=world_size,
    )

    train_idx, val_idx, _test_idx = hardcoded_split(dataset)

    num_workers = int(data_cfg.get("num_workers", 0))
    pin_memory = bool(data_cfg.get("pin_memory", device.type == "cuda"))
    persistent_workers = bool(data_cfg.get("persistent_workers", num_workers > 0))
    prefetch_factor = data_cfg.get("prefetch_factor", 2)
    loader_kwargs = {
        "batch_size": data_cfg.get("batch_size", 8),
        "num_workers": num_workers,
        "collate_fn": loader.collate_fn,
        "pin_memory": pin_memory,
    }
    if num_workers > 0:
        loader_kwargs["persistent_workers"] = persistent_workers
        if prefetch_factor is not None:
            loader_kwargs["prefetch_factor"] = int(prefetch_factor)

    train_ds = Subset(dataset, train_idx)
    val_ds = Subset(dataset, val_idx)
    if distributed:
        train_sampler = DistributedSampler(
            train_ds, num_replicas=world_size, rank=rank, shuffle=True
        )
        val_sampler = DistributedSampler(
            val_ds, num_replicas=world_size, rank=rank, shuffle=False
        )
        train_loader = DataLoader(
            train_ds,
            shuffle=False,
            sampler=train_sampler,
            **loader_kwargs,
        )
        val_loader = DataLoader(
            val_ds,
            shuffle=False,
            sampler=val_sampler,
            **loader_kwargs,
        )
    else:
        train_loader = DataLoader(
            train_ds,
            shuffle=True,
            **loader_kwargs,
        )
        val_loader = DataLoader(
            val_ds,
            shuffle=False,
            **loader_kwargs,
        )

    in_ch = len(dataset.fire_vars) + len(dataset.static_vars)
    hourly_mode = dataset._hourly_mode()
    if hourly_mode == "mean":
        in_ch += len(dataset.hourly_vars)
    else:
        in_ch += len(dataset.hourly_vars) * int(dataset.step)
    
    loss_type = train_cfg.get("loss_type", "diffusion")
    if loss_type == "diffusion":
        # diffusion takes noisy output + explicit masked indicator as input
        in_ch = in_ch + 2

    model = ViT(
        in_ch=in_ch,
        out_ch=int(model_cfg.get("out_ch", 2)),
        patch=int(model_cfg.get("patch", 2)),
        dim=int(model_cfg.get("dim", 128)),
        depth=int(model_cfg.get("depth", 6)),
        heads=int(model_cfg.get("heads", 4)),
        mlp_dim=int(model_cfg.get("mlp_dim", 256)),
    ).to(device)

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=float(train_cfg.get("lr", 1e-3)),
        weight_decay=float(train_cfg.get("weight_decay", 0.0)),
    )
    return train_loader, val_loader, model, optimizer
