import random

import matplotlib
import numpy as np
import torch
from torch.utils.data import DataLoader, Subset
from torch.utils.data.distributed import DistributedSampler

import yaml

from ml.data import build_onestep_loader
from ml.model import ViT

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Caldor fire forced to be out of the training set.
FORCED_TEST_EVENT_IDS = {"CA3858612053820210815"}


def load_config(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def set_seed(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def apply_diffusion_mask(t, y, mask):
    t = t.to(y.device).float()
    pr = 1 - t[:, None, None, None]
    samples_mask = torch.rand_like(y) < pr
    valid_mask = mask > 0
    diffusion_mask = (samples_mask & valid_mask).float()
    y_masked = y.clone()
    y_masked[diffusion_mask.bool()] = 2.0
    weight = 1.0 / torch.clamp_min(1 - pr, 1e-6)
    return y_masked, diffusion_mask, weight


def _farea_channel_index(batch):
    fire_vars = batch.get("fire_vars", [])
    if "fire_spread/farea" in fire_vars:
        return fire_vars.index("fire_spread/farea")
    return 1 if len(fire_vars) > 1 else None


def _add_persistence_prior(logits, x, batch):
    farea_idx = _farea_channel_index(batch)
    if farea_idx is None or farea_idx >= x.shape[1]:
        return logits
    prev_fire_area = x[:, farea_idx].clamp(min=1e-3, max=1 - 1e-3)
    return logits + prev_fire_area.log().unsqueeze(1)


def _diffusion_logits(model, x, y_masked, mask, batch):
    masked_indicator = (y_masked > 1.5).float()
    x_in = torch.cat([y_masked, masked_indicator, x], dim=1)
    logits = model(x_in, mask=mask)
    return _add_persistence_prior(logits, x, batch)


def all_masked_probability_map(batch, model, device):
    x = batch["x"].to(device, non_blocking=True)
    y = batch["y"].to(device, non_blocking=True)
    mask = batch["mask"].to(device, non_blocking=True).float()

    y_masked = torch.full_like(y, 2.0)
    y_masked = torch.where(mask > 0, y_masked, torch.zeros_like(y_masked))
    logits = _diffusion_logits(model, x, y_masked, mask, batch)
    return torch.sigmoid(logits), y_masked


def plot_probability_map(batch, probs, y_masked):
    i = random.randint(0, batch["y"].shape[0] - 1)
    x = batch["x"]
    y = batch["y"]
    farea_idx = _farea_channel_index(batch)
    if farea_idx is not None and farea_idx < x.shape[1]:
        prev = x[i, farea_idx].detach().cpu().numpy()
    else:
        prev = np.zeros_like(y[i, 0].detach().cpu().numpy())
    gt = y[i, 0].detach().cpu().numpy()
    masked = y_masked[i, 0].detach().cpu().numpy()
    pr = probs[i, 0].detach().cpu().numpy()

    fig, axes = plt.subplots(1, 4, figsize=(12, 3))
    axes[0].imshow(prev, cmap="gray", vmin=0.0, vmax=1.0)
    axes[0].set_title("prev farea")
    axes[0].axis("off")
    axes[1].imshow(gt, cmap="gray", vmin=0.0, vmax=1.0)
    axes[1].set_title("target")
    axes[1].axis("off")
    axes[2].imshow(masked, cmap="gray", vmin=0.0, vmax=2.0)
    axes[2].set_title("all masked")
    axes[2].axis("off")
    axes[3].imshow(pr, cmap="viridis", vmin=0.0, vmax=1.0)
    axes[3].set_title("p(farea=1)")
    axes[3].axis("off")
    fig.tight_layout()
    return fig


def run_batch_diffusion(
    step,
    batch,
    model,
    optimizer,
    device,
    training=True,
    visualize=False,
    grad_accum_steps=1,
    accum_step_idx=0,
    max_grad_norm=0.0,
):
    x = batch["x"].to(device, non_blocking=True)
    y = batch["y"].to(device, non_blocking=True)
    mask = batch["mask"].to(device, non_blocking=True).float()

    t = torch.rand(y.shape[0], device=device)
    y_masked, diffusion_mask, _weight = apply_diffusion_mask(t, y, mask)
    logits = _diffusion_logits(model, x, y_masked, mask, batch)

    target = (y > 0.5).float()
    per_pixel = torch.nn.functional.binary_cross_entropy_with_logits(
        logits, target, reduction="none"
    )
    per_pixel = per_pixel * diffusion_mask
    loss = (per_pixel * mask).sum() / mask.sum().clamp_min(1.0)

    grad_accum_steps = max(1, int(grad_accum_steps))
    accum_step_idx = int(accum_step_idx)
    accum_start = (accum_step_idx % grad_accum_steps) == 0
    accum_end = ((accum_step_idx + 1) % grad_accum_steps) == 0
    optimizer_step = bool(training and accum_end)

    if training:
        if accum_start:
            optimizer.zero_grad(set_to_none=True)
        (loss / grad_accum_steps).backward()
        if accum_end:
            grad_ok = True
            if max_grad_norm and max_grad_norm > 0:
                grad_norm = torch.nn.utils.clip_grad_norm_(
                    model.parameters(), max_grad_norm, error_if_nonfinite=False
                )
                grad_ok = bool(torch.isfinite(grad_norm))
            if grad_ok:
                optimizer.step()
            else:
                optimizer.zero_grad(set_to_none=True)
                optimizer_step = False

    fig = None
    if visualize:
        probs, all_masked = all_masked_probability_map(batch, model, device)
        fig = plot_probability_map(batch, probs.detach().cpu(), all_masked.detach().cpu())

    if optimizer_step:
        step += 1

    return step, {
        "loss": float(loss.detach().cpu().item()),
        "optimizer_step": optimizer_step,
        "figure": fig,
    }


def _dataset_event_name(dataset, idx):
    if hasattr(dataset, "sample_index") and hasattr(dataset, "base") and hasattr(dataset.base, "fire_events"):
        event_idx, _t = dataset.sample_index[idx]
        return dataset.base.fire_events[event_idx][0]
    if hasattr(dataset, "base") and hasattr(dataset.base, "fire_events"):
        return dataset.base.fire_events[idx][0]
    if hasattr(dataset, "fire_events"):
        return dataset.fire_events[idx][0]
    return dataset[idx]["event_name"]


def hardcoded_split(dataset, forced_test_ids=None):
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

    if len(remaining_events) == 0:
        return [], [], test_idx
    if len(remaining_events) == 1:
        return remaining, remaining, test_idx

    rng = random.Random(split_seed)
    rng.shuffle(remaining_events)
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
        out_type="flattened",
        required_vars=data_cfg.get("required_vars"),
        target_var=data_cfg.get("target_var", "fire_spread/farea"),
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
    if len(train_idx) == 0:
        raise RuntimeError("Training split is empty.")

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
    in_ch += 2

    model = ViT(
        in_ch=in_ch,
        out_ch=int(model_cfg.get("out_ch", 1)),
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
