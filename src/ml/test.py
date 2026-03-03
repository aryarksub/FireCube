import torch
import torch.nn as nn

from ml.train import apply_diffusion_mask, run_batch_diffusion
from ml.model import ViT

class DummyModel(nn.Module):
    def __init__(self, in_ch):
        super().__init__()
        self.net = nn.Conv2d(in_ch, 1, kernel_size=1)

    def forward(self, x, mask=None):
        return self.net(x)


def _make_batch(batch_size=2, channels=3, h=8, w=8):
    x = torch.rand(batch_size, channels, h, w)
    y = (torch.rand(batch_size, 1, h, w) > 0.5).float()
    mask = torch.ones(batch_size, 1, h, w)
    # Hide a corner to validate masking behavior.
    mask[:, :, :2, :2] = 0
    return {"x": x, "y": y, "mask": mask}


def test_apply_diffusion_mask_shapes():
    batch = _make_batch(batch_size=4, h=6, w=5)
    t = torch.rand(batch["y"].shape[0])
    y_masked, diffusion_mask, weight = apply_diffusion_mask(t, batch["y"], batch["mask"])

    assert y_masked.shape == batch["y"].shape
    assert diffusion_mask.shape == batch["y"].shape
    assert weight.shape == (batch["y"].shape[0], 1, 1, 1)
    assert torch.all(weight > 0)


def test_apply_diffusion_mask_respects_visibility_mask():
    batch = _make_batch(batch_size=3, h=10, w=10)
    t = torch.ones(batch["y"].shape[0]) * 0.5
    _, diffusion_mask, _ = apply_diffusion_mask(t, batch["y"], batch["mask"])

    # Diffusion mask must be zero where visibility mask is zero.
    invalid_region = batch["mask"] == 0
    assert torch.count_nonzero(diffusion_mask[invalid_region]) == 0


def test_run_batch_diffusion_increments_step_and_updates_weights():
    torch.manual_seed(0)
    batch = _make_batch(batch_size=2, channels=4, h=8, w=8)
    model = DummyModel(in_ch=batch["x"].shape[1] + 1)  # +1 for y_masked channel
    optimizer = torch.optim.SGD(model.parameters(), lr=1e-2)

    before = model.net.weight.detach().clone()
    step_out = run_batch_diffusion(
        step=7,
        batch=batch,
        model=model,
        optimizer=optimizer,
        wandb_run=None,
        device=torch.device("cpu"),
        log_every=1000,
        log_images_every=0,
        training=True,
    )
    after = model.net.weight.detach().clone()

    assert step_out == 8
    assert torch.any(torch.ne(before, after))


def test_run_batch_diffusion_eval_does_not_change_weights_or_step():
    torch.manual_seed(0)
    batch = _make_batch(batch_size=2, channels=4, h=8, w=8)
    model = DummyModel(in_ch=batch["x"].shape[1] + 1)
    optimizer = torch.optim.SGD(model.parameters(), lr=1e-2)

    before = model.net.weight.detach().clone()
    step_out = run_batch_diffusion(
        step=12,
        batch=batch,
        model=model,
        optimizer=optimizer,
        wandb_run=None,
        device=torch.device("cpu"),
        log_every=1000,
        log_images_every=0,
        training=False,
    )
    after = model.net.weight.detach().clone()

    assert step_out == 12
    assert torch.allclose(before, after)

def test_vit_diffusion():
    torch.manual_seed(0)
    batch = _make_batch(batch_size=2, channels=4, h=8, w=8)
    model = ViT(
        in_ch=batch["x"].shape[1] + 1,
        out_ch=1,
        patch=2,
        dim=128,
        depth=6,
        heads=4,
        mlp_dim=256,
    )

    optimizer = torch.optim.SGD(model.parameters(), lr=1e-2)

    step_out = run_batch_diffusion(
        step=12,
        batch=batch,
        model=model,
        optimizer=optimizer,
        wandb_run=None,
        device=torch.device("cpu"),
        log_every=1000,
        log_images_every=0,
        training=False,
    )
    assert step_out == 12

def run_all():
    tests = [
        test_apply_diffusion_mask_shapes,
        test_apply_diffusion_mask_respects_visibility_mask,
        test_run_batch_diffusion_increments_step_and_updates_weights,
        test_run_batch_diffusion_eval_does_not_change_weights_or_step,
        test_vit_diffusion
    ]
    for test_fn in tests:
        test_fn()
        print(f"PASS: {test_fn.__name__}")
    print("All diffusion tests passed.")


if __name__ == "__main__":
    run_all()
