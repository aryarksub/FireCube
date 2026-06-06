# --------------------------------------------------------
# Source: https://github.com/whlzy/FiT/blob/master/fit/model/rope.py
#
# Based on the following repositories:
# https://github.com/lucidrains/rotary-embedding-torch
# https://github.com/jquesnelle/yarn/blob/HEAD/scaled_rope
# --------------------------------------------------------

import torch
from torch import nn


def rotate_half(x):
    x = x.reshape(*x.shape[:-1], -1, 2)
    x1, x2 = x.unbind(dim=-1)
    return torch.stack((-x2, x1), dim=-1).flatten(-2)


class VisionRotaryEmbedding(nn.Module):
    def __init__(self, head_dim: int, theta: int = 10000, max_cached_len: int = 2048):
        super().__init__()
        dim = head_dim // 2
        if dim % 2 != 0:
            raise ValueError("head_dim must be divisible by 4 for 2D RoPE.")

        freqs = 1.0 / (theta ** (torch.arange(0, dim, 2).float() / dim))
        positions = torch.arange(max_cached_len).float()
        freqs_cached = torch.einsum("n,f->nf", positions, freqs)
        freqs_cached = torch.repeat_interleave(freqs_cached, repeats=2, dim=-1)

        self.register_buffer("freqs_h_cached", freqs_cached, persistent=False)
        self.register_buffer("freqs_w_cached", freqs_cached.clone(), persistent=False)

    def get_cached_2d_rope_from_grid(self, grid: torch.Tensor):
        """
        Args:
            grid: Tensor with shape (B, 2, N). The first coordinate is width,
                the second is height.
        """
        freqs_w = self.freqs_w_cached[grid[:, 0]]
        freqs_h = self.freqs_h_cached[grid[:, 1]]
        freqs = torch.cat([freqs_h, freqs_w], dim=-1)
        return freqs.cos(), freqs.sin()

    def forward(self, x, grid):
        """
        Args:
            x: Tensor with shape (B, heads, N, head_dim).
            grid: Tensor with shape (B, 2, N).
        """
        freqs_cos, freqs_sin = self.get_cached_2d_rope_from_grid(grid)
        freqs_cos = freqs_cos.unsqueeze(1)
        freqs_sin = freqs_sin.unsqueeze(1)
        return x * freqs_cos + rotate_half(x) * freqs_sin
