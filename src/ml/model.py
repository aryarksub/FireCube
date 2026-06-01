import torch
import torch.nn as nn
import torch.nn.functional as F
from .rope_2d import VisionRotaryEmbedding


class ViT(nn.Module):
    def __init__(
        self,
        in_ch,
        out_ch=1,
        patch=2,
        dim=64,
        depth=2,
        heads=4,
        mlp_dim=128,
        rope_max_cached_len=2048,
    ):
        super().__init__()
        self.patch = patch
        self.out_ch = out_ch
        self.heads = heads

        self.proj = nn.Linear(in_ch * patch * patch, dim)
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=dim, nhead=heads, dim_feedforward=mlp_dim, batch_first=True
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=depth)
        self.head = nn.Linear(dim, out_ch * patch * patch)

        if dim % heads != 0:
            raise ValueError("dim must be divisible by heads for RoPE.")
        self.rope = VisionRotaryEmbedding(
            head_dim=dim // heads,
            theta=10000,
            max_cached_len=rope_max_cached_len,
        )

    def _mask_to_tokens(self, mask, H, W):
        B = mask.shape[0]
        ph = pw = self.patch
        m = mask.reshape(B, 1, H // ph, ph, W // pw, pw)
        m = m.permute(0, 2, 4, 1, 3, 5)
        m = m.reshape(B, (H // ph) * (W // pw), ph * pw)
        return m.any(dim=-1)

    def _pad_to_patch(self, x, mask):
        B, C, H, W = x.shape
        ph = pw = self.patch
        pad_h = (ph - (H % ph)) % ph
        pad_w = (pw - (W % pw)) % pw
        if pad_h == 0 and pad_w == 0:
            return x, mask, H, W
        x = F.pad(x, (0, pad_w, 0, pad_h))
        if mask is not None:
            mask = F.pad(mask, (0, pad_w, 0, pad_h), value=0)
        return x, mask, H, W

    def forward(self, x, mask=None):
        x, mask, H0, W0 = self._pad_to_patch(x, mask)
        B, C, H, W = x.shape

        ph = pw = self.patch
        x = x.reshape(B, C, H // ph, ph, W // pw, pw)
        x = x.permute(0, 2, 4, 1, 3, 5)
        x = x.reshape(B, (H // ph) * (W // pw), C * ph * pw)

        x = self.proj(x)
        Ht, Wt = H // ph, W // pw
        h = torch.arange(Ht, device=x.device)
        w = torch.arange(Wt, device=x.device)
        hh, ww = torch.meshgrid(h, w, indexing="ij")
        grid = torch.stack([ww.reshape(-1), hh.reshape(-1)], dim=0)  # (2, N)
        grid = grid.unsqueeze(0).expand(B, -1, -1)  # (B, 2, N)
        x_heads = x.reshape(B, -1, self.heads, x.shape[-1] // self.heads)
        x_heads = x_heads.permute(0, 2, 1, 3)  # (B, heads, N, head_dim)
        x_heads = self.rope(x_heads, grid)
        x = x_heads.permute(0, 2, 1, 3).reshape(B, -1, x.shape[-1])

        key_padding_mask = None
        if mask is not None:
            token_valid = self._mask_to_tokens(mask, H, W)
            key_padding_mask = ~token_valid

        x = self.encoder(x, src_key_padding_mask=key_padding_mask)
        x = self.head(x)

        x = x.reshape(B, H // ph, W // pw, self.out_ch, ph, pw)
        x = x.permute(0, 3, 1, 4, 2, 5)
        x = x.reshape(B, self.out_ch, H, W)

        if H != H0 or W != W0:
            x = x[:, :, :H0, :W0]
        return x
