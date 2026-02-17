import torch

def get_revealed_positions(n, diff_mask, valid_mask, strategy = "random"):
    Br, _, H, W = diff_mask.shape
    eligible = diff_mask & valid_mask
    if strategy == "random":
        u = torch.rand_like(diff_mask.float())
        u = torch.where(eligible, u, torch.full_like(u, -1e9))
    else:
        # TODO: confidence-based unmasking strategy
        raise NotImplementedError()
    flat_score = u.view(Br, -1) # [Br, HW]
    flat_eligible = eligible.view(Br, -1) # [Br, HW]

    order = torch.argsort(flat_score, dim=1, descending=True)
    rank = torch.argsort(order, dim=1)

    reveal = rank < n[:, None]
    reveal = reveal & flat_eligible

    return reveal.view(Br, 1, H, W)

    

def get_diffusion_samples(n_samples, valid_mask, x, model, n_steps, device):
    B, C, H, W = x.shape

    # Repeat inputs per sample
    x = x.to(device)
    valid_mask = valid_mask.to(device).bool()
    x = torch.repeat_interleave(x, n_samples, dim=0) # [B*n_samples, C, H, W]
    valid_mask = torch.repeat_interleave(valid_mask, n_samples, dim=0)  # [B*n_samples, 1, H, W]
    Br = x.shape[0]

    # Time grid for denoising
    t_steps = torch.linspace(0, 1, steps=n_steps, device=device)
    
    # masked state: 2 means "masked token"
    diff_mask = valid_mask.clone() # True where still masked
    y_t = torch.full((Br, 1, H, W), 2.0, device=device)
    y_t = torch.where(valid_mask, y_t, torch.zeros_like(y_t)) # invalid locations forced to 0

    # number of valid pixels per sample
    n_pixels = valid_mask.view(Br, -1).sum(dim=1) # [Br]

    for t in t_steps:
        masked_indicator = diff_mask.float()
        logits = model(torch.cat([y_t, masked_indicator, x], dim=1), mask=valid_mask.float())
        logits = logits + x[:, 1].clamp(min=1e-3, max=1 - 1e-3).log().unsqueeze(1)
        probs = torch.sigmoid(logits)
        sampled = torch.bernoulli(probs)
        
        # Unmask new pixels
        n_target = torch.clamp_min((t * n_pixels.float()).long(), 0)
        n_revealed = ((~diff_mask) & valid_mask).view(Br, -1).sum(dim=1)
        n_eligible = diff_mask.view(Br, -1).sum(dim=1)
        n_new = torch.clamp_min(n_target - n_revealed, 0)
        n_new = torch.minimum(n_new, n_eligible)
        reveal_new = get_revealed_positions(n_new, diff_mask, valid_mask)
        
        # Preserve already revealed pixels; only fill newly revealed positions.
        y_t = torch.where(reveal_new, sampled, y_t)
        y_t = torch.where(valid_mask, y_t, torch.zeros_like(y_t))
        diff_mask = diff_mask & (~reveal_new)

    # reshape back to [B, n_samples, 1, H, W]
    y_t = y_t.view(B, n_samples, 1, H, W)
    return y_t


        

        
