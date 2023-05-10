import numpy as np
import torch
def reg_R_fn(x, y):
    """Calculate pearons'r in batch, for both numpy and torch
    Args:
    x: torch.Tensor, shape (batch_size, num_features)
    y: torch.Tensor, shape (batch_size, num_features)
    Returns:
    corrs: torch.Tensor, shape (batch_size,)
    """
    assert x.shape == y.shape, "x and y should have the same shape"
    x_mean = x.mean(axis=-1, keepdims=True)
    y_mean = y.mean(axis=-1, keepdims=True)
    num = ((x- x_mean)*(y-y_mean)).sum(axis=-1)
    if isinstance(x, np.ndarray):
        den = np.sqrt(((x- x_mean)**2).sum(axis=-1)*((y-y_mean)**2).sum(axis=-1))
    else:
        den = torch.sqrt(((x- x_mean)**2).sum(axis=-1)*((y-y_mean)**2).sum(axis=-1))
    corrs = num/den
    return corrs

