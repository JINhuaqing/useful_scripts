def lin_R_fn(x, y):
    """
    For both torch and np
    Calculate the linear correlation coefficient (Lin's R) between x and y.
    
    Args:
    x: torch.Tensor, shape (batch_size, num_features)
    y: torch.Tensor, shape (batch_size, num_features)
    
    Returns:
    ccc: torch.Tensor, shape (batch_size,)
    """
    assert x.shape == y.shape, "x and y should have the same shape"
    x_bar = x.mean(axis=-1, keepdims=True)
    y_bar = y.mean(axis=-1, keepdims=True)
    num = ((x-x_bar)*(y-y_bar)).sum(axis=-1);
    den = (x**2).sum(axis=-1) + (y**2).sum(axis=-1) - (2 * x.shape[-1] * x_bar * y_bar).squeeze()
    ccc = num/den;
    return ccc
