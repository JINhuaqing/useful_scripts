import numpy as np
def mag2db(y):
    """Convert magnitude response to decibels for a simple array.

    Args:
        y (numpy array): Power spectrum, raw magnitude response.

    Returns:
        dby (numpy array): Power spectrum in dB

    """
    # if y is in raw mag, 20, if squared, 10
    dby = 20 * np.log10(y)
    return dby
