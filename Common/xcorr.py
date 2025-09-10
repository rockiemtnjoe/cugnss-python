import numpy as np
from numpy.fft import fft, ifft

def find_transform_length(m):
    m = 2 * m
    while True:
        r = m
        for p in [2, 3, 5, 7]:
            while r > 1 and r % p == 0:
                r //= p
        if r == 1:
            return m
        m += 1

def scale_output(scale, c, x, y=None):
    if scale == 'none':
        return c
    m = len(x)
    if scale == 'biased':
        return c / m
    elif scale == 'unbiased':
        L = (len(c) - 1) // 2
        scale_factors = np.maximum(1, m - np.abs(np.arange(-L, L + 1)))
        return c / scale_factors
    elif scale in ['normalized', 'coeff']:
        if y is not None:
            cxx0 = np.sum(np.abs(x) ** 2)
            cyy0 = np.sum(np.abs(y) ** 2)
            return c / np.sqrt(cxx0 * cyy0)
        else:
            mid = len(c) // 2
            return c / c[mid]
    return c

def crosscorr(x, y, maxlag):
    nx, ny = len(x), len(y)
    m = max(nx, ny)
    m2 = find_transform_length(m)
    X = fft(x, m2)
    Y = fft(y, m2)
    c = ifft(X * np.conj(Y))
    return np.concatenate([c[-maxlag:], c[:maxlag + 1]]).real

def autocorr(x, maxlag):
    m, n = x.shape
    m2 = find_transform_length(m)
    X = fft(x, m2, axis=0)
    c = ifft(np.abs(X)**2, axis=0)
    return np.concatenate([c[-maxlag:], c[:maxlag + 1]]).real

def xcorr(x, y=None, maxlag=None, scale='none'):
    x = np.atleast_2d(x)
    if x.shape[0] == 1:
        x = x.T
    m = x.shape[0]

    if y is not None:
        y = np.atleast_1d(y)
        if y.ndim > 1:
            y = y.squeeze()
        if maxlag is None:
            maxlag = max(len(x), len(y)) - 1
        c = crosscorr(x[:, 0], y, maxlag)
        c = scale_output(scale, c, x[:, 0], y)
    else:
        if maxlag is None:
            maxlag = m - 1
        c = autocorr(x, maxlag)
        c = scale_output(scale, c, x[:, 0])

    lags = np.arange(-maxlag, maxlag + 1)
    return c, lags
