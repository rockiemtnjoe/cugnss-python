import numpy as np

def norm(ta, opt=2):
    """
    NORM tall vector and matrix norms.
    Supported matrix X syntaxes:
        C = NORM(X) 
        C = NORM(X,2)
        C = NORM(X,1)
        C = NORM(X,np.inf)
        C = NORM(X,'fro')
    Supported vector X syntaxes:
        C = NORM(X)
        C = NORM(X,P)
        C = NORM(X,np.inf)
        C = NORM(X,-np.inf)
        C = NORM(X,'fro')
    Supported N-D array X syntax:
        C = NORM(X,'fro')
    """
    # Argument checking
    if opt is None:
        opt = 2

    # Handle string options
    if isinstance(opt, str):
        if opt.lower().startswith('inf'):
            opt = np.inf
        elif opt.lower().startswith('fro'):
            opt = 'fro'
        else:
            raise ValueError('MATLAB:norm:unknownNorm')
    elif not is_numeric_scalar(opt) and opt != 'fro':
        raise ValueError('MATLAB:norm:unknownNorm')

    ta = validate_type(ta, (np.float64, np.float32))

    if opt == 'fro':
        tc = np.linalg.norm(ta, 'fro')
    elif opt == 1:
        tc = np.nanmax(np.sum(np.abs(ta), axis=0))
    elif opt == np.inf:
        tc = inf_norm(ta)
    elif opt == 2:
        tc = two_norm(ta)
    elif opt == -np.inf:
        validate_vector(ta)
        tc = minus_inf_norm(ta)
        tc = handle_empty(tc, ta.shape, opt)
    else:
        validate_vector(ta)
        tc = np.linalg.norm(ta.ravel(), opt)
        tc = handle_empty(tc, ta.shape, opt)

    return tc

def is_numeric_scalar(opt):
    return np.isscalar(opt) and (isinstance(opt, (int, float, np.integer, np.floating, bool))) and np.isreal(opt)

def handle_empty(x, sz, opt):
    if np.prod(sz) == 0 and opt == -np.inf:
        return np.zeros_like(x)
    elif np.prod(sz) == 0 and np.isnan(opt):
        return np.full_like(x, np.nan)
    return x

def inf_norm(x):
    # Compute maximum row sum for matrices, or max(abs(x)) for vectors
    if x.ndim == 1:
        return np.max(np.abs(x))
    else:
        return np.max(np.sum(np.abs(x), axis=1))

def two_norm(x):
    # For vectors, just norm(x,2)
    if x.ndim == 1 or x.shape[1] == 1:
        return np.linalg.norm(x, 2)
    else:
        # For matrices, compute the largest singular value (matrix 2-norm)
        return np.linalg.norm(x, 2)

def minus_inf_norm(x):
    # norm(x,-inf) = min(abs(x)) where x is a vector.
    if x.size == 0:
        return np.full((1,), np.inf, dtype=x.dtype)
    else:
        return np.min(np.abs(x))

def validate_type(x, types):
    if not isinstance(x, np.ndarray):
        x = np.array(x)
    if x.dtype not in [np.dtype(t) for t in types]:
        raise TypeError('Input must be double or single precision array')
    return x

def validate_vector(x):
    if x.ndim > 1 and 1 not in x.shape:
        raise ValueError('MATLAB:norm:unknownNorm')