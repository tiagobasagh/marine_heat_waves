import numpy as np

def maxmin(r):
    return np.nanmin(r), np.nanmax(r)

def limited_matrix(M, vmin, vmax):
    M[M > vmax] = vmax
    M[M < vmin] = vmin
    return M

def false_to_nan(M):
    M = M.astype(np.float64)
    M[M==0] = np.nan
    return M

def get_argmax(index_lat, index_lon, index_t0, index_tf, var): 
    """
    """
    _var = var[index_t0:index_tf, index_lat, index_lon].copy()
    arg_max = np.argmax(_var)
    vmax = np.nanmax(_var)
    vmin = np.nanmin(_var)
    vmean = np.nanmean(_var)
    vstd = np.nanstd(_var)
    vint = np.nansum(_var)
    return index_t0 + arg_max, vmax, vmin, vmean, vstd, vint

def get_season(n):
    if n in (1, 2, 3):
        s = "summer"
    elif n in (4, 5, 6):
        s = "autumn"
    elif n in (7, 8, 9):
        s = "winter"
    elif n in (10, 11, 12):
        s = "spring"
    else:
        s = "error"
    return s