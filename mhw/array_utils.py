import numpy as np
from typing import Tuple

def maxmin(r: np.ndarray) -> Tuple[float, float]:
    """
    Devuelve el valor mínimo y máximo de un array ignorando los NaNs.
    """
    return float(np.nanmin(r)), float(np.nanmax(r))


def limited_matrix(M: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    """
    Limita los valores de una matriz a un rango específico [vmin, vmax].
    """
    M = M.copy()
    M[M > vmax] = vmax
    M[M < vmin] = vmin
    return M


def false_to_nan(M: np.ndarray) -> np.ndarray:
    """
    Convierte los valores 0 (o False) de una matriz numérica a np.nan.
    """
    M_float = M.astype(np.float64)
    M_float[M_float == 0] = np.nan
    return M_float


def redim_matrix(matrix_2d: np.ndarray, time_size: int) -> np.ndarray:
    """
    Redimensiona una matriz espacial 2D (latitud, longitud) a una matriz 3D 
    (tiempo, latitud, longitud) repitiendo los valores en el nuevo eje temporal.

    Parameters
    ----------
    matrix_2d : np.ndarray
        Matriz espacial 2D (ej. una climatología de un día particular).
    time_size : int
        Cantidad de pasos temporales para los cuales repetir la matriz.

    Returns
    -------
    np.ndarray
        Matriz 3D con forma (time_size, latitud, longitud).
    """
    # Agrega una dimensión al principio y repite la matriz time_size veces
    return np.repeat(matrix_2d[np.newaxis, :, :], time_size, axis=0)


def matrices_son_iguales(A, B):
    # 1. Chequear primero si las dimensiones son iguales
    if A.shape != B.shape:
        return False
    
    # 2. Crear máscaras de dónde están los NaNs en ambas
    nan_A = np.isnan(A)
    nan_B = np.isnan(B)
    
    # 3. Si las posiciones de los NaNs no coinciden, son distintas
    if not np.array_equal(nan_A, nan_B):
        return False
    
    # 4. Comparar el resto de los valores (ignorando los NaNs)
    # Usamos isclose por si hay diferencias mínimas de punto flotante
    return np.all(np.isclose(A[~nan_A], B[~nan_B]))
