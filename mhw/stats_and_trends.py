import numpy as np
import pandas as pd
from scipy import stats
from typing import Tuple

def get_main_vars_matrix(
    df: pd.DataFrame, 
    size: Tuple[int, int, int] = (1, 1, 1), 
    t0: int = 1982
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Construye matrices 3D (tiempo, latitud, longitud) con las métricas 
    principales de las olas de calor marinas agregadas por año.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame que contiene las estadísticas anuales por pixel. Debe incluir 
        las columnas: 'year', 'index_lat', 'index_lon', 'duration_count', 
        'duration_mean', 'anom_max_max' y 'anom_mean_mean'.
    size : Tuple[int, int, int], optional
        Dimensiones de las matrices de salida en el orden (años, latitud, longitud).
    t0 : int, optional
        Año base que corresponde al índice 0 en la dimensión temporal. Por defecto 1982.

    Returns
    -------
    count_waves : np.ndarray
        Matriz 3D con la cantidad de eventos por año y por pixel.
    mean_duration : np.ndarray
        Matriz 3D con la duración media de los eventos por año y por pixel.
    max_intensity : np.ndarray
        Matriz 3D con la anomalía máxima absoluta por año y por pixel.
    mean_intensity : np.ndarray
        Matriz 3D con la anomalía media por año y por pixel.
    """
    nanmatrix = np.zeros(size) * np.nan

    # Arreglo del bug: los nombres ahora coinciden con los del bucle
    count_waves = nanmatrix.copy()
    mean_duration = nanmatrix.copy()
    max_intensity = nanmatrix.copy()
    mean_intensity = nanmatrix.copy()
    
    for index, row in df.iterrows():
        # Se fuerza a entero para asegurar una indexación correcta en las matrices espaciales
        _y = int(row.year) - t0
        _lat = int(row.index_lat)
        _lon = int(row.index_lon)
        
        count_waves[_y, _lat, _lon] = row.duration_count
        mean_duration[_y, _lat, _lon] = row.duration_mean
        max_intensity[_y, _lat, _lon] = row.anom_max_max
        mean_intensity[_y, _lat, _lon] = row.anom_mean_mean

    return count_waves, mean_duration, max_intensity, mean_intensity


def trend_matrix(
    df: pd.DataFrame, 
    varname: str, 
    timename: str, 
    lat_size: int, 
    lon_size: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calcula la tendencia lineal temporal para una variable específica 
    en cada pixel de la grilla espacial utilizando regresión por mínimos cuadrados.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame con los datos agregados. Debe contener columnas de coordenadas 
        ('index_lat', 'index_lon'), la variable temporal y la variable a evaluar.
    varname : str
        Nombre de la columna en el DataFrame que contiene la variable dependiente.
    timename : str
        Nombre de la columna en el DataFrame que contiene el eje temporal.
    lat_size : int
        Tamaño de la dimensión latitudinal de la grilla.
    lon_size : int
        Tamaño de la dimensión longitudinal de la grilla.

    Returns
    -------
    trend : np.ndarray
        Matriz 2D (latitud, longitud) con la pendiente (slope) de la regresión.
    trend_pvalue : np.ndarray
        Matriz 2D (latitud, longitud) con el valor p (p-value) de la tendencia.
    trend_intercept : np.ndarray
        Matriz 2D (latitud, longitud) con la ordenada al origen de la regresión.
    """
    nanmatrix = np.zeros((lat_size, lon_size)) * np.nan
    trend = nanmatrix.copy()
    trend_pvalue = nanmatrix.copy()
    trend_intercept = nanmatrix.copy()
    
    for _lat in range(lat_size):
        for _lon in range(lon_size):
            pixel = df[(df.index_lon == _lon) & (df.index_lat == _lat)].copy()
            if pixel.shape[0] > 2:
                X = pixel[timename].to_numpy()
                y = pixel[varname].to_numpy()
                
                X2 = np.array([t for t in range(int(X.min()), int(X.max()) + 1)])
                y2 = np.array([0 if X[X == t].shape[0] == 0 else y[X == t][0] for t in X2])
                
                xmean = np.nanmean(X2)
                X2 = X2 - xmean
                
                slope, intercept, r_value, p_value, std_err = stats.linregress(X2, y2)
                trend[_lat, _lon] = slope
                trend_pvalue[_lat, _lon] = p_value
                trend_intercept[_lat, _lon] = intercept
                
    return trend, trend_pvalue, trend_intercept


def get_argmax(
    index_lat: int, 
    index_lon: int, 
    index_t0: int, 
    index_tf: int, 
    var: np.ndarray
) -> Tuple[int, float, float, float, float, float]:
    """
    Extrae los estadísticos principales de un evento específico dentro 
    de una matriz 3D, acotando la búsqueda a su duración temporal y coordenada.

    Parameters
    ----------
    index_lat : int
        Índice espacial de la latitud.
    index_lon : int
        Índice espacial de la longitud.
    index_t0 : int
        Índice temporal de inicio del evento.
    index_tf : int
        Índice temporal de fin del evento.
    var : np.ndarray
        Matriz 3D (tiempo, latitud, longitud) de donde extraer los datos 
        (usualmente las anomalías calculadas previamente).

    Returns
    -------
    Tuple[int, float, float, float, float, float]
        Retorna en orden:
        - Índice temporal del pico máximo (argmax)
        - Valor máximo (vmax)
        - Valor mínimo (vmin)
        - Valor medio (vmean)
        - Desviación estándar (vstd)
        - Valor integrado o suma acumulada (vint)
    """
    _var = var[index_t0:index_tf, index_lat, index_lon].copy()
    
    arg_max = np.argmax(_var)
    vmax = float(np.nanmax(_var))
    vmin = float(np.nanmin(_var))
    vmean = float(np.nanmean(_var))
    vstd = float(np.nanstd(_var))
    vint = float(np.nansum(_var))
    
    return index_t0 + int(arg_max), vmax, vmin, vmean, vstd, vint