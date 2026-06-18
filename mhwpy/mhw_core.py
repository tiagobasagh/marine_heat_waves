import numpy as np
import pandas as pd
import xarray as xr
from typing import Tuple
from tqdm import tqdm

import logging

# Configuración básica
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.StreamHandler()
    ]
)

logger = logging.getLogger(__name__)



def get_climatology_and_anomalies(
    sst_da: xr.DataArray, 
    clim_year_init: int = 1991, 
    clim_year_end: int = 2020, 
    window_time: int = 5,
    percentil: float = .9,
    coldwaves=False,
) -> xr.Dataset:
    """
    Calcula la climatología media diaria, el percentil 90 y las anomalías 
    de SST utilizando Xarray para vectorización espacial y temporal.

    Parameters
    ----------
    sst_da : xr.DataArray
        DataArray de SST con dimensión temporal 'time' (y típicamente 'lat', 'lon').
    clim_year_init : int, optional
        Año de inicio para el período de climatología base. Por defecto 1991.
    clim_year_end : int, optional
        Año de fin para el período de climatología base. Por defecto 2020.
    window_time : int, optional
        Ancho de la ventana temporal en días (hacia adelante y hacia atrás). 
        Por defecto 5 (crea una ventana total de 11 días).
    percentil: float, optional
        Percentil a calcularse. Por defecto es 0.9,
    coldwaves: boolean, optional
        Tipo de ola a calcular. Por defecto, olas de calor
    Returns
    -------
    xr.Dataset
        Dataset que contiene 'climatology', 'sst_threshold', 'anomalias', e 'posible_mhw_mask'.
        Mantiene las coordenadas y dimensiones originales.
    """
    # Removemos dimensiones menores (depth/profunidad)
    sst_da = sst_da.squeeze()
    #  Construir la ventana móvil SOBRE TODA LA SERIE (para no perder el efecto borde)
    window_size = window_time * 2 + 1
    info_coldwaves = "Olas de Frio" if coldwaves else "Olas de Calor"
    q_val = (1 - percentil) if coldwaves else percentil
    
    info = f""" Parametros del proceso:
        Climatologia: desde {clim_year_init} a {clim_year_end}
        Eventos a calcular: {info_coldwaves}
        Percentil: {q_val}
        Ventana: {window_size} dias
    """
    logger.info(info )

    # Asignando ventana movil para calculo climatologia
    sst_rolled = sst_da.rolling(time=window_size, center=True).construct("window_dim")
    
    # Eliminando valores fuera del periodo de la climatologia esperada
    base_rolled = sst_rolled.sel(time=slice(str(clim_year_init), str(clim_year_end)))
    
    # Agrupando por mes-dia
    month_day = base_rolled.time.dt.strftime('%m-%d')
    climatology = base_rolled.groupby(month_day).mean(dim=["time", "window_dim"], skipna=True)
    
    
    
    logger.info(f"Calculando percentil {q_val} de SST. El proceso puede durar algunos minutos.")
    threshold = base_rolled.groupby(month_day).quantile(
            q_val, 
            dim=["time", "window_dim"], 
            skipna=True
        ).squeeze(drop=True)
    # Limpieza des dimensiones auxiales de thrshold antes que se hereden
    threshold = threshold.drop_vars('quantile', errors='ignore')
    logger.info(f"Calculando anomalias y posibles mhw")    
    sst_grouped = sst_da.groupby(sst_da.time.dt.strftime('%m-%d'))
    anomalies = sst_grouped - climatology
    posible_mhw_mask = sst_grouped < threshold if coldwaves else sst_grouped > threshold

    # Limpieza de variables y dimensiones auxiliares
    anomalies = anomalies.drop_vars("dayofyear", errors="ignore")
    posible_mhw_mask = posible_mhw_mask.drop_vars("dayofyear", errors="ignore")
   
    logger.info(f"Creando dataset")    
    ds = xr.Dataset({
        "climatology": climatology,
        "sst_threshold": threshold,
        "anomalies": anomalies,
        "posible_mhw_mask": posible_mhw_mask
    })

    ds.attrs["percentil_base"] = q_val
    ds.attrs["is_coldwaves"] = int(coldwaves)
    ds.attrs["window_time"] = window_time
    logger.info(f"Dataset creado")    
    return ds


def _anomalies_brief(
        wave_anom: xr.DataArray,
        coldwave: bool=False
    ) -> dict:
    """Calcula las métricas de intensidad y tasas dinámicas (crecimiento y

    decaimiento) para el ciclo de vida de una Ola de Calor Marina (MHW).

    Parámetros
    ----------
    wave_anom : np.ndarray o pd.Series
        Vector con las anomalías de temperatura de la superficie del mar (SSTa)
        durante la duración exacta del evento.
    coldwae: Boolean
        Booleano que indica si las metricas se calculan para una ola de frio
        o una olar de calor. Por defecto: False (Ola de calor)

    Devuelve
    --------
    dict
        Métricas clave de la MHW (máxima, media, acumulada, tasas, etc.).
    """
    #  Métricas Clásicas
    if coldwave:
        max_anom = min(wave_anom)  # Intensidad máxima (pico del evento)
        min_anom = max(wave_anom)  # Anomalía mínima registrada en el período
        idx_max = np.argmin(wave_anom) # Día en el que se alcanzó el pico
    else:
        max_anom = max(wave_anom)  # Intensidad máxima (pico del evento)
        min_anom = min(wave_anom)  # Anomalía mínima registrada en el período
        idx_max = np.argmax(wave_anom) # Día en el que se alcanzó el pico
    
    mean_anom = wave_anom.mean()  # Intensidad media de la ola
    std_anom = wave_anom.std()  # Variabilidad interna de la anomalía

    # Intensidad acumulada (°C-días): Medida del estrés térmico total
    acc_anom = wave_anom.sum()

    # Análisis Dinámico (Evolución temporal)
    n_days = wave_anom.shape[0]  # Duración total de la ola en días

    # Tasa de crecimiento: velocidad de calentamiento inicial hasta el pico (°C/día)
    growth_anom = (
        (max_anom - wave_anom[0]) / idx_max if idx_max > 0 else np.nan
    )

    # Tasa de decaimiento: velocidad de enfriamiento desde el pico hasta el final (°C/día)
    days_to_decline = n_days - 1 - idx_max
    decline_anom = (
        (wave_anom[-1] - max_anom) / days_to_decline
        if days_to_decline > 0
        else np.nan
    )

    return dict(
        max_anom=max_anom,
        min_anom=min_anom,
        mean_anom=mean_anom,
        std_anom=std_anom,
        acc_anom=acc_anom,
        growth_anom=growth_anom,
        decline_anom=decline_anom,
    )


def waves_to_dataframe(
    waves_da: xr.DataArray,
    anomalies: xr.DataArray = None,
    lat_name: str = "latitude",
    lon_name: str = "longitude",
    time_name: str = "time",
    coldwaves: bool=False,
) -> pd.DataFrame:
    """Transforma la matriz 3D de Olas de Calor Marinas (MHW) en un DataFrame compilado

    utilizando iteración secuencial clásica por píxel.

    Parámetros
    ----------
    waves_da : xr.DataArray
        Matriz 3D (time, lat, lon) con 1 en días con presencia de MHW y 0/NaN
        en el resto.
    anomalies : xr.DataArray, opcional
        Matriz 3D con las anomalías de temperatura de la superficie del mar
        (SSTa).
    lat_name : str, opcional
        Nombre de la dimensión de latitud. Por defecto "latitude".
    lon_name : str, opcional
        Nombre de la dimensión de longitud. Por defecto "longitude".
    time_name : str, opcional
        Nombre de la dimensión de tiempo. Por defecto "time".
    coldwaves: bool, opcional
        Indica si son olas de calor/frio. Por defecto False(olas de calor)
    Devuelve
    --------
    pd.DataFrame
        Catálogo/Inventario indexado donde cada fila representa un evento único.
    """
    logger.info(f"Preparando matrices para crear DataFrame") 
    # Preparación de datos base 
    data = waves_da.fillna(0).values.astype(np.int8)

    # Controlar la existencia de anomalías 
    tiene_anomalias = anomalies is not None
    if tiene_anomalias:
        logger.info(f"Se agregaran estadisticos de anomalias al inventario")    
        # Se elimina .values dentro del bucle porque 'anom' ya será un array NumPy
        anom = np.squeeze(anomalies.fillna(0).values)
        anom_stats = []

    # Extraer vectores de coordenadas originales
    original_lats = waves_da[lat_name].values
    original_lons = waves_da[lon_name].values
    original_time = waves_da[time_name].values
    size_time = original_time.shape[0]

    # Listas contenedoras para armar el DataFrame estructurado
    lats, lons, t0s, tfs, duration = [], [], [], [], []

    #  Bucle espacial anidado (Píxel por píxel)
    logger.info(f"Recorriendo pixel por pixel (puede demorar unos instantes)")    
    for ilat, lat in tqdm(enumerate(original_lats), desc="Latitud"):
        for ilon, lon in enumerate(original_lons):
            serie = data[:, ilat, ilon]
            t = 0

            # Escaneo temporal lineal dentro del píxel activo
            while t < size_time:
                if serie[t]:  # Si se detecta el inicio de un evento (1)
                    t0 = t
                    while t < size_time and serie[t]:
                        t += 1
                    tf = t - 1

                    # Almacenar propiedades básicas del evento detectado
                    lats.append(lat)
                    lons.append(lon)
                    t0s.append(original_time[t0])
                    tfs.append(original_time[tf])
                    duration.append((tf - t0) + 1)  # Duración inclusiva

                    # Extraer métricas estadísticas si se pasaron anomalías
                    if tiene_anomalias:
                        # Extraemos el sub-vector temporal de anomalías para este evento
                        segmento_anom =  anom[t0 : tf + 1, ilat, ilon]
                        anom_stats.append(_anomalies_brief(segmento_anom, 
                                                           coldwave=coldwaves
                                                           ))
                else:
                    t += 1

    logger.info(f"Creando DataFrame")    
    if len(lats) == 0:
        # Resguardo si no se detectó ningún evento en toda la grilla
        cols = ["lat", "lon", "t0", "tf", "duration"]
        if tiene_anomalias:
            cols += [
                "max_anom",
                "min_anom",
                "mean_anom",
                "std_anom",
                "acc_anom",
                "growth_anom",
                "decline_anom",
            ]
        logger.info(f"DataFrame creado") 
        return pd.DataFrame(columns=cols)

    df_base = pd.DataFrame(
        {"lat": lats, "lon": lons, "t0": t0s, "tf": tfs, "duration": duration}
    )

    if tiene_anomalias:
        # Desempaquetar la lista de diccionarios estadísticos en columnas del DataFrame
        df_anom = pd.DataFrame(anom_stats)
        logger.info(f"DataFrame creado") 
        return pd.concat([df_base, df_anom], axis=1)
    
    logger.info(f"DataFrame creado") 
    return df_base
    


def _procesar_pixel(serie_da: np.ndarray,
                       min_size_mhw: int = 5, 
                       max_gap: int = 2) -> np.ndarray:
    """Funcion interna detecta Olas de Calor Marinas (MHW) aplicando 
    el criterio de Hobday et al. (2016).

    Parámetros
    ----------
    serie_da: xr.DataArray
        Array booleano (o numérico) 1D (time). True/1 si supera el p90.
    min_size_mhw : int, opcional
        Duración mínima del evento en días. Por defecto 5.
    max_gap : int, opcional
        Brecha máxima de días consecutivos permitida por debajo del p90
        para unir dos eventos contiguos. Por defecto 2.

    Devuelve
    --------
    xr.DataArray
        Matriz 1D con la misma estructura original:
        - 1.0 = Presencia de MHW
        - 0.0 = Océano sin MHW
        - NaN = Tierra firme
    """
    serie_da = serie_da.copy()

    # si arranca una cadena de 1, change == 1
    # si termina una cade de 1, change == -1
    change_da = np.diff(serie_da.astype(int), prepend=0, append=0)
    starts = np.where(change_da == 1)[0]   
    ends = np.where(change_da == -1)[0]

    if len(starts) > 0:
        size_posible_waves = starts.shape[0]
        i = 1
        new_starts = []
        new_ends = []

        # calculo la duracion del primer evento
        len_chain = ends[0] - starts[0]
        while i < size_posible_waves - 1:
            # calculo duracion del evento que sigue
            len_next_chain = ends[i] - starts[i]
            # si el evento actual es una mhw
            if len_chain >= min_size_mhw:
                len_gap = starts[i] - ends[i - 1]
                # si el gap entre eventos es menor a lo esperado
                # y el siguiente evento evento es una mhw
                if len_gap <=max_gap and len_next_chain >= min_size_mhw:
                    # transformo mi gap en mhw
                    new_starts.append(starts[i - 1])
                    new_ends.append(starts[i])
                else:
                    # mi mhw no cambia
                    new_starts.append(starts[i - 1])
                    new_ends.append(ends[i - 1])
            # tomo el siguiente evento como punto de partida
            len_chain = len_next_chain
            # me paro donde inicia ese vento
            i += 1

        # me fijo si el ultimo evento es o no mhw
        if len_chain >= min_size_mhw:
            new_starts.append(starts[i - 1])
            new_ends.append(ends[i - 1])

        # Crear máscara final
        serie_final = np.zeros_like(serie_da , dtype=bool)
        for start, end in zip(new_starts, new_ends):
            serie_final[start:end] = True
    else:
        serie_final = np.zeros_like(serie_da , dtype=bool)
    
    return serie_final


def get_marine_heat_waves(
    posible_mhw_mask: xr.DataArray, min_size_mhw: int = 5, max_gap: int = 2
) -> xr.DataArray:
    """Detecta Olas de Calor Marinas (MHW) aplicando el criterio de Hobday et al. (2016).

    Parámetros
    ----------
    posible_mhw_mask: : xr.DataArray
        Array booleano (o numérico) 3D (time, lat, lon). True/1 si supera un dado percentil.
    min_size_mhw : int, opcional
        Duración mínima del evento en días. Por defecto 5.
    max_gap : int, opcional
        Brecha máxima de días consecutivos permitida por debajo del p90
        para unir dos eventos contiguos. Por defecto 2.

    Devuelve
    --------
    xr.DataArray
        Matriz 3D con la misma estructura original:
        - 1.0 = Presencia de MHW
        - 0.0 = Océano sin MHW
        - NaN = Tierra firme
    """
    logger.info("Iniciando proceso")
    # Tierra (donde todo el registro temporal sea NaN)
    es_tierra = posible_mhw_mask.isnull().all(dim="time")

    # Limpiar NaNs temporales convirtiéndolos a 0 para el análisis
    data_bool = np.nan_to_num(posible_mhw_mask.values, nan=0).astype(bool)
    logger.info("Procesando pixel por pixel (puede demorar algunos minutos)")
    # Aplicar la función _procesar_pixel a lo largo del eje del tiempo (axis=0) de manera eficiente
    waves_np = np.apply_along_axis(_procesar_pixel, 
                                   axis=0, 
                                   arr=data_bool, 
                                   min_size_mhw=min_size_mhw,
                                   max_gap=max_gap)

    # Reconstruir el DataArray e inyectar correctamente la máscara de tierra original
    logger.info("Agregando dimensiones y coordenadas al DataArray")
    waves_da = xr.DataArray(
        waves_np, 
        dims=posible_mhw_mask.dims, 
        coords=posible_mhw_mask.coords,
        name="waves"
    )
    logger.info("DataArray tipo mascara de olas de calor creado")
    return waves_da.where(~es_tierra, np.nan)

