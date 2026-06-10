import numpy as np
import pandas as pd
import xarray as xr
from typing import Tuple


def get_climatology_and_anomalies(
    sst_da: xr.DataArray, 
    clim_year_init: int = 1991, 
    clim_year_end: int = 2020, 
    window_time: int = 5,
    percentil: float = .9
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
        Percentil a calcularse. Por defecto es 0.9

    Returns
    -------
    xr.Dataset
        Dataset que contiene 'climatology', 'percentil90', 'anomalias', e 'is_sup_p90'.
        Mantiene las coordenadas y dimensiones originales.
    """
    # 1. Construir la ventana móvil SOBRE TODA LA SERIE (para no perder el efecto borde)
    window_size = window_time * 2 + 1
    # Esto asegura que el 1 de enero de 1991 pueda "ver" a diciembre de 1990
    sst_rolled = sst_da.rolling(time=window_size, center=True).construct("window_dim")
    
    # 2. Recortar el período base climatológico YA con las ventanas construidas
    base_rolled = sst_rolled.sel(time=slice(str(clim_year_init), str(clim_year_end)))
    
    # 3. Agrupar estrcitamente por mes y día (vital para alinear los bisiestos)
    month_day = base_rolled.time.dt.strftime('%m-%d')
    climatology = base_rolled.groupby(month_day).mean(dim=["time", "window_dim"], skipna=True)
    percentil90 = base_rolled.groupby(month_day).quantile(percentil, dim=["time", "window_dim"], skipna=True)
    
    # 4. Calcular anomalías para TODA la serie original (Xarray alinea el dayofyear automáticamente)
    sst_grouped = sst_da.groupby(sst_da.time.dt.strftime('%m-%d'))
    anomalias = sst_grouped - climatology
    is_sup_p90 = sst_grouped > percentil90
    
    # 5. Limpieza de variables auxiliares que genera groupby
    anomalias = anomalias.drop_vars("dayofyear", errors="ignore")
    is_sup_p90 = is_sup_p90.drop_vars("dayofyear", errors="ignore")
    
    # Retornamos todo empaquetado prolijamente en un Dataset
    return xr.Dataset({
        "climatology": climatology,
        "percentil90": percentil90,
        "anomalias": anomalias,
        "is_sup_p90": is_sup_p90
    })


def _anomalies_brief(wave_anom: xr.DataArray) -> dict:
    """Calcula las métricas de intensidad y tasas dinámicas (crecimiento y

    decaimiento) para el ciclo de vida de una Ola de Calor Marina (MHW).

    Parámetros
    ----------
    wave_anom : np.ndarray o pd.Series
        Vector con las anomalías de temperatura de la superficie del mar (SSTa)
        durante la duración exacta del evento.

    Devuelve
    --------
    dict
        Métricas clave de la MHW (máxima, media, acumulada, tasas, etc.).
    """
    # 1. Métricas de Intensidad Clásicas
    max_anom = max(wave_anom)  # Intensidad máxima (pico del evento)
    min_anom = min(wave_anom)  # Anomalía mínima registrada en el período
    mean_anom = wave_anom.mean()  # Intensidad media de la ola
    std_anom = wave_anom.std()  # Variabilidad interna de la anomalía

    # Intensidad acumulada (°C-días): Medida directa del estrés térmico total
    acc_anom = wave_anom.sum()

    # 2. Análisis Dinámico (Evolución temporal)
    idx_max = np.argmax(wave_anom)  # Día en el que se alcanzó el pico
    n_days = wave_anom.shape[0]  # Duración total de la ola en días

    # Tasa de crecimiento: velocidad de calentamiento inicial hasta el pico (°C/día)
    # Se agrega salvaguarda por si el pico ocurre el día de inicio (evita división por 0)
    growth_anom = (
        (max_anom - wave_anom[0]) / idx_max if idx_max > 0 else np.nan
    )

    # Tasa de decaimiento: velocidad de enfriamiento desde el pico hasta el final (°C/día)
    # Se agrega salvaguarda por si el pico ocurre el último día (evita división por 0)
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

    Devuelve
    --------
    pd.DataFrame
        Catálogo/Inventario indexado donde cada fila representa un evento único.
    """
    # 1. Preparación de datos base (conversión a matrices NumPy limpias)
    data = waves_da.fillna(0).values.astype(np.int8)
    data = np.squeeze(data)

    # Controlar la existencia de anomalías de forma segura para Xarray
    tiene_anomalias = anomalies is not None
    if tiene_anomalias:
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

    # 2. Bucle espacial anidado (Píxel por píxel)
    for ilat, lat in enumerate(original_lats):
        for ilon, lon in enumerate(original_lons):
            serie = data[:, ilat, ilon]
            t = 0

            # 3. Escaneo temporal lineal dentro del píxel activo
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
                        segmento_anom = anom[t0 : tf + 1, ilat, ilon]
                        anom_stats.append(_anomalies_brief(segmento_anom))
                else:
                    t += 1

    # 4. Construcción y retorno del DataFrame de salida
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
        return pd.DataFrame(columns=cols)

    df_base = pd.DataFrame(
        {"lat": lats, "lon": lons, "t0": t0s, "tf": tfs, "duration": duration}
    )

    if tiene_anomalias:
        # Desempaquetar la lista de diccionarios estadísticos en columnas del DataFrame
        df_anom = pd.DataFrame(anom_stats)
        return pd.concat([df_base, df_anom], axis=1)

    return df_base
    


def _procesar_pixel(serie: np.ndarray,
                       min_dur: int = 5, 
                       max_gap: int = 2) -> np.ndarray:
    """Funcion interna detecta Olas de Calor Marinas (MHW) aplicando 
    el criterio de Hobday et al. (2016).

    Parámetros
    ----------
    serie : xr.DataArray
        Array booleano (o numérico) 1D (time). True/1 si supera el p90.
    min_dur : int, opcional
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
    size_t = serie.shape[0]
    output = np.zeros(size_t, dtype=np.float32)
    t = 0

    while t <= size_t - min_dur:
        # Condición de inicio: requiere 'min_dur' días consecutivos sobre el percentil
        if serie[t] and np.all(serie[t : t + min_dur]):
            t0 = t
            t += min_dur

            # Extender el evento tolerando caídas breves (gaps cortos)
            while t < size_t:
                if not serie[t]:
                    # Si los próximos días sin percentil superan el max_gap, se corta el evento
                    if t + max_gap < size_t and not np.any(
                        serie[t : t + max_gap + 1]
                    ):
                        break
                t += 1

            output[t0:t] = 1.0
        else:
            t += 1

    return output


def get_marine_heat_waves(
    is_sup_p90: xr.DataArray, min_dur: int = 5, max_gap: int = 2
) -> xr.DataArray:
    """Detecta Olas de Calor Marinas (MHW) aplicando el criterio de Hobday et al. (2016).

    Parámetros
    ----------
    is_sup_p90 : xr.DataArray
        Array booleano (o numérico) 3D (time, lat, lon). True/1 si supera el p90.
    min_dur : int, opcional
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
    # 1. Eliminar dimensiones pequenas (zlev, depth)
    is_sup_p90 = np.squeeze(is_sup_p90)
    
    # 2. Identificar la máscara estática de tierra (donde todo el registro temporal sea NaN)
    es_tierra = is_sup_p90.isnull().all(dim="time")

    # 3. Limpiar NaNs temporales convirtiéndolos a 0 para el análisis numérico de la serie
    data_bool = np.nan_to_num(is_sup_p90.values, nan=0).astype(bool)
    size_t = data_bool.shape[0]

    # 4. Aplicar la función _procesar_pixel a lo largo del eje del tiempo (axis=0) de manera eficiente
    waves_np = np.apply_along_axis(_procesar_pixel, 
                                   axis=0, 
                                   arr=data_bool, 
                                   min_dur=min_dur,
                                   max_gap=max_gap)

    # 5. Reconstruir el DataArray e inyectar correctamente la máscara de tierra original
    waves_da = xr.DataArray(
        waves_np, 
        dims=is_sup_p90.dims, 
        coords=is_sup_p90.coords,
        name="waves"
    )

    return waves_da.where(~es_tierra, np.nan)

