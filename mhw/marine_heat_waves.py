import numpy as np
import pandas as pd

from scripts_santi.mhw_tools import days_from_numpytimeserie, create_index_time, redim_matrix

def get_climatology_and_anomalies(sst, time, clim_year_init=1991, clim_year_end=2020, window_time=5):
    # constantes 
    window_time = np.timedelta64(window_time, "D")
    n_mes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    dias_mes = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    
    # dimensiones
    days_of_year = sum(dias_mes)
    size_lat = sst.shape[1]
    size_lon = sst.shape[2]
    size_time = sst.shape[0]
    
    # variables
    climatology = np.zeros((days_of_year, size_lat, size_lon)) * np.nan
    percentil90 = np.zeros((days_of_year, size_lat, size_lon)) * np.nan
    anomalias = np.zeros((size_time, size_lat, size_lon)) * np.nan
    is_sup_p90 = np.zeros((size_time, size_lat, size_lon)) * np.nan
    
    day_of_year = 0
    for month in n_mes:
        for day in range(1, dias_mes[month - 1] + 1):
            """
            index_time_window = indice que toma de TODOS los anios el mes/dia +-5 dia
            index_days = indece que toma el mes-dia de TODOS los anios
            """
            if (n_mes==2) and (day==29):
                climatology[day_of_year, :, :] = climatology[day_of_year - 1, :, :].copy()
                anomalias[index_days, :, :] =  anomalias[index_days - 1, :, :].copy()
                is_sup_p90[index_days, :, :] = is_sup_p90[index_days - 1, :, :].copy()
            else:
                index_time_window, index_days = create_index_time(time.copy(), clim_year_init, clim_year_end, month, day, window_time,)
                # print(month, day, np.sum(index_time_window), np.sum(index_days))
                _sst = sst[index_time_window, : , :].copy()
                climatology[day_of_year, :, :] = np.nanmean(_sst, axis=0)
                percentil90[day_of_year, :, :] = np.nanpercentile(_sst, 90, axis=0)
        
                redim_climatology = redim_matrix(climatology[day_of_year, :, :], index_days[index_days==True].shape[0])
                redim_percentil90 = redim_matrix(percentil90[day_of_year, :, :], index_days[index_days==True].shape[0])
                anomalias[index_days, :, :] =  sst[index_days, :, :] - redim_climatology
                is_sup_p90[index_days, :, :] = sst[index_days, :, :] > redim_percentil90
            day_of_year += 1
    return climatology, percentil90, anomalias, is_sup_p90



def posible_mhw(M, t, size_time, is_wave=False):
    """
    """
    t0 = t
    while (M[t] and (t < size_time - 1)):
        t += 1 
    return t0, t


def join_waves(waves):
    """
    """
    i = 0
    nr = 1 
    while (i < waves.shape[0] - 1) and (i + nr < waves.shape[0] - 1):
        d = waves.t0.loc[i + nr] - waves.tf.loc[i]
        if d <= 2: 
            waves.tf.loc[i] = waves.tf.loc[i + nr]
            waves.drop(i + nr, inplace=True)
            nr += 1
        else:
            i += nr
            nr = 1
         
    waves.reset_index(inplace=True, drop=True)
    return waves



def get_marine_heat_waves(is_sup_p90):
    """
    """
    _lats = []
    _lons = []
    _t0s = []
    _tfs = []

    waves = np.zeros(is_sup_p90.shape) * np.nan
    size_time = is_sup_p90.shape[0]
    size_lat = is_sup_p90.shape[1]
    size_lon = is_sup_p90.shape[2]
    
    for _i in range(size_lat):
        for _j in range(size_lon):
            point_is_sup_p90 = is_sup_p90[:, _i, _j].copy()
            _t = 0
            while (_t < size_time - 4):
                _s = _t
                t0, tf = posible_mhw(point_is_sup_p90, _s, size_time)
                if (tf - t0) >= 5:
                    _lats.append(_i)
                    _lons.append(_j)
                    _t0s.append(t0)
                    _tfs.append(tf)
                waves[t0:tf + 1, _i, _j] = (tf - t0) >= 5
                _t += 1 + (tf - t0) 
    
    df_waves = pd.DataFrame({"index_lat": _lats, "index_lon": _lons, "t0":_t0s, "tf":_tfs})

    dfs = []
    coords = set(df_waves[["index_lat", "index_lon"]].itertuples(index=False, name=None))
    for _lat, _lon in coords:
        waves = df_waves[(df_waves.index_lon==_lon) & (df_waves.index_lat==_lat)].copy()
        waves.reset_index(inplace=True)
        aux = join_waves(waves.copy())
        dfs.append(aux)
    
    df_waves = pd.concat(dfs)
    
    return waves, df_waves