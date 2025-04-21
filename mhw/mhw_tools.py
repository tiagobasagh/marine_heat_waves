import numpy as np 
import datetime as dt


def days_from_numpytimeserie(time, year_o=1970):
    """
    """
    years = time.astype('datetime64[Y]').astype(int) + year_o
    months = time.astype('datetime64[M]').astype(int) % 12 + 1
    K = 10**(-9)/60/60/24
    days = []
    doy = []
    
    for i, t in enumerate(time): 
        t0 = np.datetime64(dt.datetime(years[i], months[i], 1)) 
        ti = np.datetime64(dt.datetime(years[i], 1, 1))
        delta = (t - t0)
        _delta_doy = (t - ti)
        days.append(np.trunc(delta.astype(int) * K ) + 1)
        doy.append(np.trunc(_delta_doy.astype(int) * K ) + 1)
    
    days = np.asarray(days)
    doy = np.asarray(doy)
    return days, doy


def create_index_time(time, year_init, year_end, month, day, window_time):
    """
    Esto me permite crear un indice de dias que me interasan para crear la climatologia  
    """
    years = time.astype('datetime64[Y]').astype(int) + 1970
    months = time.astype('datetime64[M]').astype(int) % 12 + 1
    days, _ = days_from_numpytimeserie(time)
    cond_time_window = np.full(time.shape, False, dtype=bool)
    cond_years = (year_init <= years) & (year_end >= years)
    cond_month = months==month
    cond_day = days==day 
    aux_time = time[cond_month & cond_day & cond_years]
    delta = np.timedelta64(window_time, "D")
    for t in aux_time:
        t0 = t - delta
        tf = t + delta
        cond_aux_time = (t0 <= time) & (tf >= time)
        cond_time_window = cond_time_window | cond_aux_time
    return cond_time_window, cond_month & cond_day


def redim_matrix(M, size):
    return np.repeat(M[np.newaxis, :, :], size, axis=0)


def get_main_vars_matrix(df, size=(1, 1, 1), t0=1982):
    nanmatrix = np.zeros(size) * np.nan

    count_waves_spring = nanmatrix.copy()
    mean_duration_spring = nanmatrix.copy()
    max_intensity_spring = nanmatrix.copy()
    mean_intensity = nanmatrix.copy()
    for index, row in df.iterrows():
        count_waves[row.year - t0, row.index_lat, row.index_lon] = row.duration_count
        mean_duration[row.year - t0, row.index_lat, row.index_lon] = row.duration_mean
        max_intensity[row.year - t0, row.index_lat, row.index_lon] = row.anom_max_max
        mean_intensity[row.year - t0, row.index_lat, row.index_lon] = row.anom_mean_mean

    return count_waves, mean_duration, max_intensity, mean_intensity


def trend_matrix(df, varname, timename, lat_size, lon_size):
    nanmatrix = np.zeros((lat_size, lon_size)) * np.nan
    trend = nanmatrix.copy()
    trend_pvalue = nanmatrix.copy()
    trend_intercept = nanmatrix.copy()
    for _lat in range(lat_size):
        for _lon in range(lon_size):
            pixel = df[(df.index_lon==_lon) & (df.index_lat==_lat)].copy()
            if pixel.shape[0] > 2:
                X = pixel[timename].to_numpy()
                y = pixel[varname].to_numpy()
                
                X2 = [t for t in range(X.min(), X.max() + 1)]
                y2 = [ 0 if X[X==t].shape[0]==0  else y[X==t][0] for t in X2]
                xmean = np.nanmean(X2)
                X2 = X2 - xmean
                #X2 = np.array([np.ones(X.shape), X]).T
                #r = linalg.lstsq(X2, y)
                #print(np.sum(np.isnan(y)))
                slope, intercept, r_value, p_value, std_err = stats.linregress(X2, y2)
                trend[_lat, _lon] = slope
                trend_pvalue[_lat, _lon] = p_value
                trend_intercept[_lat, _lon] = intercept
                #print(r[0][1] -  slope)
    return trend, trend_pvalue, trend_intercept