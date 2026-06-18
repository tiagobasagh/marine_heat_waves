import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats

from typing import Tuple

def trend_matrix(
    da: xr.DataArray, 
    time_dim: str = "time"
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    """Calcula la tendencia lineal (slope, p-value, intercept) 
   """

    # 1. Definimos la función interna que operará sobre el vector temporal de cada píxel
    def _pixel_trend(y, x):
        y_nans = np.isnan(y)

        # si es tierra o si hay menos de 3 años con datos, devolvemos NaN
        if y_nans.all() or np.sum(~y_nans) < 3:
            return np.nan, np.nan, np.nan

        # Nos quedamos solo con los años válidos (removemos nans=True)
        y_valid = y[~y_nans]
        x_valid = x[~y_nans]

        # centramos 
        x_centered = x_valid - np.mean(x_valid)

        slope, intercept, r_value, p_value, std_err = stats.linregress(
            x_centered, y_valid
        )

        return slope, p_value, intercept

    # apply_ufunc se encarga de aplicar la función "pixel por pixel"
    trend, trend_pvalue, trend_intercept = xr.apply_ufunc(
        _pixel_trend,
        da,
        da[time_dim].values,
        input_core_dims=[[time_dim], [time_dim]],  # Sincroniza las dimensiones
        output_core_dims=[[], [], []],  # Tres mapas de salida escalares
        vectorize=True,  # Bucle interno veloz en NumPy
    )

    return trend, trend_pvalue, trend_intercept
