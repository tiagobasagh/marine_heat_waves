import pandas as pd 
import numpy as np

from typing import Tuple

def extract_year(date)->int:
    """
    
    """
    if isinstance(date, pd.Timestamp):
        return date.year
    elif isinstance(date, str):
        return int(date.split("-")[0])
    else:
        raise TypeError(f"date es {type(date)} se espera pd.Timestamp/str")


def extract_month(date)->int:
    """
    
    """
    if isinstance(date, pd.Timestamp):
        return date.month
    elif isinstance(date, str):
        return int(date.split("-")[1])
    else:
        raise TypeError(f"date es {type(date)} se espera pd.Timestamp/str")


def extract_season(date, hemisferio="sur" ) -> str:
    """
    Devuelve la estación del año dependiendo el hemiferio dada una fecha.

    Parameters
    ----------
    date

    Returns
    -------
    str
        Nombre de la estación ('summer', 'autumn', 'winter', 'spring') o 'error'.
    """
    month = extract_month(date)
    if hemisferio.lower()=="sur" or hemisferio.lower()=="south":
        if month in (1, 2, 3):
            return "summer"
        elif month in (4, 5, 6):
            return "autumn"
        elif month in (7, 8, 9):
            return "winter"
        elif month in (10, 11, 12):
            return "spring"
        raise "Error! Mes no valido"