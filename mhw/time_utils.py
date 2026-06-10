import numpy as np
import datetime as dt
from typing import Tuple

def extract_year_from_str_date(date: str)->int:
    """
    
    """
    return int(date.split("-")[0])


def extract_season_from_str_date(date: str) -> str:
    """
    Devuelve la estación del año (hemisferio sur) correspondiente a un mes numérico.

    Parameters
    ----------
    n : int
        Número del mes (1 a 12).

    Returns
    -------
    str
        Nombre de la estación ('summer', 'autumn', 'winter', 'spring') o 'error'.
    """
    month = int(date.split("-")[1])
    if month in (1, 2, 3):
        return "summer"
    elif month in (4, 5, 6):
        return "autumn"
    elif month in (7, 8, 9):
        return "winter"
    elif month in (10, 11, 12):
        return "spring"
    return "Mes no valido"