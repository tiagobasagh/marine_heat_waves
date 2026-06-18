# Procesamiento y Análisis de Olas de Calor Marinas (MHW)

Este repositorio contiene un conjunto de herramientas en Python diseñado para la detección, procesamiento estadístico y análisis de tendencias espaciotemporales de Olas de Calor Marinas (MHW, por sus siglas en inglés) en datos de temperatura superficial del mar (SST). El algoritmo implementado sigue Hobday et al.(2016) y fue utilizado para la publicacion M. Luz Clara et al (2026).

---

## Estructura del Repositorio

El proyecto se organiza en dos carpetas principales que separan la lógica de cálculo numérico de los entornos de visualización y aprendizaje:

### 1. Módulos de Cálculo (`./mhwpy`)
Contiene los scripts de Python con las funciones científicas y matemáticas fundamentales. Están escritos de forma modular para que puedan ser importados desde cualquier script o notebook:

* **`mhw_core.py`**: Es el núcleo del repositorio. Contiene las funciones para calcular la climatología media diaria y el percentil 90 móvil dentro de una ventana temporal (por defecto de 11 días), identificar anomalías de temperatura, y ejecutar el algoritmo de detección de eventos que filtra y une las olas de calor/frio según su duración mínima y los días intermedios de tregua.
* **`stats_and_trends.py`**: Contiene las herramientas para el análisis posterior a la detección. Incluye funciones vectorizadas con Xarray para calcular tendencias lineales de largo plazo (pendientes, interceptos y significancia estadística a través de p-valores) píxel por píxel sobre la grilla marina, además de extraer estadísticos clave (picos máximos, intensidades medias e integradas) para cada evento detectado.
* **`datasets_utils.py`**: Centraliza las rutas de almacenamiento y provey funciones para guardar de forma estandarizada los resultados en formatos tipo NetCDF (`.nc`) o tipo CSV (`.csv`).
* **`time_utils.py`**: Contiene funciones auxiliares para el procesamiento cronológico, permitiendo la extracción de componentes anuales, mensuales y la clasificación estacional.

### 2. Guías de Flujo de Trabajo (`./notebooks`)
Incluye cuadernos interactivos de Jupyter (Jupyter Notebooks) diseñados a modo de tutorial secuencial. Están pensados para guiar al investigador paso a paso, facilitando la comprensión del flujo de trabajo científico sin requerir un conocimiento avanzado de programación:

* **`Tutorial I - Liberia MHW.ipynb`**: Ilustra cómo cargar los datos de temperatura superficial. Calcular la climatología base, aplicar el percentil umbral, calculo de anomalias. Muestra como identificar una ola de calor y exportar tanto la matriz tridimensional (time, lat, lon) resultante como un catálogo inicial de eventos en formato de tabla. Muestra ejemplo de calculo de estadisticos, agrupando los eventos detectados por año calendario/estacion/mes para construir indicadores climáticos (del periodo) por cada píxel del mapa.
* **`Tutorial II - Ejemplos Plots M. Luz Clara et al.ipynb`**: Ejemplo de visualizaciones y analisis posibles a partir de los dataset generados. Graficos presentes en la publicacion M. Luz Clara et al [2026].
---

## Requisitos del Sistema y Dependencias

El repositorio aprovecha el ecosistema científico estándar de Python. 

* **Xarray**: 
* **NumPy**: 
* **Pandas**: 
* **SciPy**: 
* **Matplotlib**:
* **Cartopy**:
* **CMocean**:
* **Seaborn**
---

## Instrucciones de Instalación

* **clonar el repositorio**
git clone git@github.com:tiagobasagh/mhwpy.git
* **Ingresar a la carpeta**
cd mhwpy
* **installar usnado** 
pip install -e .

## Como citar


## Rerefencias

