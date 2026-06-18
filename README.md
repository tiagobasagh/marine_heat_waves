# Procesamiento y Análisis de Olas de Calor Marinas (MHW)

Este repositorio contiene un conjunto de herramientas en Python diseñado para la detección, procesamiento estadístico y análisis de tendencias espaciotemporales de Olas de Calor Marinas (MHW, por sus siglas en inglés) en datos de temperatura superficial del mar (SST). El algoritmo implementado sigue Hobday et al.(2016) y fue utilizado para la publicacion M.L.Clara et Al (2026).

---

## Estructura del Repositorio

El proyecto se organiza en dos carpetas principales que separan la lógica de cálculo numérico de los entornos de visualización y aprendizaje:

### 1. Módulos de Cálculo (`./mhw`)
Contiene los scripts de Python con las funciones científicas y matemáticas fundamentales. Están escritos de forma modular para que puedan ser importados desde cualquier script o notebook:

* **`mhw_core.py`**: Es el núcleo del repositorio. Contiene las funciones para calcular la climatología media diaria y el percentil 90 móvil dentro de una ventana temporal (por defecto de 11 días), identificar anomalías de temperatura, y ejecutar el algoritmo de detección de eventos que filtra y une las olas de calor/frio según su duración mínima y los días intermedios de tregua.
* **`stats_and_trends.py`**: Contiene las herramientas para el análisis posterior a la detección. Incluye funciones vectorizadas con Xarray para calcular tendencias lineales de largo plazo (pendientes, interceptos y significancia estadística a través de p-valores) píxel por píxel sobre la grilla marina, además de extraer estadísticos clave (picos máximos, intensidades medias e integradas) para cada evento detectado.
* **`datasets_utils.py`**: Gestiona la persistencia de la información, centralizando las rutas de almacenamiento y proveyendo funciones para guardar de forma estandarizada los resultados en formatos científicos de tipo NetCDF (`.nc`) o tablas estructuradas en CSV (`.csv`).
* **`time_utils.py`**: Contiene funciones auxiliares para el procesamiento cronológico, permitiendo la extracción de componentes anuales y la clasificación estacional.

### 2. Guías de Flujo de Trabajo (`./notebooks`)
Incluye cuadernos interactivos de Jupyter (Jupyter Notebooks) diseñados a modo de tutorial secuencial. Están pensados para guiar al investigador paso a paso, facilitando la comprensión del flujo de trabajo científico sin requerir un conocimiento avanzado de programación:

* **`Tutorial I - Liberia MHW M.L Clara et Al.ipynb`**: Primer paso del flujo. Ilustra cómo cargar los campos de temperatura, calcular la climatología base, aplicar el percentil umbral, identificar las ventanas temporales que constituyen una ola de calor y exportar tanto la matriz tridimensional resultante como un catálogo inicial de eventos en formato de tabla. Permite una segunda etapa de ejemplo donde agrupa los eventos detectados por año calendario/estacion/mes para construir indicadores climáticos (del periodo) por cada píxel del mapa (frecuencia de eventos, duración media anual, anomalía máxima absoluta y anomalia maxima media).
* **`Tutorial II - Ejemplos Plots M.L Clara et Al.ipynb`**: Ejemplo de visualizaciones y analisis posibles a partir de los dataset generados.

---

## Requisitos del Sistema y Dependencias

El repositorio aprovecha el ecosistema científico estándar de Python. 

* **Xarray**: Utilizada para la gestión de variables georreferenciadas con coordenadas de latitud, longitud y tiempo. Permite realizar operaciones matemáticas en bloque sobre toda la grilla marina de forma eficiente.
* **NumPy**: Soporte subyacente para el cálculo matricial y la manipulación de arreglos numéricos optimizados.
* **Pandas**: Manejo y filtrado de los catálogos de eventos en formato tabular (similar a planillas de cálculo de Excel).
* **SciPy**: Módulo estadístico empleado específicamente para resolver las regresiones lineales y determinar la significancia física de las tendencias analizadas.
* **Matplotlib**:
* **Cartopy**:
* **CMocean**:
---

## Instrucciones de Instalación

???

## Como citar

???

## Rerefencias

???