# Procesamiento y Análisis de Olas de Calor Marinas (MHW)

Este repositorio contiene un conjunto de herramientas en Python diseñado para la detección, procesamiento estadístico y análisis de tendencias espaciotemporales de Olas de Calor Marinas (MHW, por sus siglas en inglés) en datos de temperatura superficial del mar (SST). El algoritmo implementado sigue los criterios metodológicos estandarizados a nivel internacional por Hobday et al. (2016).

El código está estructurado para trabajar de forma eficiente con grandes volúmenes de datos georreferenciados (como archivos NetCDF provenientes de sensores satelitales o modelos climáticos) utilizando bibliotecas de alto rendimiento para el análisis de matrices multidimensionales.

---

## Estructura del Repositorio

El proyecto se organiza en dos carpetas principales que separan la lógica de cálculo numérico de los entornos de visualización y aprendizaje:

### 1. Módulos de Cálculo (`./mhw`)
Contiene los scripts de Python con las funciones científicas y matemáticas fundamentales. Están escritos de forma modular para que puedan ser importados desde cualquier script o notebook:

* **`mhw_core.py`**: Es el núcleo del repositorio. Contiene las funciones para calcular la climatología media diaria y el percentil 90 móvil dentro de una ventana temporal (por defecto de 11 días), identificar anomalías de temperatura, y ejecutar el algoritmo de detección de eventos que filtra y une las olas de calor según su duración mínima y los días intermedios de tregua.
* **`stats_and_trends.py`**: Contiene las herramientas para el análisis posterior a la detección. Incluye funciones vectorizadas con Xarray para calcular tendencias lineales de largo plazo (pendientes, interceptos y significancia estadística a través de p-valores) píxel por píxel sobre la grilla marina, además de extraer estadísticos clave (picos máximos, intensidades medias e integradas) para cada evento detectado.
* **`array_utils.py`**: Proveedor de utilidades numéricas de bajo nivel para la manipulación segura de matrices en NumPy, tales como la conversión de máscaras lógicas a valores flotantes nulos (NaN), delimitación de rangos numéricos y reescalado dimensional de campos climatológicos espaciales al eje temporal.
* **`datasets_utils.py`**: Gestiona la persistencia de la información, centralizando las rutas de almacenamiento y proveyendo funciones para guardar de forma estandarizada los resultados en formatos científicos de tipo NetCDF (`.nc`) o tablas estructuradas en CSV (`.csv`).
* **`time_utils.py`**: Contiene funciones auxiliares para el procesamiento cronológico, permitiendo la extracción de componentes anuales y la clasificación estacional de las fechas analizadas adaptada al Hemisferio Sur.

### 2. Guías de Flujo de Trabajo (`./notebooks`)
Incluye cuadernos interactivos de Jupyter (Jupyter Notebooks) diseñados a modo de tutorial secuencial. Están pensados para guiar al investigador paso a paso, facilitando la comprensión del flujo de trabajo científico sin requerir un conocimiento avanzado de programación:

* **`marine_heat_waves.ipynb`**: Primer paso del flujo. Ilustra cómo cargar los campos de temperatura, calcular la climatología base, aplicar el percentil umbral, identificar las ventanas temporales que constituyen una ola de calor y exportar tanto la matriz tridimensional resultante como un catálogo inicial de eventos en formato de tabla.
* **`anual_statistics.ipynb`**: Segundo paso del flujo. Agrupa los eventos detectados por año calendario para construir indicadores climáticos anuales por cada píxel del mapa (frecuencia anual de eventos, duración media anual, anomalía máxima absoluta y carga térmica acumulada).
* **`anual_trend.ipynb`**: Tercer paso del flujo. Aplica análisis de regresión lineal sobre los indicadores anuales generados en el paso anterior. Muestra cómo mapear espacialmente la tasa de cambio espacial (por ejemplo, el incremento o descenso en la cantidad de eventos por década) y filtrar cartográficamente aquellas regiones con significancia estadística comprobada.

---

## Requisitos del Sistema y Dependencias

El repositorio aprovecha el ecosistema científico estándar de Python. No es necesario conocer la arquitectura interna de estas librerías para utilizar las herramientas, ya que la complejidad matemática se maneja internamente:

* **Xarray**: Utilizada para la gestión de variables georreferenciadas con coordenadas de latitud, longitud y tiempo. Permite realizar operaciones matemáticas en bloque sobre toda la grilla marina de forma eficiente.
* **NumPy**: Soporte subyacente para el cálculo matricial y la manipulación de arreglos numéricos optimizados.
* **Pandas**: Manejo y filtrado de los catálogos de eventos en formato tabular (similar a planillas de cálculo de Excel).
* **SciPy**: Módulo estadístico empleado específicamente para resolver las regresiones lineales y determinar la significancia física de las tendencias analizadas.
* **Matplotlib**: Utilizada en las notebooks para la generación de la cartografía y la representación visual de los mapas de tendencias.

---

## Instrucciones de Instalación

Para asegurar el correcto funcionamiento del repositorio y evitar conflictos con otras herramientas de software instaladas en su computadora, se recomienda realizar la instalación dentro de un entorno virtual limpio de Python.

### Paso 1: Clonar o descargar el repositorio
Descargue el código fuente en su máquina local utilizando Git o descargando el archivo comprimido desde la plataforma de alojamiento:

```bash
git clone <url-del-repositorio>
cd <nombre-de-la-carpeta-del-repositorio>