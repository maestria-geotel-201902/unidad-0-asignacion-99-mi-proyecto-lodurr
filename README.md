# Desarrollo de Proyecto Final

Basado en lo visto durante la asignatura Análisis Espacial se pretende responder a las siguientes preguntas 

* Cuál es la distribución espacial de los puntos de calor - Hot Spots para los incendios forestales en la República Dominicana basados en la correlación entre las variables climáticas y la ubicación espacial de estos?

Para él o las áreas con mayor índice de influencia.
* Qué relación existe entre este resultado y las elevaciones del área?
* Qué relación existe entre este resultado y las pendientes del área?
* Qué relación existe entre este resultado y la cobertura arbórea del área?
* Qué relación existe entre este resultado y el uso de suelo del área?

## Data a usar
 Para el desarrollo de este trabajo se estará usando la siguiente información
 
 * Incendios forestales *"https://firms.modaps.eosdis.nasa.gov/download/create.php"* (Información de los ultimos 10 años para el pais)
 * Variables climaticas *"http://www.worldclim.org/"*  (10m)
   * Temperatura Minima       *(Cº)*
   * Temperatura Maxima       *(Cº)*
   * Temperatura Promedio     *(Cº)*
   * Precipitacion            *(mm)*
   * Radiación Solar          *(Kj m-2 day-1)*
   * Velocidad del viento     *(m s-1)*
   * Presion de Vapor de agua *(KPa)*
 * Data de OpenStreetMaps *https://download.geofabrik.de/*
   * LandUse (Industrial, Quarry
   * Buildings

Es de resaltar que para las siguiente data se realizará una selección basada en los resultados del análisis de hot-spots realizado.
 
 * Elevaciones *"https://www.eorc.jaxa.jp/ALOS/en/aw3d30/index.htm"*
 * Pendientes -> Se realizará un análisis de pendientes basado el modelo DEM (Digital Elevation Model) descargado de la JAXA
 * Cobertura arbórea *http://earthenginepartners.appspot.com/science-2013-global-forest*
 * Uso de suelo *http://due.esrin.esa.int/page_globcover.php*
 
## Desarrollo
 * El archivo [proyecto.Rmd](proyecto.Rmd) contendrá la información referente al desarrollo del mismo.
