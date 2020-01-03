library(sf)
library(sp)
library(raster)
library(tidyverse)
library(tmap)
library(RColorBrewer)

#carga de capas y transformación
incendios <- st_read(dsn = "data/Incendios/FiredataBuffer.shp")
mun <- st_read(dsn = 'data/DivisionRD/divisionRD.gpkg', layer = 'MUNCENSO2010')
mun4326 <- st_transform(mun, crs = 4326)
prov <- st_read(dsn = 'data/DivisionRD/divisionRD.gpkg', layer = 'PROVCENSO2010')
prov4326 <- st_transform(prov, crs = 4326)
usoSuelo <- raster('data/UsoSuelo/GLOBCOVER_RD.color.tif')

#extracción de datos raster
incUsoSuelo <- raster::extract(usoSuelo, incendios,sp=TRUE)
incForestales <- subset(incUsoSuelo,incUsoSuelo$GLOBCOVER_RD.color>10 & incUsoSuelo$GLOBCOVER_RD.color<130)
incForestales$GLOBCOVER_RD.color
incForestales.df <- data.frame(incForestales)
colnames(incForestales.df)
summary(incForestales.df[[18]])

#
incForestales.sf <- st_as_sf(incForestales)
incendiosForestales <- st_intersection(incForestales.sf,mun4326)
plot(incendiosForestales['GLOBCOVER_RD.color'])
table(incendiosForestales$ENLACE)

#Conversion a Poligono
munInc <-arrange(mun4326, ENLACE)
munInc$NumIncendios <- c(table(incendiosForestales$ENLACE))
plot(munInc['NumIncendios'])

#Vecindad
munInc.sp <- as_Spatial(munInc)
colnames(munInc.sp@data)
row.names(munInc.sp) <- as.character(munInc.sp$TOPONIMIA)
munInc.nb <- poly2nb(munInc.sp, queen = TRUE)
summary(munInc.nb)
plot(munInc.sp, border="grey", lwd=0.5)
plot(munInc.nb, coordinates(munInc.sp), add=T)

#Num Vecinos
coords <- coordinates(munInc.sp)
ident <- row.names(munInc.sp)
munInc.nb.k1 <- knn2nb(knearneigh(coords, k = 1), row.names = ident)
summary(munInc.nb.k1)

card(munInc.nb.k1)
plot(munInc.sp, border="grey", lwd=0.5)
plot(munInc.nb.k1, coordinates(munInc.sp), add=T)
is.symmetric.nb(munInc.nb.k1)
dist <- unlist(nbdists(munInc.nb.k1, coords))
summary(dist)
hist(dist)
boxplot(dist)
(distmin <- min(dist)) 
(distmax <- max(dist))
indicemin <- which(distmin==dist)
ident[indicemin]
indicemax <- which(distmax==dist)
ident[indicemax]
ident[order(dist)]
 
#Ponderadoes espaciales
munInc.w.W <- nb2listw(munInc.nb)
munInc.w.W

munInc.w.B <- nb2listw(munInc.nb, style = 'B')
munInc.w.B


#Correlacion Incendios
sum(munInc$NumIncendios)
munIncPerc <- munInc %>% mutate('IncPercentage' = munInc$NumIncendios/sum(munInc$NumIncendios)*100,
         'IncPercentage_log' = log1p(munInc$NumIncendios/sum(munInc$NumIncendios)*100))

#Mapa Porcentajes
p1 <- tm_shape(munIncPerc) +
  tm_fill(col = "IncPercentage", style = 'jenks',
          palette = brewer.pal(9, name = 'Reds'), title = 'Porcentaje Incendios') +
  tm_borders(lwd = 0.5)

p2 <- tm_shape(munIncPerc) +
  tm_fill(col = "IncPercentage_log", style = 'jenks',
          palette = brewer.pal(9, name = 'Reds'), midpoint = NA, title = 'Porcentaje Incendios') +
  tm_borders(lwd = 0.5)
tmap_arrange(p1)
plot(munIncPerc['IncPercentage'])


