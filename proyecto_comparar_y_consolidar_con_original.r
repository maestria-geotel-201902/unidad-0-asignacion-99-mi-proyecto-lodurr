library(sf)
library(sp)
library(raster)
library(tidyverse)
library(tmap)
library(RColorBrewer)
library(lmtest)
library(spdep)
library(parallel)

#carga de capas y transformación
incendios <- st_read(dsn = "data/Incendios/FiredataBuffer.shp")
mun <- st_read(dsn = 'data/DivisionRD/divisionRD.gpkg', layer = 'MUNCENSO2010')
mun4326 <- st_transform(mun, crs = 4326)
usoSuelo <- raster('data/UsoSuelo/GLOBCOVER_RD.color.tif')

#extracción de datos raster
incUsoSuelo <- raster::extract(usoSuelo, incendios,sp=TRUE)
incForestales <- subset(incUsoSuelo,incUsoSuelo$GLOBCOVER_RD.color>20 & incUsoSuelo$GLOBCOVER_RD.color<130)
incForestales.df <- data.frame(incForestales)
summary(incForestales.df[[18]])
colnames(incForestales.df)

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
munIncPerc <- munInc %>%st_centroid() %>%  mutate( 
  'IncPercentage' = munInc$NumIncendios/sum(munInc$NumIncendios)*100,
  'IncPercentage_log' = log1p(munInc$NumIncendios/sum(munInc$NumIncendios)*100),
  'IncPercentage_tukey' = rcompanion::transformTukey(NumIncendios/sum(munInc$NumIncendios)*100, plotit = F),
  'IncPercentage_tukey_lambda' = rcompanion::transformTukey(NumIncendios/sum(munInc$NumIncendios)*100, returnLambda = T),
  'AreaKm2' = as.numeric((st_area(munInc)/1000000)),
  'IncXArea' = (munInc$NumIncendios/AreaKm2),
  'IncXArea_log' = log1p(munInc$NumIncendios/AreaKm2),
  'IncXArea_tukey' = rcompanion::transformTukey(NumIncendios/AreaKm2, plotit = F),
  'IncXArea_tukey_lambda' = rcompanion::transformTukey(NumIncendios/AreaKm2, returnLambda = T),
  
  x=unlist(map(geom,1)), y=unlist(map(geom,2)))  %>% 
  st_drop_geometry() 
#     lambda      W Shapiro.p.value
# 409    0.2 0.9944          0.8194
# 
# if (lambda >  0){TRANS = x ^ lambda} 
# if (lambda == 0){TRANS = log(x)} 
# if (lambda <  0){TRANS = -1 * x ^ lambda} 

#Join 
munIncPercPol <- munInc %>%
  merge(munIncPerc, all.y=TRUE)
plot(munIncPercPol['ENLACE'])

#Mapa Porcentajes
p1 <- tm_shape(munIncPercPol) +
  tm_fill(col = "IncPercentage", style = 'jenks',
          palette = brewer.pal(9, name = 'Reds'), title = 'Porcentaje Incendios') +
  tm_borders(lwd = 0.5)

p2 <- tm_shape(munIncPercPol) +
  tm_fill(col = "IncPercentage_log", style = 'jenks',
          palette = brewer.pal(9, name = 'Reds'), midpoint = NA, title = 'Porcentaje Incendios') +
  tm_borders(lwd = 0.5)

tmap_arrange(p1,p2)

#qq 
qqnorm(munIncPerc$IncPercentage)
shapiro.test(munIncPerc$IncPercentage)

qqnorm(munIncPerc$IncPercentage_log)
shapiro.test(munIncPerc$IncPercentage_log)

qqnorm(munIncPerc$IncPercentage_tukey)
shapiro.test(munIncPerc$IncPercentage_tukey)

munIncPerc %>% lm(IncPercentage ~ x, .) %>% bptest()
munIncPerc %>% lm(IncPercentage ~ y, .) %>% bptest()
munIncPerc %>% lm(IncPercentage_log ~ x, .) %>% bptest()
munIncPerc %>% lm(IncPercentage_log ~ y, .) %>% bptest()
munIncPerc %>% lm(IncPercentage_tukey ~ x, .) %>% bptest()
munIncPerc %>% lm(IncPercentage_tukey ~ y, .) %>% bptest()


match(attr(munInc.w.W$neighbours, "region.id"), munIncPerc$TOPONIMIA)==1:155

(gmoranw <- moran.test(x = munIncPerc$'IncPercentage_log', listw = munInc.w.W ))
(gmoranb <- moran.test(x = munIncPerc$'IncPercentage_log', listw = munInc.w.B))

(gmoranw <- moran.test(x = munIncPerc$IncXArea_tukey, listw = munInc.w.W ))
(gmoranb <- moran.test(x = munIncPerc$IncXArea_tukey, listw = munInc.w.B))


moran.plot(x = munIncPerc$IncPercentage_log, listw = munInc.w.W)
moran.plot(x = munIncPerc$IncPercentage_tukey, listw = munInc.w.W)

source('lisaclusters.R')
lisamap(objesp = munIncPercPol,
        var = 'IncPercentage_log',
        pesos = munInc.w.W,
        tituloleyenda = 'Significancia\n("x-y", léase\ncomo "x"\nrodeado de "y"',
        leyenda = T,
        anchuratitulo = 1000,
        tamanotitulo = 16,
        fuentedatos = '',
        titulomapa = paste0('Clusters LISA de Porcentaje de Incendios Forestales'))

lisamap(objesp = munIncPercPol,
        var = 'IncPercentage_tukey',
        pesos = munInc.w.W,
        tituloleyenda = 'Significancia\n("x-y", léase\ncomo "x"\nrodeado de "y"',
        leyenda = T,
        anchuratitulo = 1000,
        tamanotitulo = 16,
        fuentedatos = '',
        titulomapa = paste0('Clusters LISA de Porcentaje de Incendios Forestales (Tukey Trans.)'))


# Mapa Porcentaje por Km2
p3 <- tm_shape(munIncPercPol) +
  tm_fill(col = "IncXArea", style = 'jenks',
          palette = brewer.pal(9, name = 'Reds'), title = 'Porcentaje Incendios por Km2') +
  tm_borders(lwd = 0.5)
p4 <- tm_shape(munIncPercPol) +
  tm_fill(col = "IncXArea_log", style = 'jenks',
          palette = brewer.pal(9, name = 'Reds'), title = 'Porcentaje Incendios por Km2') +
  tm_borders(lwd = 0.5)

tmap_arrange(p3,p4)

#qq 
qqnorm(munIncPerc$IncXArea)
shapiro.test(munIncPerc$IncXArea)

qqnorm(munIncPerc$IncXArea_log)
shapiro.test(munIncPerc$IncXArea_log)

munIncPerc %>% lm(IncXArea ~ x, .) %>% bptest()
munIncPerc %>% lm(IncXArea ~ y, .) %>% bptest()
munIncPerc %>% lm(IncXArea_log ~ x, .) %>% bptest()
munIncPerc %>% lm(IncXArea_log ~ y, .) %>% bptest()

match(attr(munInc.w.W$neighbours, "region.id"), munIncPerc$TOPONIMIA)==1:155

(gmoranw <- moran.test(x = munIncPerc$'IncXArea_log', listw = munInc.w.W ))
(gmoranb <- moran.test(x = munIncPerc$'IncXArea_log', listw = munInc.w.B))

(gmoranw <- moran.test(x = munIncPerc$IncXArea_tukey, listw = munInc.w.W ))
(gmoranb <- moran.test(x = munIncPerc$IncXArea_tukey, listw = munInc.w.B))

moran.plot(x = munIncPerc$IncXArea_log, listw = munInc.w.W)
moran.plot(x = munIncPerc$IncPercentage_tukey, listw = munInc.w.W)

source('lisaclusters.R')
lisamap(objesp = munIncPercPol,
        var = 'IncXArea_log',
        pesos = munInc.w.W,
        tituloleyenda = 'Significancia\n("x-y", léase\ncomo "x"\nrodeado de "y"',
        leyenda = T,
        anchuratitulo = 1000,
        tamanotitulo = 16,
        fuentedatos = '',
        titulomapa = paste0('Clusters LISA de Porcentaje de Incendios Forestales'))

lisamap(objesp = munIncPercPol,
        var = 'IncXArea_tukey',
        pesos = munInc.w.W,
        tituloleyenda = 'Significancia\n("x-y", léase\ncomo "x"\nrodeado de "y"',
        leyenda = T,
        anchuratitulo = 1000,
        tamanotitulo = 16,
        fuentedatos = '',
        titulomapa = paste0('Clusters LISA de Porcentaje de Incendios Forestales (Tukey Trans.)'))

#correlación variables world clim

# * Stack
wclayerspath <- list.files(path = 'data/WorldClim/', pattern = '*.tif', recursive = T, full.names = T)
wcstack <- stack(wclayerspath)

# * Add month field
incendiosForestales$month <- format(as.Date(incendiosForestales$ACQ_DATE), "%m")

# * Add unique field
incendiosForestales$unique <- 1:nrow(incendiosForestales)

# * Extract values from WC corresponding to the month of each point
system.time(
  foo <- sapply(1:20, function(x) {
    sp <- incendiosForestales[x,]
    m <- ifelse(nchar(sp$month)==1, paste0('0', sp$month), sp$month)
    e <- raster::extract(wcstack[[grep(paste0(m,'$'), names(wcstack))]], sp, sp=T)
    d <- e@data[,c('unique', 'month', grep('RD_wc2.*', colnames(e@data), value = T))]
    return(d)
  }, simplify = F)
)
#Laptop 8 cores, 8 GB
#   user  system elapsed 
# 20.801   0.008  20.807
#Workstation, 8 cores, 64 GB
#  user  system elapsed 
# 12.036   0.012  12.047 
bar <- inner_join(incendiosForestales, bind_rows(foo), by = c('unique', 'month'))
bar

#Parallel
UseCores <- detectCores() - 1
cl <- makeCluster(UseCores, outfile=paste0('info_parallel.log'))
clusterExport(cl, list('incendiosForestales', 'wcstack'))
system.time(
  foo <- parSapply(cl, 1:20, function(x) {
    sp <- sf::as_Spatial(incendiosForestales[x,])
    m <- ifelse(nchar(sp$month)==1, paste0('0', sp$month), sp$month)
    e <- raster::extract(wcstack[[grep(paste0(m,'$'), names(wcstack))]], sp, sp=T)
    d <- e@data[,c('unique', 'month', grep('RD_wc2.*', colnames(e@data), value = T))]
    return(d)
  }, simplify = F)
)
#Laptop 8 cores, 8 GB
#  user  system elapsed 
# 0.009   0.008   2.775
#Workstation, 8 cores, 64 GB
#  user  system elapsed 
# 0.009   0.000   1.552 
#Estimated entire job:
1.552/20*nrow(incendiosForestales)/60
# [1] 60.96256
stopCluster(cl)
bar <- inner_join(incendiosForestales, bind_rows(foo), by = c('unique', 'month'))
bar
#

#Parallel
UseCores <- detectCores() - 1
cl <- makeCluster(UseCores, outfile=paste0('info_parallel.log'))
clusterExport(cl, list('incendiosForestales', 'wcstack'))

## Not run: 
system.time(
  foo <- parSapply(cl, 1:nrow(incendiosForestales), function(x) {
    sp <- sf::as_Spatial(incendiosForestales[x,])
    m <- ifelse(nchar(sp$month)==1, paste0('0', sp$month), sp$month)
    e <- raster::extract(wcstack[[grep(paste0(m,'$'), names(wcstack))]], sp, sp=T)
    d <- e@data[,c('unique', 'month', grep('RD_wc2.*', colnames(e@data), value = T))]
    return(d)
  }, simplify = F)
)
## End(Not run)
#Workstation, 8 cores, 64 GB



stopCluster(cl)
bar <- inner_join(incendiosForestales, bind_rows(foo), by = c('unique', 'month'))
bar
#