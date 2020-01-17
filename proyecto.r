library(sf)
library(sp)
library(raster)
library(tidyverse)
library(tmap)
library(RColorBrewer)
library(lmtest)
library(spdep)
library(parallel)
library(ggplot2)
library(gstat)
library(stars)


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

#Adición de la cobertura 
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
munIncPercGeom <- munInc %>%st_centroid() %>%  mutate( 
  'IncPercentage' = munInc$NumIncendios/sum(munInc$NumIncendios)*100,
  'IncPercentage_log' = log1p(munInc$NumIncendios/sum(munInc$NumIncendios)*100),
  'IncPercentage_tukey' = rcompanion::transformTukey(NumIncendios/sum(munInc$NumIncendios)*100, plotit = F),
  'IncPercentage_tukey_lambda' = rcompanion::transformTukey(NumIncendios/sum(munInc$NumIncendios)*100, returnLambda = T),
  'AreaKm2' = as.numeric((st_area(munInc)/1000000)),
  'IncXArea' = (munInc$NumIncendios/AreaKm2),
  'IncXArea_log' = log1p(munInc$NumIncendios/AreaKm2),
  'IncXArea_tukey' = rcompanion::transformTukey(NumIncendios/AreaKm2, plotit = F),
  'IncXArea_tukey_lambda' = rcompanion::transformTukey(NumIncendios/AreaKm2, returnLambda = T),
  x=unlist(map(geom,1)), y=unlist(map(geom,2))) 

munIncPerc <-  munIncPercGeom %>%st_drop_geometry() 

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
# * Add month and Year field
incendiosForestales$month <- format(as.Date(incendiosForestales$ACQ_DATE), "%m")
incendiosForestales$year <- strtoi(format(as.Date(incendiosForestales$ACQ_DATE), "%Y"))
hist(incendiosForestales$year)
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

puntos_calor_fuegos_y_worldclim_viento_temp_radiacion


########### EASYMODE
#extracción de datos raster
st_write(incendiosForestales, "data/Incendios/incendiosForestales.shp")

tablaWorldClim <- read.table("data/Incendios/Test2.csv", 
                 header = TRUE,
                 sep = ",")

###########END EASYMODE


incForWorldClim <- left_join(incendiosForestales, select(tablaWorldClim, c(unique,month,wind,temp,srad,prec)), by = 'unique')

#quitar valores 0
incForestWorldClim <- incForWorldClim[incForWorldClim$srad != 0 & incForWorldClim$prec != 0, ] 

#Segmentacion por año
incForestWC2017 <- incForWorldClim[incForWorldClim$year == 2017, ] 
incForWC2017 <- incForestWC2017[incForestWC2018$srad != 0 & incForestWC2018$prec != 0, ] 
incForestWC2016 <- incForWorldClim[incForWorldClim$year == 2016, ] 
incForWC2016 <- incForestWC2016[incForestWC2018$srad != 0 & incForestWC2018$prec != 0, ] 
incForestWC2018 <- incForWorldClim[incForWorldClim$year == 2018, ] 
incForWC2018 <- incForestWC2018[incForestWC2018$srad != 0 & incForestWC2018$prec != 0, ] 


#Año 2018

#EDA
nrow(incForWC2018)
summary(incForWC2018$prec)


#Temperatura
hist(incForWC2018$prec)
shapiro.test(incForWC2018$temp)
shapiro.test(log(incForWC2018$temp))
ggplot() +
  geom_sf(data = mun4326, fill = 'white') +
  geom_sf(data = incForWC2018, aes(col = temp), size = 3) +
  scale_colour_gradient(low="#f7dede", high="#ff0011") +
  geom_sf_text(data = incForWC2018, aes(label=MUN), check_overlap = T, size = 0.5) +
  theme_bw()

#Radiacion Solar
shapiro.test(incForWC2018$srad)
shapiro.test(log(incForWC2018$srad))
ggplot() +
  geom_sf(data = mun4326, fill = 'white') +
  geom_sf(data = incForWC2018, aes(col = srad), size = 3) +
  scale_colour_gradient(low="#f7dede", high="#ff0011") +
  geom_sf_text(data = incForWC2018, aes(label=MUN), check_overlap = T, size = 0.5) +
  theme_bw()

#Viento
shapiro.test(incForWC2018$wind)
shapiro.test(log(incForWC2018$wind))
ggplot() +
  geom_sf(data = mun4326, fill = 'white') +
  geom_sf(data = incForWC2018, aes(col = wind), size = 3) +
  scale_colour_gradient(low="#deebf7", high="#3182bd") +
  geom_sf_text(data = incForWC2018, aes(label=MUN), check_overlap = T, size = 0.5) +
  theme_bw()

#Precipitacion
shapiro.test(incForWC2018$prec)
shapiro.test(log(incForWC2018$prec))
ggplot() +
  geom_sf(data = mun4326, fill = 'white') +
  geom_sf(data = incForWC2018, aes(col = prec), size = 3) +
  scale_colour_gradient(low="#deebf7", high="#3182bd") +
  geom_sf_text(data = incForWC2018, aes(label=MUN), check_overlap = T, size = 0.5) +
  theme_bw()


#variograma 
#Radiacion Solar
vsrad18 <- variogram(srad~1, incForWC2018)
plot(vsrad18, plot.numbers = T)
vsrad18_m1 <- fit.variogram(vsrad18, vgm(model = "Exp", range = 5))
vsrad18_m1
plot(vsrad18, vsrad18_m1, plot.numbers = T)

#Temperatura
vtemp18 <- variogram(temp~1, incForWC2018)
plot(vtemp18, plot.numbers = T)
vtemp18_m1 <- fit.variogram(vtemp18, vgm(model = "Exp", range = 50))
vtemp18_m1
plot(vtemp18, vtemp18_m1, plot.numbers = T)
attr(vtemp18_m1, 'SSErr')

#Viento
vwind18 <- variogram(wind~1, incForWC2018)
plot(vwind18, plot.numbers = T)
vwind18_m1 <- fit.variogram(vwind18, vgm(model = "Lin", range = 50))
vwind18_m1
plot(vwind18, vwind18_m1, plot.numbers = T)
vwind18_m2 <- fit.variogram(vwind18, vgm(model = "Exp", range = 50))
vwind18_m2
plot(vwind18, vwind18_m2, plot.numbers = T)
vwind18_m3 <- fit.variogram(vwind18, vgm(model = "Pow", range = 1))
vwind18_m3
plot(vwind18, vwind18_m3, plot.numbers = T)

attr(vwind18_m1, 'SSErr')
attr(vwind18_m2, 'SSErr')
attr(vwind18_m3, 'SSErr') #***

#Precipitacion
vprec18 <- variogram(prec~1, incForWC2018)
plot(vprec18, plot.numbers = T)
vprec18_m1 <- fit.variogram(vprec18, vgm(model = "Pen", range = 5000))
vprec18_m1
plot(vprec18, vprec18_m1, plot.numbers = T)

var1 <- fit.variogram()

## Loading required package: abind
grd <- st_bbox(mun) %>%
  st_as_stars(dx = 10000) %>% #10000 metros=10km de resolución espacial
  st_crop(mun)
grd
grd4326 <- st_transform(grd, crs = 4326)
plot(grd4326)

incForWC2018.32619 <- st_transform(incForWC2018, crs = 32619)

#krigging
k <- krige(formula = wind~1, locations = incForWC2018.32619, newdata = grd, model = vwind18_m3)
plot(k)


ggplot() +
  geom_stars(data = k, aes(fill = var1.pred, x = x, y = y)) + 
  scale_fill_gradient(low="#deebf7", high="#3182bd") +
  geom_sf(data = st_cast(mun, "MULTILINESTRING")) +
  geom_sf(data = incForWC2018) +
  geom_sf_text(data = mun, aes(label=TOPONIMIA), check_overlap = T, size = 0.8) +
  theme_bw()


#DEM

nrow(munIncPerc)
summary(munIncPerc$NumIncendios)
hist(munIncPerc$NumIncendios)
qqnorm(munIncPerc$NumIncendios)

hist(log(munIncPerc$NumIncendios))
qqnorm(log(munIncPerc$NumIncendios))
shapiro.test(munIncPerc$NumIncendios)
shapiro.test(log(munIncPerc$NumIncendios))

ggplot() +
  geom_sf(data = mun4326, fill = 'white') +
  geom_sf(data = munIncPercGeom, aes(col = NumIncendios), size = 6) + 
  scale_colour_gradientn(colours = rev(brewer.pal(9, name = 'RdBu'))) +
  geom_sf_text(data = munIncPercGeom, aes(label=TOPONIMIA), check_overlap = T, size = 2) +
  theme_bw()



#dem
dem <- read_stars('data/DTM/dem_srtm_remuestreado.tif')
names(dem) <- 'ele'
dem4326 <- st_transform(dem, crs = 4326)
plot(dem4326)

grdcovars <- aggregate(dem, grd, mean, na.rm=T) %>% st_transform(crs = 4326)
plot(grdcovars)

munIncPercGeom$ele <- st_as_sf(aggregate(grdcovars, munIncPercGeom, mean))[[1]]
munIncPercGeom <- munIncPercGeom [!is.na(munIncPercGeom$ele),]
plot(munIncPercGeom$NumIncendios, munIncPercGeom$ele)

munIncPercGeom_lm <- lm(NumIncendios ~ ele, munIncPercGeom)
summary(munIncPercGeom_lm)
plot(munIncPercGeom_lm)

vt <- variogram(NumIncendios ~ ele, munIncPercGeom)
vt
plot(vt)

vt_m <- fit.variogram(vt, vgm(model = "Exp", range = 50))
vt_m
plot(vt, vt_m, plot.numbers = T)

vt_m2 <- fit.variogram(vt, vgm(model = "Pen", range = 50))
vt_m2
plot(vt, vt_m2, plot.numbers = T)

vt_m3 <- fit.variogram(vt, vgm(model = "Sph", range = 50))
vt_m3
plot(vt, vt_m3, plot.numbers = T)

k_u <- krige(NumIncendios ~ ele , munIncPercGeom, st_rasterize(st_as_sf(grdcovars)), vt_m)


ggplot() +
  geom_stars(data = k_u, aes(fill = var1.pred, x = x, y = y)) + 
  scale_fill_gradientn(colours = rev(brewer.pal(9, name = 'RdBu'))) +
  geom_sf(data = st_cast(mun4326, "MULTILINESTRING")) +
  geom_sf(data = munIncPercGeom) +
  geom_sf_text(data = mun4326, aes(label=TOPONIMIA), check_overlap = T, size = 1) +
  theme_bw()