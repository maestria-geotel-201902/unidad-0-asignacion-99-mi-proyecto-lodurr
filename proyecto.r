library(sf)
library(spdep)

mun <- st_read(dsn = 'data/DivisionRD/divisionRD.gpkg', layer = 'MUNCENSO2010')
prov <- st_read(dsn = 'data/DivisionRD/divisionRD.gpkg', layer = 'PROVCENSO2010')
mun4326 <- st_transform(mun, crs = 4326)
prov4326 <- st_transform(prov, crs = 4326)
colnames(mun4326)

incendios <- st_read(dsn = "data/Incendios/FiredataBuffer.shp")
incendiosMunRD <- st_intersection(incendios,mun4326)
incendiosProvRD <- st_intersection(incendios,prov4326)
plot(incendiosProvRD['ENLACE'])

munIncendios <- mun4326 %>% inner_join(incendiosMunRD, by = 'ENLACE')

match (incendiosMunRD$ENLACE, mun4326$ENLACE)

incendiosMunRD.sp <- as_Spatial(incendiosMunRD)
colnames(incendiosMunRD.sp@data)
incendiosProvRD.sp <- as_Spatial(incendiosProvRD) 
colnames(incendiosProvRD.sp@data)

incendiosMunRD.nb <- poly2nb(incendiosMunRD.sp, queen=TRUE)
