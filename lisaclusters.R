lisamap <- function(objesp = NULL, var = 'mivariable_pct_log', pesos = NULL,
                    titulomapa = NULL, anchuratitulo = NULL, tamanotitulo = NULL,
                    tituloleyenda = NULL, leyenda = T, fuentedatos = NULL) {
  require(tidyverse)
  require(spdep)
  require(sf)
  require(ggplot2)
  
  # Variable en forma vectorial
  varvectorial <- objesp[var] %>% st_drop_geometry %>% dplyr::select(var) %>% pull
  
  # Moral local
  lomo <- localmoran(varvectorial, listw = pesos)
  
  # Puntuaciones z
  objesp$puntuacionz <- varvectorial %>% scale %>% as.vector
  
  # Crear variable con rezago
  objesp$lagpuntuacionz <- lag.listw(pesos, objesp$puntuacionz)
  
  # Variable nueva sobre significancia de la correlación local, rellena con NAs
  objesp$quad_sig <- NA
  
  # Cuadrante high-high quadrant
  objesp[(objesp$puntuacionz >= 0 & 
            objesp$lagpuntuacionz >= 0) & 
           (lomo[, 5] <= 0.05), "quad_sig"] <- "high-high"
  # Cuadrante low-low
  objesp[(objesp$puntuacionz <= 0 & 
            objesp$lagpuntuacionz <= 0) & 
           (lomo[, 5] <= 0.05), "quad_sig"] <- "low-low"
  # Cuadrante high-low
  objesp[(objesp$puntuacionz >= 0 & 
            objesp$lagpuntuacionz <= 0) & 
           (lomo[, 5] <= 0.05), "quad_sig"] <- "high-low"
  # Cuadrante low-high
  objesp[(objesp$puntuacionz <= 0 
          & objesp$lagpuntuacionz >= 0) & 
           (lomo[, 5] <= 0.05), "quad_sig"] <- "low-high"
  # No significativas
  objesp[(lomo[, 5] > 0.05), "quad_sig"] <- "not signif."  
  
  #Convertir a factorial
  objesp$quad_sig <- as.factor(objesp$quad_sig)
  
  # Mapa
  p <- objesp %>% 
    ggplot() +
    aes(fill = quad_sig) + 
    geom_sf(color = "black", size = .05)  +
    theme_bw() +
    scale_fill_manual(values = c("high-high" = "red", "low-low" = "blue", "high-low" = "lightblue", "low-high" = "pink", "not signif." = "light grey")) +
    labs(title = stringr::str_wrap(titulomapa, width = anchuratitulo),
         fill = tituloleyenda,
         caption = paste('Fuente de datos:', fuentedatos)) +
    scale_x_continuous(breaks = -68:-72) +
    scale_y_continuous(breaks = 18:20) +
    theme(plot.title = element_text(size = tamanotitulo), legend.position = {`if`(leyenda, 'right', 'none')}) 
  return(list(grafico = p, objeto = objesp))
}