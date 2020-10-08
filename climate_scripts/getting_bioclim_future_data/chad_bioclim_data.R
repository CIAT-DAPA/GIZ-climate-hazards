# Chad data for Ani
# By: H. Achicanoy
# Alliance Bioversity-CIAT, 2020

library(raster)
library(tidyverse)

chad <- raster::getData(name = 'GADM', country = 'TCD', level = 1, path = 'C:/Users/haachicanoy/Downloads/CMIP5')

root <- 'C:/Users/haachicanoy/Downloads/CMIP5'
fldr <- list.dirs(root, recursive = F)

fldr %>%
  purrr::map(.f = function(fld){
    stck <- list.files(fld, full.names = T) %>%
      raster::stack() %>%
      raster::crop(x = ., y = raster::extent(chad)) %>%
      raster::mask(x = ., mask = chad)
    nms <- c('bio_01.tif','bio_12.tif')
    lsts <- stck %>% raster::unstack()
    1:length(lsts) %>%
      purrr::map(.f = function(i){
        raster::writeRaster(lsts[[i]], paste0(fld,'/',nms[i]))
      })
  })
