# Create AEZ regions
# By: H. Achicanoy
# Alliance Bioversity-CIAT, 2020

library(tidyverse)
library(raster)
library(sf)

bio12 <- raster("//dapadfs.cgiarad.org/data_cluster_4/observed/gridded_products/worldclim/Global_30s_v2/bio_12.tif")
shpjd <- raster::shapefile("//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data/shps/Jordan/Jordan_wgs84_pol.shp")
shpjd$ID <- rownames(shpjd@data)
jordan_valley <- shpjd[shpjd@data$ID == 2,]
crs_obj <- crs(shpjd)

bio12 <- bio12 %>% raster::crop(x = ., y = raster::extent(shpjd)) %>% raster::mask(x = ., mask = shpjd)

# Desert (annual prec < 150 mm)
desert <- raster::rasterToPolygons(bio12, fun = function(x){x < 150})
desert <- as(desert, "sf")
desert <- sf::st_union(desert)
desert <- as(desert, "Spatial")
crs(desert) <- crs_obj
desert <- rgeos::gDifference(desert, jordan_valley)
raster::shapefile(desert, '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/Jordan_desert.shp')

# Agro-pastoral (annual prec >= 150 mm & < 250 mm)
agrpas <- raster::rasterToPolygons(bio12, fun = function(x){x >= 150 & x < 250})
agrpas <- as(agrpas, "sf")
agrpas <- sf::st_union(agrpas)
agrpas <- as(agrpas, "Spatial")
crs(agrpas) <- crs_obj
agrpas <- rgeos::gDifference(agrpas, jordan_valley)

raster::shapefile(agrpas, '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/Jordan_agro-pastoral.shp')

# Rainfed (annual prec > 250 mm)
rainfd <- raster::rasterToPolygons(bio12, fun = function(x){x >= 250})
rainfd <- as(rainfd, "sf")
rainfd <- sf::st_union(rainfd)
rainfd <- as(rainfd, "Spatial")
crs(rainfd) <- crs_obj
rainfd <- rgeos::gDifference(rainfd, jordan_valley)

raster::shapefile(rainfd, '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/Jordan_rainfed.shp')

jordan_valley@data
desert@data

jordan_shp <- raster::bind(jordan_valley,desert,agrpas,rainfd)
jordan_shp@data$SHAPE_Leng <- NULL
jordan_shp@data$SHAPE_Area <- NULL
names(jordan_shp@data)[1] <- 'Region'
jordan_shp@data$Region <- c('Irrigated','Desert','Agro-pastoral','Rainfed')
jordan_shp@data$ID <- 1:4
jordan_shp@data$Description <- c('Annual prec >= 250 mm (mostly Jordan Valley)', 'Annual prec < 150 mm', 'Annual prec >= 150 mm & < 250 mm', 'Annual prec >= 250 mm')

raster::shapefile(jordan_shp, '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/Jordan_AEZ_defined.shp')
