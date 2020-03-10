# Soil capacity calculation: GIZ climate-hazards
# By: H. Achicanoy
# CIAT, 2020

options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
pacman::p_load(raster, tidyverse, fst, GSIF)

root_depth <- 60 # cm

coords <- fst::read_fst('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data/Chirps_Chirts/D_1985.01.01.fst')
coords <- coords %>% dplyr::select(id,x,y)

soils_root <- '//catalogue/BaseLineData_cluster04/GLOBAL/Biofisico/SoilGrids250m'
orc <- raster::stack(list.files(paste0(soils_root,'/Chemical soil properties/Soil organic carbon content'), pattern = '.tif$', full.names = T) %>% sort())
cec <- raster::stack(list.files(paste0(soils_root,'/Chemical soil properties/Cation exchange capacity (CEC)'), pattern = '.tif$', full.names = T) %>% sort())
phx <- raster::stack(list.files(paste0(soils_root,'/Chemical soil properties/Soil ph in H2O'), pattern = '.tif$', full.names = T) %>% sort())
snd <- raster::stack(list.files(paste0(soils_root,'/Physical soil properties/Sand content'), pattern = '.tif$', full.names = T) %>% sort())
slt <- raster::stack(list.files(paste0(soils_root,'/Physical soil properties/Silt content'), pattern = '.tif$', full.names = T) %>% sort())
cly <- raster::stack(list.files(paste0(soils_root,'/Physical soil properties/Clay content (0-2 micro meter) mass fraction'), pattern = '.tif$', full.names = T) %>% sort())
bld <- raster::stack(list.files(paste0(soils_root,'/Physical soil properties/Bulk density (fine earth)'), pattern = '.tif$', full.names = T) %>% sort())

soil <- raster::stack(orc,cec,phx,snd,slt,cly,bld)
soil_data <- cbind(coords, raster::extract(soil, coords[,c('x','y')]))

soil_data2 <- soil_data %>%
  tidyr::gather(key = 'var', value = 'val', -(1:3)) %>%
  tidyr::separate(col = 'var', sep = '_M_', into = c('var','depth')) %>%
  tidyr::spread(key = 'var', value = 'val') %>%
  dplyr::arrange(id)
soil_data2$depth <- gsub('_250m_ll','',soil_data2$depth)

fst::write_fst(soil_data2, '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data/soil_data.fst')

soil_data2$d <- lapply(1:nrow(soil_data2), function(i){
  y <- GSIF::AWCPTF(SNDPPT = soil_data2$SNDPPT[i],
                    SLTPPT = soil_data2$SLTPPT[i],
                    CLYPPT = soil_data2$CLYPPT[i],
                    ORCDRC = soil_data2$ORCDRC[i],
                    BLD = soil_data2$BLDFIE[i],
                    CEC = soil_data2$CECSOL[i],
                    PHIHOX = soil_data2$PHIHOX[i])
  y <- y$AWCh2 * 100
  return(y)
}) %>% base::unlist()

soil_data3 <- soil_data2 %>%
  dplyr::select(id,x,y,depth,d) %>%
  tidyr::spread(key='depth',value='d')

names(soil_data3)[4:ncol(soil_data3)] <- paste0('d.',c(0, 5, 15, 30, 60, 100, 200))
soil_data3$rdepth <- root_depth
soil_data3 <- soil_data3 %>%
  dplyr::select('id','x','y','rdepth',paste0('d.',c(0, 5, 15, 30, 60, 100, 200)))

soilcap_calc_mod <- function(x, minval, maxval) {
  if(!is.na(x[4])){
    rdepth <- max(c(x[4],minval)) #cross check
    rdepth <- min(c(rdepth,maxval)) #cross-check
    wc_df <- data.frame(depth=c(2.5,10,22.5,45,80,150),wc=(x[5:10])*.01)
    if (!rdepth %in% wc_df$depth) {
      wc_df1 <- wc_df[which(wc_df$depth < rdepth),]
      wc_df2 <- wc_df[which(wc_df$depth > rdepth),]
      y1 <- wc_df1$wc[nrow(wc_df1)]; y2 <- wc_df2$wc[1]
      x1 <- wc_df1$depth[nrow(wc_df1)]; x2 <- wc_df2$depth[1]
      ya <- (rdepth-x1) / (x2-x1) * (y2-y1) + y1
      wc_df <- rbind(wc_df1,data.frame(depth=rdepth,wc=ya),wc_df2)
    }
    wc_df <- wc_df[which(wc_df$depth <= rdepth),]
    wc_df$soilthick <- wc_df$depth - c(0,wc_df$depth[1:(nrow(wc_df)-1)])
    wc_df$soilcap <- wc_df$soilthick * wc_df$wc
    soilcp <- sum(wc_df$soilcap) * 10 #in mm
    return(soilcp)
  } else {
    soilcp <- NA
    return(soilcp)
  }
}

# Calculate soil water holding capacity in mm, minval and maxval taken from
# Fatondji et al. (2012) --in: Kihara, J. et al. Improving soil fert. recommendation using DSSAT
soil_data$soilcp <- apply(soil_data, 1, FUN=soilcap_calc_mod, minval=45, maxval=100)
save(soil_data, file='//dapadfs/workspace_cluster_8/Kenya_KACCAL/data/input_tables/soil_data.RData')

# # ===================================================================== #
# # AfSIS process
# # ===================================================================== #
# 
# # Load packages
# options(warn=-1); options(scipen = 999)
# suppressMessages(if(!require(raster)){install.packages('raster'); library(raster)} else {library(raster)})
# suppressMessages(if(!require(ncdf)){install.packages('ncdf'); library(ncdf)} else {library(ncdf)})
# suppressMessages(if(!require(ncdf4)){install.packages('ncdf4'); library(ncdf4)} else {library(ncdf4)})
# suppressMessages(if(!require(maptools)){install.packages('maptools'); library(maptools)} else {library(maptools)})
# suppressMessages(if(!require(ff)){install.packages('ff'); library(ff)} else {library(ff)})
# suppressMessages(if(!require(data.table)){install.packages('data.table'); library(data.table)} else {library(data.table)})
# suppressMessages(if(!require(miscTools)){install.packages('miscTools'); library(miscTools)} else {library(miscTools)})
# 
# 
# 
# # Soils directory
# isric_dir <- '//dapadfs/workspace_cluster_8/Kenya_KACCAL/data/ISRIC_soil'
# 
# # Rasterized Kenya shapefile
# rs_adm <- raster('//dapadfs/workspace_cluster_8/Kenya_KACCAL/data/Kenya_counties_rst/Kenya_base.tif')
# 
# # Kenya shapefile
# shp <- readShapeSpatial('//dapadfs/workspace_cluster_8/Kenya_KACCAL/data/Kenya_counties_shp/County.shp', proj4string=CRS("+proj=longlat +datum=WGS84"))
# 
# # Extract all coordinates different from NULL
# loc_xy <- as.data.frame(xyFromCell(rs_adm, which(!is.na(rs_adm[]))))
# names(loc_xy) <- c("lon","lat")
# 
# # Extract soil data, and calculate soil water holding capacity for selected sites
# # soilcap = sum(af_AWCh2__M_sd[1-x]_1km) * x
# root_depth <- raster(paste(isric_dir, "/af_ERDICM__M_1km.tif", sep=""))
# rs_prj <- projectRaster(rs_adm, crs=root_depth@crs)
# root_depth <- crop(root_depth, extent(rs_prj))
# root_depth <- projectRaster(root_depth, crs=shp@proj4string)
# 
# # Calculate average of root depth for each big grid cell
# rs_res <- resample(rs_adm, root_depth, method="ngb")
# rs_ids <- rs_adm; rs_ids[which(!is.na(rs_ids[]))] <- which(!is.na(rs_ids[]))
# rs_pts <- as.data.frame(xyFromCell(rs_res, which(!is.na(rs_res[]))))
# rs_pts$id_coarse <- extract(rs_ids, cbind(x=rs_pts$x, y=rs_pts$y))
# rs_pts$rdepth <- extract(root_depth, data.frame(x=rs_pts$x, y=rs_pts$y))
# rs_pts <- aggregate(rs_pts[,c("x","y","rdepth")], by=list(id_coarse=rs_pts$id_coarse), FUN=function(x){mean(x,na.rm=T)})
# rs_pts$x <- rs_pts$y <- NULL
# 
# # Put root_depth data on soil_data data.frame
# soil_data <- loc_xy
# soil_data$id_coarse <- extract(rs_ids, data.frame(x=soil_data$lon, y=soil_data$lat))
# soil_data <- merge(soil_data, rs_pts, by="id_coarse")
# rm(rs_pts)
# 
# # Extract soil water holding capacity data on soil_data data.frame
# depths <- c(25,100,225,450,800,1500)
# for (s_i in 1:6) {
#   #s_i <- 1
#   tdepth <- depths[s_i]
#   cat("...extracting depth=",tdepth*.1,"cm\n")
#   rs <- raster(paste(isric_dir,"/af_AWCh2__M_sd",s_i,"_1km.tif",sep=""))
#   rs <- crop(rs, extent(rs_prj))
#   rs <- projectRaster(rs, crs=shp@proj4string)
#   rs_res <- resample(rs_adm, rs, method="ngb")
#   rs_pts <- as.data.frame(xyFromCell(rs_res, which(!is.na(rs_res[]))))
#   rs_pts$id_coarse <- extract(rs_ids, cbind(x=rs_pts$x, y=rs_pts$y))
#   rs_pts$value <- extract(rs, data.frame(x=rs_pts$x, y=rs_pts$y))
#   rs_pts <- aggregate(rs_pts[,c("x","y","value")],by=list(id_coarse=rs_pts$id_coarse),FUN=function(x) {mean(x,na.rm=T)})
#   rs_pts$x <- rs_pts$y <- NULL
#   soil_data <- merge(soil_data, rs_pts, by="id_coarse")
#   names(soil_data)[ncol(soil_data)] <- paste("d.",tdepth,sep="")
#   rm(list=c("rs","rs_pts","rs_res"))
# }
# 
# soilcap_calc_mod <- function(x, minval, maxval) {
#   if(!is.na(x[4])){
#     rdepth <- max(c(x[4],minval)) #cross check
#     rdepth <- min(c(rdepth,maxval)) #cross-check
#     wc_df <- data.frame(depth=c(2.5,10,22.5,45,80,150),wc=(x[5:10])*.01)
#     if (!rdepth %in% wc_df$depth) {
#       wc_df1 <- wc_df[which(wc_df$depth < rdepth),]
#       wc_df2 <- wc_df[which(wc_df$depth > rdepth),]
#       y1 <- wc_df1$wc[nrow(wc_df1)]; y2 <- wc_df2$wc[1]
#       x1 <- wc_df1$depth[nrow(wc_df1)]; x2 <- wc_df2$depth[1]
#       ya <- (rdepth-x1) / (x2-x1) * (y2-y1) + y1
#       wc_df <- rbind(wc_df1,data.frame(depth=rdepth,wc=ya),wc_df2)
#     }
#     wc_df <- wc_df[which(wc_df$depth <= rdepth),]
#     wc_df$soilthick <- wc_df$depth - c(0,wc_df$depth[1:(nrow(wc_df)-1)])
#     wc_df$soilcap <- wc_df$soilthick * wc_df$wc
#     soilcp <- sum(wc_df$soilcap) * 10 #in mm
#     return(soilcp)
#   } else {
#     soilcp <- NA
#     return(soilcp)
#   }
# }
# 
# # Calculate soil water holding capacity in mm, minval and maxval taken from
# # Fatondji et al. (2012) --in: Kihara, J. et al. Improving soil fert. recommendation using DSSAT
# soil_data$soilcp <- apply(soil_data, 1, FUN=soilcap_calc_mod, minval=45, maxval=100)
# save(soil_data, file='//dapadfs/workspace_cluster_8/Kenya_KACCAL/data/input_tables/soil_data.RData')
# 
# # Create soil_data table for each county
# lapply(1:nrow(countyList), function(i)
# {
#   cat('Load soil properties in Kenya\n')
#   load('//dapadfs/workspace_cluster_8/Kenya_KACCAL/data/input_tables/soil_data.RData')
#   
#   cat('Load shapefile for:', countyList$County[[i]], '\n')
#   shp_ras <- raster(paste('//dapadfs/workspace_cluster_8/Kenya_KACCAL/data/Kenya_counties_rst/', gsub(pattern=' ', replacement='_', countyList$County[[i]]), '_base.tif', sep=''))
#   shp_tab <- as.data.frame(xyFromCell(object=shp_ras, cell=1:ncell(shp_ras)))
#   shp_tab <- cbind(cellID=1:ncell(shp_ras), shp_tab)
#   shp_tab <- shp_tab[which(!is.na(shp_ras[])),]
#   names(shp_tab)[2:3] <- c('lon', 'lat')
#   
#   soil_data$cellID <- cellFromXY(object=shp_ras, xy=as.matrix(soil_data[, c('lon','lat')]))
#   soil_data_county <- merge(x=shp_tab, y=soil_data, by='cellID')
#   
#   counDir <- paste('//dapadfs/workspace_cluster_8/Kenya_KACCAL/data/input_tables/', gsub(pattern=' ', replacement='_', countyList$County[[i]]), '/soil', sep='')
#   if(!dir.exists(counDir)){ dir.create(counDir, recursive = T) } else { cat('County folder exists\n') }
#   save(soil_data_county, file=paste(counDir, '/soil_data.RData', sep=''))
#   
#   return(cat('Process has been done for:', countyList$County[[i]], '\n'))
#   
# })