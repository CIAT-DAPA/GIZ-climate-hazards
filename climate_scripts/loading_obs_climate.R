# Load climate files per county
# By: H. Achicanoy & A. Esquivel
# CIAT, 2020

options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, vroom, sp, fst))

# country <- 'Ethiopia'
# county  <- 'Arsi'

country <- 'Mali'
iso3c   <- 'MLI'
adm_lvl <- 1
counties <- c('Koulikoro','Sikasso','Kayes','Segou','Mopti')
county <- counties[1]

# Paths
root <- '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles'

# Scripts
source(paste0(root,'/scripts/win_parallelization.R'))

# Load county shapefile
shp <- raster::shapefile(paste0(root,'/data/shps/',country,'/',iso3c,'_adm',adm_lvl,'.shp'))
glue::glue('shp <- shp[shp@data$NAME_{adm_lvl} == county,]') %>%
  as.character %>%
  parse(text = .) %>%
  eval(expr = ., envir = .GlobalEnv)

# Load id coords
crd <- vroom(paste0(root,'/data/id_country.csv'), delim = ',')
crd <- crd %>%
  dplyr::filter(Country == country)
pnt <- crd %>% dplyr::select(x,y) %>% sp::SpatialPoints(coords = .)
crs(pnt) <- crs(shp)
# Filter coordinates that are present in the county
pnt <- sp::over(pnt, shp) %>% data.frame %>% dplyr::select(ISO) %>% complete.cases() %>% which()
crd <- crd[pnt,]
crd <<- crd

# His obs
if(!file.exists(paste0(root,'/data/observational_data/',tolower(country),'/',tolower(county),'_prec_temp.RDS'))){
  ts <- list.files(path=paste0(root,'/data/Chirps_Chirts'),full.names=F,pattern='*.fst$') %>% sort()
  cl <- createCluster(30, export = list("ts","root","crd"), lib = list("tidyverse","fst"))
  # Temperatures, precipitation, and solar radiation data
  temp_prec <- ts %>% parallel::parLapply(cl, ., function(i){
    df <- fst::read_fst(paste0(root,'/data/Chirps_Chirts/',i))
    df <- df[df$id %in% crd$id,]
    return(df)
  }) %>%
    do.call(rbind, .)
  parallel::stopCluster(cl)
  temp_prec2 <- temp_prec %>% dplyr::group_by(id) %>% dplyr::arrange(Date) %>% dplyr::group_split(id)
  
  his_obs <- tibble::tibble(id      = crd$id,
                            x       = crd$x,
                            y       = crd$y,
                            ISO3    = crd$ISO3,
                            Country = crd$Country,
                            Climate = temp_prec2 %>% purrr::map(function(tb){tbl <- tb %>% dplyr::select(id,Date,prec,tmax,tmin); return(tbl)}))
  outDir <- paste0(root,'/data/observational_data/',tolower(country))
  if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
  out <- paste0(outDir,'/',tolower(county),'_prec_temp.RDS')
  saveRDS(his_obs, out)
} else {
  out      <- paste0(root,'/data/observational_data/',tolower(country),'/',tolower(county),'_prec_temp.RDS')
  his_obs  <- readRDS(out)
  if(!file.exists(paste0(root,'/data/observational_data/',tolower(country),'/',tolower(county),'.RDS'))){
    srad_fls <- paste0(root,"/data/NASA/",crd$id,'.fst')
    if(sum(file.exists(srad_fls)) == nrow(crd)){
      srad <- srad_fls %>% purrr::map(.f = function(x){
        df <- fst::read_fst(x)
      })
      all_clim <- purrr::map2(.x = his_obs$Climate, .y = srad, .f = function(x, y){
        z <- dplyr::left_join(x = x, y = y %>% dplyr::select(id, Date, srad), by = c('id','Date'))
        return(z)
      })
      his_obs$Climate <- all_clim
      outDir <- paste0(root,'/data/observational_data/',tolower(country))
      if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
      out <- paste0(outDir,'/',tolower(county),'.RDS')
      saveRDS(his_obs, out)
    } else {
      cat('SRAD files are not available\n')
    }
  } else {
    cat('File exists\n')
  }
  
}
