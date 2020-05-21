# Calculate all agro-climatic indices
# A. Esquivel and H. Achicanoy
# CIAT, 2020

options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, ncdf4, sf, future, furrr, lubridate, glue, vroom, sp, fst, compiler))

OSys <- Sys.info()[1]
root <<- switch(OSys,
                'Linux'   = '/home/jovyan/work/cglabs',
                'Windows' = '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles')

calc_indices <- function(country = 'Pakistan',
                         county  = 'Mithi',
                         iso3c   = 'PAK',
                         adm_lvl = 3,
                         gcm     = NULL,
                         period  = '1985_2015',
                         time    = 'past'){
  
  country <<- country
  county  <<- county
  
  # Load county shapefile
  shp <<- raster::shapefile(paste0(root,'/data/shps/',country,'/',iso3c,'_adm',adm_lvl,'.shp'))
  glue::glue('shp <<- shp[shp@data$NAME_{adm_lvl} == county,]') %>%
    as.character %>%
    parse(text = .) %>%
    eval(expr = ., envir = .GlobalEnv)
  
  # Load id coords
  crd <- vroom(paste0(root,'/data/id_all_country.csv'), delim = ',')
  crd <- crd %>%
    dplyr::filter(Country == country)
  pnt <- crd %>% dplyr::select(x,y) %>% sp::SpatialPoints(coords = .)
  crs(pnt) <- crs(shp)
  # Filter coordinates that are present in the county
  pnt <- sp::over(pnt, shp) %>% data.frame %>% dplyr::select(ISO) %>% complete.cases() %>% which()
  crd <- crd[pnt,]
  crd <<- crd
  
  timList <- c('past','future')
  periodList <- c('2021_2045','2041_2065')
  
  source(paste0(root,'/scripts/indices.R'))
  
  # Paths
  obsDir <- paste0(root,'/data/observational_data/',tolower(country))
  futDir <- paste0(root,'/data/bc_quantile_0_05deg_lat_county/',tolower(country),'/',gcm,'/',period)
  outDir <- ifelse(test = time == 'past',
                   yes  = paste0(root,'/results/',country,'/',time),
                   no   = paste0(root,'/results/',country,'/',time,'/',gcm,'/',period))
  
  out <- paste0(outDir,'/',county,'_',period,'_all_new.fst')
  if(!file.exists(out)){
    # Soil data
    Soil <- fst::read.fst(paste0(root,'/data/soilcp_data.fst')) %>%
      tibble::as_tibble() %>%
      dplyr::select(id, soilcp) %>%
      dplyr::filter(id %in% dplyr::pull(crd, id))
    
    # Read climate data
    if(time == 'past'){
      if(!file.exists(paste0(obsDir,'/',tolower(county),'.RDS'))){
        clim_data <- readRDS(paste0(obsDir,'/',tolower(county),'_prec_temp.RDS'))
      } else {
        clim_data <- readRDS(paste0(obsDir,'/',tolower(county),'.RDS'))
      }
    } else {
      if(!file.exists(paste0(futDir,'/',tolower(county),'.RDS'))){
        clim_data <- readRDS(paste0(futDir,'/',tolower(county),'_prec_temp.RDS'))
      } else {
        clim_data <- readRDS(paste0(futDir,'/',tolower(county),'.RDS'))
      }
    }
    
    # Impute missing data
    impute_missings <- function(tbl = clim_data){
      Climate <- 1:nrow(tbl) %>%
        purrr::map(.f = function(i){
          df <- tbl$Climate[[i]]
          if(sum(is.na(df$tmax)) > 0){
            df$tmax[which(is.na(df$tmax))] <- median(df$tmax, na.rm = T)
          }
          if(sum(is.na(df$tmin)) > 0){
            df$tmin[which(is.na(df$tmin))] <- median(df$tmin, na.rm = T)
          }
          if(sum(is.na(df$srad)) > 0){
            df$srad[which(is.na(df$srad))] <- median(df$srad, na.rm = T)
          }
          if(sum(is.na(df$prec)) > 0){
            df$prec[which(is.na(df$prec))] <- median(df$prec, na.rm = T)
          }
          return(df)
        })
      
      return(Climate)
    }
    clim_data$Climate <- impute_missings(tbl = clim_data)
    
    run_pixel <- function(id = 72516){
      
      cat(' --- Obtain complete time series per pixel\n')
      tbl <- clim_data$Climate[[which(clim_data$id == id)]]
      tbl <- tbl %>%
        dplyr::mutate(year  = lubridate::year(as.Date(Date)),
                      month = lubridate::month(as.Date(Date)))
      years <- tbl$year %>% unique
      
      cat(' --- Calculate water balance for complete time series\n')
      soilcp <- Soil$soilcp[Soil$id == id]
      watbal_loc <- watbal_wrapper(out_all = tbl, soilcp = soilcp)
      watbal_loc$IRR <- watbal_loc$Etmax - watbal_loc$prec
      
      tbl <- tbl %>%
        dplyr::mutate(ERATIO = watbal_loc$ERATIO,
                      IRR = watbal_loc$IRR)
      
      cat(' --- Subsetting by Kafir and Rabi seasons\n')
      tbl_kafir <- tbl %>%
        dplyr::filter(month %in% 5:11)
      pairs     <- NA; for(i in 1:length(years)-1){pairs[i] <- paste0(years[i:(i+1)], collapse = '-')}
      tbl_list  <- lapply(1:(length(years)-1), function(i){
        df <- tbl %>%
          dplyr::filter(year %in% years[i:(i+1)])
        df$pairs <- paste0(years[i:(i+1)], collapse = '-')
        df1 <- df %>%
          dplyr::filter(year == years[i] & month %in% 11:12)
        df2 <- df %>%
          dplyr::filter(year == years[i+1] & month %in% 1:4)
        df <- rbind(df1, df2); rm(df1, df2)
        return(df)
      })
      tbl_rabi  <- dplyr::bind_rows(tbl_list); rm(tbl_list, pairs)
      
      idx_kafir <- tbl_kafir %>%
        dplyr::group_split(year) %>%
        purrr::map(.f = function(df){
          idx <- tibble::tibble(CDD  = calc_cddCMP(PREC = df$prec),
                                P5D  = calc_p5dCMP(PREC = df$prec),
                                P95  = calc_p95CMP(PREC = df$prec),
                                NT35 = calc_htsCMP(tmax = df$tmax, t_thresh = 35),
                                ndws = calc_wsdays(df$ERATIO, season_ini=1, season_end=length(df$ERATIO), e_thresh=0.5),
                                IRR  = sum(df$IRR, na.rm = T))
          return(idx)
        })
      idx_kafir <- dplyr::bind_rows(idx_kafir)
      idx_kafir$year <- tbl_kafir$year %>% unique
      idx_kafir$season <- 'Kafir'
      
      idx_rabi <- tbl_rabi %>%
        dplyr::group_split(pairs) %>%
        purrr::map(.f = function(df){
          idx <- tibble::tibble(CDD  = calc_cddCMP(PREC = df$prec),
                                P5D  = calc_p5dCMP(PREC = df$prec),
                                P95  = calc_p95CMP(PREC = df$prec),
                                NT35 = calc_htsCMP(tmax = df$tmax, t_thresh = 35),
                                ndws = calc_wsdays(df$ERATIO, season_ini=1, season_end=length(df$ERATIO), e_thresh=0.5),
                                IRR  = sum(df$IRR, na.rm = T))
          return(idx)
        })
      idx_rabi <- dplyr::bind_rows(idx_rabi)
      idx_rabi$year <- tbl_rabi$pairs %>% unique %>% substring(text = ., first = 1, last = 4) %>% as.numeric
      idx_rabi$season <- 'Rabi'
      
      indices <- rbind(idx_kafir, idx_rabi)
      indices$id <- id
      
      return(indices)
      
    }
    plan(multiprocess)
    index_by_pixel <- clim_data %>%
      dplyr::pull(id) %>% 
      furrr::future_map(.x = ., .f = run_pixel) %>% 
      dplyr::bind_rows()
    gc()
    gc(reset = T)
    
    index_by_pixel$Country <- country
    index_by_pixel$county <- county
    index_by_pixel$ISO3 <- iso3c
    index_by_pixel <- dplyr::left_join(x = index_by_pixel, y = clim_data %>% dplyr::select(id, x, y) %>% unique, by = 'id')
    index_by_pixel <- index_by_pixel %>% dplyr::select(id, x, y, ISO3, Country, county, year, season, dplyr::everything(.))
    
    if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
    fst::write_fst(x = index_by_pixel, path = out)
    cat('>>> File created successfully ...\n')
  } else {
    cat('>>> File exists it is not necessary to create it again\n')
  }
  
}

countyList <- c('Muzaffargarh',
                'Rajan Pur',
                'Jhang',
                'Ghotki',
                'Kashmore',
                'Dadu',
                'Chitral',
                'Dera Ismail Khan',
                'South Waziristan',
                'North Waziristan',
                'Orakzai',
                'Kurram')
for(i in 1:length(countyList)){
  calc_indices(country = 'Pakistan',
               county  = countyList[i],
               iso3c   = 'PAK',
               adm_lvl = 3,
               gcm     = NULL,
               period  = '1985_2015',
               time    = 'past')
}

gcmList <- c("ipsl_cm5a_mr","miroc_esm_chem","ncc_noresm1_m")
periodList <- c('2021_2045','2041_2065')
for(i in 1:length(countyList)){
  for(gcm in gcmList){
    for(period in periodList){
      calc_indices(country = 'Pakistan',
                   county  = countyList[i],
                   iso3c   = 'PAK',
                   adm_lvl = 3,
                   gcm     = gcm,
                   period  = period,
                   time    = 'future')
    }
  }
}
