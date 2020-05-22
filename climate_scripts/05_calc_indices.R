# Calculate agro-climatic partial indices
# A. Esquivel and H. Achicanoy
# CIAT, 2020

options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, ncdf4, sf, future, furrr, lubridate, glue, vroom, sp, fst, compiler))

OSys <- Sys.info()[1]
root <<- switch(OSys,
                'Linux'   = '/home/jovyan/work/cglabs',
                'Windows' = '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles')

calc_indices <- function(country = 'Malawi',
                         county  = 'Central Region',
                         iso3c   = 'MWI',
                         adm_lvl = 1,
                         seasons = 1,
                         gcm     = 'ipsl_cm5a_mr',
                         period  = '2041_2065',
                         time    = 'future'){
  
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
  # Soil data
  Soil <- fst::read.fst(paste0(root,'/data/soilcp_data.fst')) %>%
    tibble::as_tibble() %>%
    dplyr::select(id, soilcp) %>%
    dplyr::filter(id %in% dplyr::pull(crd, id))
  
  # Main functions
  run_each_semester <<- function(one, semester){
    
    one2 <- one %>%
      dplyr::select(id, year) %>%
      unique()
    
    
    # one2 <- one2 %>% dplyr::mutate(season = list(rsum.lapply(x = one$prec %>% as.numeric, n = 150)))
    
    one2 <- one2 %>%
      dplyr::mutate(season = purrr::map(.x =year, .f = function(y){
        oi = list(rsum.lapply(x = filter(one, year == y)$prec %>% as.numeric, n = 150))[[1]] } ) )
    
    # one2 <- mutate(one2, sum.ind = list(season[[1]][[which.max(cumulative.r.sum(season[[1]]))]]))
    
    one2 <- one2 %>%
      mutate(sum.ind = purrr::map(.x = season, function(y){list(y[[which.max(cumulative.r.sum(y[[1]]))]])[[1]]}))
    
    # data <-  one2 %>%
    #   mutate(semester = semester) %>%
    #   dplyr::mutate(CDD = calc_cddCMP(one$prec[one2$sum.ind[[1]][[2]]])) %>%
    #   dplyr::mutate(P5D = calc_p5dCMP(one$prec[one2$sum.ind[[1]][[2]]])) %>%
    #   dplyr::mutate(P95 = calc_p95CMP(one$prec[one2$sum.ind[[1]][[2]]])) %>%
    #   dplyr::mutate(NT35 = calc_htsCMP(one$tmax[one2$sum.ind[[1]][[2]]], t_thresh = 35)) %>%
    #   dplyr::select(-season) %>%
    #   dplyr::select(id, year, semester, everything())
    
    data <- one2 %>% 
      mutate(semester = semester) %>% 
      mutate(climate = purrr::map(.x = year, .f = function(x){filter(one , year == x)})) %>%
      mutate(index = purrr::map2(.x = sum.ind, .y = climate, .f = function(x, y){
        tibble(CDD = calc_cddCMP(y$prec[x[[2]]]), 
               P5D = calc_p5dCMP(y$prec[x[[2]]]), 
               P95 = calc_p95CMP(y$prec[x[[2]]]), 
               NT35 = calc_htsCMP(y$tmax[x[[2]]], t_thresh = 35))
      })) %>% unnest(index)
    
    
    return(data)}
  reading_run_pys   <<- function(px){
    
    id <- px$id ; semester <- px$semester ; ISO3 <- px$ISO3 ; county <- px$county; soilcp <- px$soilcp
    
    if(semester == 2){
      # Si se tienen 2 se tiene un semestres...
      # semester = 2
      data_s <-  dplyr::select(px, Climate)  %>% unnest() %>% 
        mutate(year = lubridate::year(Date), month = lubridate::month(Date), 
               semester = ifelse(month < 7, 1, 2)) %>% 
        dplyr::select(-month) 
      
      one <- data_s %>% 
        dplyr::select(id, semester,  year, prec, tmax, tmin) %>%
        nest(-semester)
      
      data_base <- purrr::map2(.x = one$data, .y = one$semester, .f = run_each_semester) %>% bind_rows()
      
    }else if(semester == 1){
      
      one  <- dplyr::select(px, Climate)  %>% unnest() %>% 
        mutate(year = lubridate::year(Date)) %>%
        dplyr::select(id, year, prec, tmax, tmin)
      
      
      data_base <- run_each_semester(one = one, semester =  1)
      
    }else{ data_base <- NULL}
    
    return(data_base)}
  
  if(time == 'past'){
    if(!file.exists(paste0(obsDir,'/',tolower(county),'.fst'))){
      clim_data <- fst::read_fst(paste0(obsDir,'/',tolower(county),'_prec_temp.fst'))
      clim_data <- clim_data %>%
        tidyr::nest(Climate = c('id','Date','prec','tmax','tmin')) %>%
        dplyr::rename(id = 'id1') %>%
        dplyr::select(id, everything(.))
    } else {
      clim_data <- fst::read_fst(paste0(obsDir,'/',tolower(county),'.fst'))
      clim_data <- clim_data %>%
        tidyr::nest(Climate = c('id','Date','prec','tmax','tmin','srad')) %>%
        dplyr::rename(id = 'id1') %>%
        dplyr::select(id, everything(.))
    }
  } else {
    if(!file.exists(paste0(futDir,'/',tolower(county),'.fst'))){
      clim_data <- fst::read_fst(paste0(futDir,'/',tolower(county),'_prec_temp.fst'))
      clim_data <- clim_data %>%
        tidyr::nest(Climate = c('id','Date','prec','tmax','tmin')) %>%
        dplyr::rename(id = 'id1') %>%
        dplyr::select(id, everything(.))
    } else {
      clim_data <- fst::read_fst(paste0(futDir,'/',tolower(county),'.fst'))
      clim_data <- clim_data %>%
        tidyr::nest(Climate = c('id','Date','prec','tmax','tmin','srad')) %>%
        dplyr::rename(id = 'id1') %>%
        dplyr::select(id, everything(.))
    }
  }
  
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
  
  clim_data <- clim_data %>% 
    dplyr::mutate(county = county, semester = seasons) %>%
    dplyr::mutate(id = as.integer(id))
  clim_data <- dplyr::inner_join(clim_data, Soil) %>%
    dplyr::mutate(soilcp = ifelse(is.na(soilcp), 100, soilcp))
  
  plan(multiprocess)
  index_by_pixel <- clim_data %>%
    dplyr::group_split(id) %>%
    furrr::future_map(.f = reading_run_pys) %>% 
    dplyr::bind_rows()
  gc()
  gc(reset = T)
  
  data_with_C_index <- index_by_pixel %>% dplyr::select(-climate)
  if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
  out <- paste0(outDir,'/',county,'_',period,'_prec_temp.fst')
  if(!file.exists(out)){
    fst::write_fst(x = data_with_C_index %>% dplyr::select(-season, -sum.ind), path = out)
  }
  
}
calc_indices(country = 'Malawi',
             county  = 'Central Region',
             iso3c   = 'MWI',
             adm_lvl = 1,
             seasons = 1,
             gcm     = NULL,
             period  = '1985_2015',
             time    = 'past')
gcmList <- c("ipsl_cm5a_mr","miroc_esm_chem","ncc_noresm1_m")
periodList <- c('2021_2045','2041_2065')
for(gcm in gcmList){
  for(period in periodList){
    calc_indices(country = 'Malawi',
                 county  = 'Central Region',
                 iso3c   = 'MWI',
                 adm_lvl = 1,
                 seasons = 1,
                 gcm     = gcm,
                 period  = period,
                 time    = 'future')
  }
}
