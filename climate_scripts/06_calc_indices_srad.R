# =-------------------------------------------------------------------
# Climate Indexes Project
# Alejandra E. - Harold A. 
# April - 2020
# =-------------------------------------------------------------------

# =--------------------
rm(list = ls())
gc(reset = TRUE)# 

options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, ncdf4, sf, future, furrr, lubridate, glue, cowsay, vroom, sp, fst, compiler, ggspatial))

# =---------------
root <- '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles'

# =---------- Functions... 

# =-----------------------------------------------
# Estas estan probadas que funcionan bien...
# =-----------------------------------------------

source(paste0(root,'/scripts/indices.R'))

# Main functions
run_each_semester <<- function(one, semester){
  
  one2 <- one %>%
    dplyr::select(id, year) %>%
    unique()
  
  one2 <- one2 %>%
    dplyr::mutate(season = purrr::map(.x =year, .f = function(y){
      oi = list(rsum.lapply(x = filter(one, year == y)$prec %>% as.numeric, n = 150))[[1]] } ) )
  
  one2 <- one2 %>%
    mutate(sum.ind = purrr::map(.x = season, function(y){list(y[[which.max(cumulative.r.sum(y[[1]]))]])[[1]]}))

  return(one2)}

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

# =-----
GSEASProcess <- function(df){
  
  cons <- unique(df$id)
  path <- unique(df$path)
  
  out_all <- df
  srad_miss <- which(is.na(out_all$srad)) %>%
    purrr::map(function(i){
      out_all$srad[(i-1):(i+1)] %>% mean(., na.rm = T)
    }) %>%
    unlist()
  soilcp <- out_all$soilcp %>% unique
  out_all$srad[which(is.na(out_all$srad))] <- srad_miss
  if(!is.na(soilcp)){
    
    watbal_loc <- watbal_wrapper(out_all=out_all, soilcp=soilcp) # If we need more indexes are here
    watbal_loc$TAV <- (watbal_loc$tmin + watbal_loc$tmax)/2
    
    # Esta linea linea se debe modificar. 
    fst::write_fst(watbal_loc, glue::glue('{path}/wb_{cons}.fst'))
    
    watbal_loc <- watbal_loc[,c('Date','id','x','y','TAV','ERATIO')]
    watbal_loc$GDAY <- ifelse(watbal_loc$TAV >= 6 & watbal_loc$ERATIO >= 0.35, yes=1, no=0)
    
    watbal_loc$Year <- lubridate::year(watbal_loc$Date %>% lubridate::as_date())
    years_analysis <- watbal_loc$Year %>% unique()
    
    ### CONDITIONS TO HAVE IN ACCOUNT
    # Length of growing season per year
    # Start: 5-consecutive growing days.
    # End: 12-consecutive non-growing days.
    
    # Run process by year
    lgp_year_pixel <- lapply(1:length(years_analysis), function(k){
      
      # Subsetting by year
      watbal_year <- watbal_loc[watbal_loc$Year==years_analysis[k],]
      
      # Calculate sequences of growing and non-growing days within year
      runsDF <- rle(watbal_year$GDAY)
      runsDF <- data.frame(Lengths=runsDF$lengths, Condition=runsDF$values)
      runsDF$Condition <- runsDF$Condition %>% tidyr::replace_na(replace = 0)
      
      # Identify start and extension of each growing season during year
      if(!sum(runsDF$Lengths[runsDF$Condition==1] < 5) == length(runsDF$Lengths[runsDF$Condition==1])){
        
        LGP <- 0; LGP_seq <- 0
        for(i in 1:nrow(runsDF)){
          if(runsDF$Lengths[i] >= 5 & runsDF$Condition[i] == 1){
            LGP <- LGP + 1
            LGP_seq <- c(LGP_seq, LGP)
            LGP <- 0
          } else {
            if(LGP_seq[length(LGP_seq)]==1){
              if(runsDF$Lengths[i] >= 12 & runsDF$Condition[i] == 0){
                LGP <- 0
                LGP_seq <- c(LGP_seq, LGP)
              } else {
                LGP <- LGP + 1
                LGP_seq <- c(LGP_seq, LGP)
                LGP <- 0
              }
            } else {
              LGP <- 0
              LGP_seq <- c(LGP_seq, LGP)
            }
          }
        }
        LGP_seq <- c(LGP_seq, LGP)
        LGP_seq <- LGP_seq[-c(1, length(LGP_seq))]
        runsDF$gSeason <- LGP_seq; rm(i, LGP, LGP_seq)
        LGP_seq <- as.list(split(which(runsDF$gSeason==1), cumsum(c(TRUE, diff(which(runsDF$gSeason==1))!=1))))
        
        # Calculate start date and extension of each growing season by year and pixel
        growingSeason <- lapply(1:length(LGP_seq), function(g){
          
          LGP_ini <- sum(runsDF$Lengths[1:(min(LGP_seq[[g]])-1)]) + 1
          LGP <- sum(runsDF$Lengths[LGP_seq[[g]]])
          results <- data.frame(id=watbal_loc$id %>% unique, year=years_analysis[k], gSeason=g, SLGP=LGP_ini, LGP=LGP)
          return(results)
          
        })
        growingSeason <- do.call(rbind, growingSeason)
        if(nrow(growingSeason)>2){
          growingSeason <- growingSeason[rank(-growingSeason$LGP) %in% 1:2,]
          growingSeason$gSeason <- rank(growingSeason$SLGP)
          growingSeason <- growingSeason[order(growingSeason$gSeason),]
        }
        
      } else {
        
        growingSeason <- data.frame(id=watbal_loc$id %>% unique, year=years_analysis[k], gSeason = 1:2, SLGP = NA, LGP = NA)
        
      }
      
      print(k)
      return(growingSeason)
      
    })
    lgp_year_pixel <- do.call(rbind, lgp_year_pixel); rownames(lgp_year_pixel) <- 1:nrow(lgp_year_pixel)
    
    return(lgp_year_pixel)
    
  } else {
    return(cat('Check your inputs\n'))
  }
  
}

# =-----
# Aqui se calcula el indicador de ndws
ndws_index <- function(row_1){
  Eratio <- row_1 %>% dplyr::select(-season, -sum.ind) %>% unnest() %>% 
    .[row_1$season[[1]][[2]][[2]],] %>% pull(ERATIO)
  
  ndws <- calc_wsdays(Eratio, season_ini=1, season_end=150, e_thresh=0.5)
  
  row_1 <- mutate(row_1, ndws = ndws) %>% dplyr::select(-`row_number()`, -season, -sum.ind)
  
  return(row_1)}



# =----------------------------------------------------------------
# =----------------------------------------------------------------
cowsay::say('Reading data ---- :v', by = 'ghost')
# =----------------------------------------------------------------

# Please check the correct county name within: Country_Counts.xlsx
country <- 'Ethiopia'
county  <- 'Arsi'
iso3c   <- 'ETH'
adm_lvl <- 2


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

gcmList <- c('ipsl_cm5a_mr','miroc_esm_chem','ncc_noresm1_m')
timList <- c('past','future')
periodList <- c('2021_2045','2041_2065')



# 
calc_indices <- function(country = 'Ethiopia', county = 'Arsi', seasons = 1, gcm = 'ipsl_cm5a_mr', period = '2041_2065', time = 'future'){
  
  # country = 'Ethiopia'; county = 'Arsi'; seasons = 1; period = '1985-2015'; time = 'past'
  
  # Paths
  obsDir <- paste0(root,'/data/observational_data/',tolower(country))
  futDir <- paste0(root,'/data/bc_quantile_0_05deg_lat_county/',tolower(country),'/',gcm,'/',period)
  outDir <- ifelse(test = time == 'past',
                   yes  = paste0(root,'/results/',country,'/',time),
                   no   = paste0(root,'/results/',country,'/',time,'/',gcm,'/',period))
  Med_Dir <- ifelse(test = time == 'past',
                   yes  = paste0(root,'/data/metadata/watbal/',country,'/',time),
                   no   = paste0(root,'/data/metadata/watbal/',country,'/',time,'/',gcm,'/',period))
  # Soil data
  Soil <- fst::read.fst(paste0(root,'/data/soilcp_data.fst')) %>%
    tibble::as_tibble() %>%
    dplyr::select(id, soilcp) %>%
    dplyr::filter(id %in% dplyr::pull(crd, id))
  
  
  # Aqui debo de revisar que debo de leer --- al ser radiación deberia ser los archivos normales. 
  # Deberia poner que si no tiene la radiación mande un aviso o pare. 
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
  
  # Cambiar esto por el multicluster para ver como le va... 
  # Abrir el cluster :: Mod esa linea en caso de que se deseen menos cluster...
  # cores <- parallel::detectCores() - 1
  cores <- 15
  plan(cluster, workers = cores)
  
  # parte inicial de los indices climaticos.
  index_by_pixel <- clim_data %>% 
    dplyr::group_split(id) %>% 
    furrr::future_map(.f = reading_run_pys) %>% 
    dplyr::bind_rows()
  
  
  if(!dir.exists(Med_Dir)){dir.create(Med_Dir, recursive = T)}
  
  # 3 indices basados en radiacion. 
  first_Ind <- clim_data %>%
    mutate(path = Med_Dir) %>% 
    group_split(row_number()) %>% 
    purrr::map(.f = function(x){x %>% unnest() %>% dplyr::select(-`row_number()`)}) %>% 
    furrr::future_map(.f = GSEASProcess)
  
  # 
  watbal_data <- clim_data %>% 
    dplyr::select(id) %>%  
    mutate(path_w = glue::glue('{Med_Dir}/wb_{id}.fst')) %>% 
    mutate(watbal = furrr::future_map(.x = path_w, .f = function(x){fst::fst(x) %>% as_tibble(.) %>%  dplyr::mutate(year = lubridate::year(Date %>% lubridate::as_date())) %>% dplyr::select(-id1)})) %>% 
    dplyr::select(-path_w) %>% 
    unnest() %>% 
    dplyr::select(-id1) %>%
    nest(-id, -x, -y, -ISO3, -Country, -semester,-year)
  
  # =-
  tictoc::tic()
  full_data <- full_join(watbal_data ,  index_by_pixel) %>% 
    group_split(row_number()) %>%
    furrr::future_map(.f = ndws_index) %>%
    bind_rows()
  tictoc::toc() # 3.58 seg con 2 cores cabe aclarar las pruebas son para 10 pixels. 
  
  
  # Indices de radiacion solar. 
  rad_index <- first_Ind %>% 
    bind_rows() %>% as_tibble() %>% 
    full_join(full_data)
  
  
  climate_index <- fst::fst(paste0(outDir,'/',county,'_',period,'_prec_temp.fst') ) %>% as_tibble()
  
  all_index <- rad_index %>% 
    dplyr::select(-data) %>% 
    full_join(climate_index , . ) %>% # 
    mutate(county = county) %>% 
    dplyr::select(id, x, y, ISO3, Country, county, everything(.))
  
  
  if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
  out <- paste0(outDir,'/',county,'_',period,'.fst')
  if(!file.exists(out)){fst::write_fst(x = all_index , path = out)}
  
return(all_index)}

calc_indices(country = 'Ethiopia',
             county  = 'Arsi',
             seasons = 1,
             gcm     = NULL,
             period  = '1985_2015',
             time    = 'past')


for(gcm in gcmList){
  for(period in periodList){
    calc_indices(country = 'Ethiopia',
                 county  = 'Arsi',
                 seasons = 1,
                 gcm     = gcm,
                 period  = period,
                 time    = 'future')
  }
}

