# Calculate growing season indices
# A. Esquivel and H. Achicanoy
# CIAT, 2020

options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, ncdf4, sf, future, furrr, lubridate, glue, vroom, sp, fst, compiler))

OSys <- Sys.info()[1]
root <<- switch(OSys,
                'Linux'   = '/home/jovyan/work/cglabs',
                'Windows' = '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles')

calc_indices <- function(country = 'Jordan',
                         county  = 'jordan_rcp45',
                         iso3c   = 'JOR',
                         adm_lvl = 0,
                         seasons = list(s1 = c(10:12,1:4)), # Seasons manually defined
                         n_ssns  = NULL,    # 2-seasons automatically defined
                         n_wtts  = NULL,  # 100-wettest days
                         gcm     = NULL,
                         period  = '1980_2050',
                         time    = 'past',
                         big_cnt = TRUE,
                         ncores  = 20){
  
  country <<- country
  county  <<- county
  
  # # Load county shapefile
  # shp <<- raster::shapefile(paste0(root,'/data/shps/',country,'/',iso3c,'_adm',adm_lvl,'.shp'))
  # glue::glue('shp <<- shp[shp@data$NAME_{adm_lvl} == county,]') %>%
  #   as.character %>%
  #   parse(text = .) %>%
  #   eval(expr = ., envir = .GlobalEnv)
  # 
  # # Load id coords
  # crd <- vroom(paste0(root,'/data/id_all_country.csv'), delim = ',')
  # crd <- crd %>%
  #   dplyr::filter(Country == country)
  # pnt <- crd %>% dplyr::select(x,y) %>% sp::SpatialPoints(coords = .)
  # crs(pnt) <- crs(shp)
  # # Filter coordinates that are present in the county
  # pnt <- sp::over(pnt, shp) %>% data.frame %>% dplyr::select(ISO) %>% complete.cases() %>% which()
  # crd <- crd[pnt,]
  # crd <<- crd
  
  timList <- c('past','future')
  periodList <- c('2021_2045','2041_2065')
  
  source(paste0(root,'/scripts/indices.R'))
  
  # Paths
  obsDir <- paste0(root,'/data/observational_data/',tolower(country))
  futDir <- paste0(root,'/data/bc_quantile_0_05deg_lat_county/',tolower(country),'/',gcm,'/',period)
  outDir <- ifelse(test = time == 'past',
                   yes  = paste0(root,'/results/',country,'/',time),
                   no   = paste0(root,'/results/',country,'/',time,'/',gcm,'/',period))
  
  if(big_cnt){
    out <- paste0(outDir,'/',county,'_',period,'_gSeasons.fst')
  } else {
    out <- paste0(outDir,'/',county,'_',period,'_gSeasons_idw.fst')
  }
  
  if(!file.exists(out)){
    # Soil data
    Soil <- fst::read.fst(paste0(root,'/data/soilcp_data_jordan.fst')) %>%
      tibble::as_tibble() %>%
      dplyr::select(id, soilcp)
    
    # Read climate data
    if(time == 'past'){
      if(!file.exists(paste0(obsDir,'/',tolower(county),'.fst'))){
        clim_data <- fst::read_fst(paste0(obsDir,'/',tolower(county),'_prec_temp.fst'))
        clim_data <- clim_data %>%
          tidyr::nest(Climate = c('id','Date','prec','tmax','tmin')) %>%
          dplyr::rename(id = 'id1') %>%
          dplyr::select(id, everything(.))
      } else {
        clim_data <- fst::read_fst(paste0(obsDir,'/',tolower(county),'.fst'))
        clim_data$id1 <- clim_data$id
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
    
    if(big_cnt){
      set.seed(1235)
      sample_n  <- nrow(clim_data)*0.3
      id_sample <- sample(unique(clim_data$id), sample_n) 
      clim_data <- dplyr::filter(clim_data, id %in% id_sample)
    }
    
    run_pixel <- function(id = 362540){
      
      cat(' --- Obtain complete time series per pixel\n')
      tbl <- clim_data$Climate[[which(clim_data$id == id)]]
      tbl <- tbl %>%
        dplyr::mutate(year  = lubridate::year(as.Date(Date)),
                      month = lubridate::month(as.Date(Date)))
      years <- tbl$year %>% unique
      
      cat(' --- Calculate water balance for complete time series\n')
      soilcp <- Soil$soilcp[Soil$id == id]
      watbal_loc <- watbal_wrapper(out_all = tbl, soilcp = soilcp)
      # watbal_loc$IRR <- watbal_loc$Etmax - watbal_loc$prec
      
      tbl <- tbl %>%
        dplyr::mutate(ERATIO = watbal_loc$ERATIO,
                      TAV    = (watbal_loc$tmin + watbal_loc$tmax)/2,
                      GDAY   = ifelse(TAV >= 6 & ERATIO >= 0.35, yes=1, no=0))
      
      cat(' --- Estimate growing seasons from water balance\n')
      
      ### CONDITIONS TO HAVE IN ACCOUNT
      # Length of growing season per year
      # Start: 5-consecutive growing days.
      # End: 12-consecutive non-growing days.
      
      # Run process by year
      lgp_year_pixel <- lapply(1:length(years), function(k){
        
        # Subsetting by year
        watbal_year <- tbl[tbl$year==years[k],]
        
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
            results <- data.frame(id=tbl$id %>% unique, year=years[k], gSeason=g, SLGP=LGP_ini, LGP=LGP)
            return(results)
            
          })
          growingSeason <- do.call(rbind, growingSeason)
          if(nrow(growingSeason)>2){
            growingSeason <- growingSeason[rank(-growingSeason$LGP) %in% 1:2,]
            growingSeason$gSeason <- rank(growingSeason$SLGP)
            growingSeason <- growingSeason[order(growingSeason$gSeason),]
          }
          
        } else {
          
          growingSeason <- data.frame(id=tbl$id %>% unique, year=years[k], gSeason = 1:2, SLGP = NA, LGP = NA)
          
        }
        
        print(k)
        return(growingSeason)
        
      })
      lgp_year_pixel <- do.call(rbind, lgp_year_pixel); rownames(lgp_year_pixel) <- 1:nrow(lgp_year_pixel)
      all <- lgp_year_pixel
      
      return(all)
      
    }
    
    plan(cluster, workers = ncores)
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
    index_by_pixel <- index_by_pixel %>% dplyr::select(id, x, y, ISO3, Country, county, year, dplyr::everything(.))
    
    if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
    fst::write_fst(x = index_by_pixel, path = out)
    cat('>>> File created successfully ...\n')
  } else {
    cat('>>> File exists it is not necessary to create it again\n')
  }
  
}

# calc_indices(country = 'Jordan',
#              county  = 'jordan_rcp45',
#              iso3c   = 'JOR',
#              adm_lvl = 0,
#              seasons = list(s1 = c(10:12,1:4)), # Seasons manually defined
#              n_ssns  = NULL,    # 2-seasons automatically defined
#              n_wtts  = NULL,  # 100-wettest days
#              gcm     = NULL,
#              period  = '1980_2050',
#              time    = 'past',
#              big_cnt = TRUE,
#              ncores  = 5)

calc_indices(country = 'Jordan',
             county  = 'jordan_rcp85',
             iso3c   = 'JOR',
             adm_lvl = 0,
             seasons = list(s1 = c(10:12,1:4)), # Seasons manually defined
             n_ssns  = NULL,    # 2-seasons automatically defined
             n_wtts  = NULL,  # 100-wettest days
             gcm     = NULL,
             period  = '1980_2050',
             time    = 'past',
             big_cnt = FALSE,
             ncores  = 20)
