### Perform quantile-mapping bias correction of daily climate data
### H. Achicanoy - C. Navarro
### CIAT, 2020

# R options
options(warn = -1, scipen = 999)

# Load libraries
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(qmap, ncdf4, raster, tidyverse, compiler, vroom, gtools, fst))

# Parallelization
clusterExport <- local({
  gets <- function(n, v) { assign(n, v, envir = .GlobalEnv); NULL }
  function(cl, list, envir = .GlobalEnv) {
    ## do this with only one clusterCall--loop on slaves?
    for (name in list) {
      clusterCall(cl, gets, name, get(name, envir = envir))
    }
  }
})
createCluster <- function(noCores, logfile = "/dev/null", export = NULL, lib = NULL) {
  require(doSNOW)
  cl <- makeCluster(noCores, type = "SOCK", outfile = logfile)
  if(!is.null(export)) clusterExport(cl, export)
  if(!is.null(lib)) {
    plyr::l_ply(lib, function(dum) { 
      clusterExport(cl, "dum", envir = environment())
      clusterEvalQ(cl, library(dum, character.only = TRUE))
    })
  }
  registerDoSNOW(cl)
  return(cl)
}

# Quantile-mapping bias correction function for available pixels with solar radiation from NASA within a country
BC_Qmap <- function(country   = "Pakistan",
                    rcp       = "rcp85",
                    gcm       = "ipsl_cm5a_mr",
                    period    = "2021_2045",
                    srad_avlb = T)
{
  
  cat(paste0(' *** Performing Quantile-mapping bias correction for available pixels with SRAD (NASA) within ',country,' in the period',period,', using: ',rcp,', GCM: ',gcm,'***\n'))
  
  # Establish directories
  root       <<- "//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data"
  obsDir     <<- paste0(root,"/Chirps_Chirts")
  gcmHistDir <<- paste0(root,"/gcm_0_05deg_lat/",tolower(country),"/",gcm,"/1971_2000")
  gcmFutDir  <<- paste0(root,"/gcm_0_05deg_lat/",tolower(country),"/",gcm,"/",period,"/",rcp)
  outDir     <<- ifelse(srad_avlb,paste0(root,"/bc_quantile_0_05deg_lat"),paste0(root,"/bc_quantile_0_05deg_lat_no_srad"))
  if(!dir.exists(outDir)){dir.create(outDir)}
  
  # Identify pixels in country
  px_id <- vroom::vroom("//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data/id_country.csv", delim=',')
  px_id <- px_id %>% dplyr::filter(Country == country)
  
  if(srad_avlb){
    # Identify available pixels with solar radiation from NASA
    avlb_px <- list.files(path=paste0(root,'/NASA'),full.names=F) %>% gsub('.fst','',.)
    avlb_px <- avlb_px[avlb_px %in% px_id$id]
    avlb_px <- gtools::mixedsort(avlb_px)
  } else {
    avlb_px <- px_id$id
    avlb_px <- gtools::mixedsort(avlb_px)
  }
  
  
  process_by_px <- avlb_px %>% purrr::map(.f = function(px){
    
    px <<- px
    if(!dir.exists(paste0(outDir,'/',gcm,'/',period,'/',rcp))){dir.create(paste0(outDir,'/',gcm,'/',period,'/',rcp))}
    out <- paste0(outDir,'/',gcm,'/',period,'/',rcp,'/',px,'.fst')
    if(!file.exists(out)){
      
      cat(paste0('>>> Processing pixel: ',px,'...\n'))
      
      cat(paste0('> Load historical observational data\n'))
      # Historical time series (observations)
      ts <- list.files(path=obsDir,full.names=F,pattern='*.fst$') %>% sort()
      cl <- createCluster(30, export = list("ts","obsDir","px"), lib = list("tidyverse","fst"))
      # Temperatures, precipitation, and solar radiation data
      temp_prec <- ts %>% parallel::parLapply(cl, ., function(i){
        # temp_prec <- ts %>% purrr::map(.f = function(i){
        df <- fst::read_fst(paste0(obsDir,'/',i),from=as.numeric(px),to=as.numeric(px))
        return(df)
      }) %>%
        do.call(rbind, .)
      parallel::stopCluster(cl)
      temp_prec <- temp_prec %>% dplyr::arrange(Date)
      solr_radt <- fst::read_fst(paste0(root,"/NASA/",px,".fst"))
      
      obs_climate <- merge(x = temp_prec, y = solr_radt %>% dplyr::select('Date','srad'), by = 'Date'); rm(temp_prec, solr_radt, ts)
      px_crds <<- obs_climate[,c('x','y')] %>% unique
      
      cat(paste0('> Load historical GCMs data\n'))
      # Historical time series (GCMs)
      prec_fls <<- list.files(path=paste0(gcmHistDir,'/by-month'),pattern='^prec',full.names=T)
      cl       <- createCluster(30, export = list("prec_fls","px_crds","px"), lib = list("tidyverse","raster","fst"))
      prec <- prec_fls %>%
        parallel::parLapply(cl, ., function(x){
          date <- basename(x)
          date <- date %>% gsub('prec_','',.) %>% gsub('.nc','',.) %>% gsub('_','-',.)
          rsts <- raster::stack(x)
          prec <- raster::extract(x = rsts, y = px_crds) %>% as.numeric
          df   <- data.frame(Date = paste0(date,'-',c(paste0('0',1:9),10:length(prec))),
                             id   = px,
                             x    = px_crds$x,
                             y    = px_crds$y,
                             prec = prec)
          return(df)
        }) %>% do.call(rbind, .)
      parallel::stopCluster(cl)
      rm(prec_fls)
      prec <- prec %>% dplyr::arrange(Date)
      
      tmax_fls <<- list.files(path=paste0(gcmHistDir,'/by-month'),pattern='^tmax',full.names=T)
      cl       <- createCluster(30, export = list("tmax_fls","px_crds","px"), lib = list("tidyverse","raster","fst"))
      tmax <- tmax_fls %>%
        parallel::parLapply(cl, ., function(x){
          date <- basename(x)
          date <- date %>% gsub('tmax_','',.) %>% gsub('.nc','',.) %>% gsub('_','-',.)
          rsts <- raster::stack(x)
          tmax <- raster::extract(x = rsts, y = px_crds) %>% as.numeric
          df   <- data.frame(Date = paste0(date,'-',c(paste0('0',1:9),10:length(tmax))),
                             id   = px,
                             x    = px_crds$x,
                             y    = px_crds$y,
                             tmax = tmax)
          return(df)
        }) %>% do.call(rbind, .)
      parallel::stopCluster(cl)
      rm(tmax_fls)
      tmax <- tmax %>% dplyr::arrange(Date)
      
      tmin_fls <<- list.files(path=paste0(gcmHistDir,'/by-month'),pattern='^tmin',full.names=T)
      cl       <- createCluster(30, export = list("tmin_fls","px_crds","px"), lib = list("tidyverse","raster","fst"))
      tmin <- tmin_fls %>%
        parallel::parLapply(cl, ., function(x){
          date <- basename(x)
          date <- date %>% gsub('tmin_','',.) %>% gsub('.nc','',.) %>% gsub('_','-',.)
          rsts <- raster::stack(x)
          tmin <- raster::extract(x = rsts, y = px_crds) %>% as.numeric
          df   <- data.frame(Date = paste0(date,'-',c(paste0('0',1:9),10:length(tmin))),
                             id   = px,
                             x    = px_crds$x,
                             y    = px_crds$y,
                             tmin = tmin)
          return(df)
        }) %>% do.call(rbind, .)
      parallel::stopCluster(cl)
      rm(tmin_fls)
      tmin <- tmin %>% dplyr::arrange(Date)
      
      if(srad_avlb){
        srad_fls <<- list.files(path=paste0(gcmHistDir,'/by-month'),pattern='^rsds',full.names=T)
        cl       <- createCluster(30, export = list("srad_fls","px_crds","px"), lib = list("tidyverse","raster","fst"))
        srad <- srad_fls %>%
          parallel::parLapply(cl, ., function(x){
            date <- basename(x)
            date <- date %>% gsub('rsds_','',.) %>% gsub('.nc','',.) %>% gsub('_','-',.)
            rsts <- raster::stack(x)
            srad <- raster::extract(x = rsts, y = px_crds) %>% as.numeric
            df   <- data.frame(Date = paste0(date,'-',c(paste0('0',1:9),10:length(srad))),
                               id   = px,
                               x    = px_crds$x,
                               y    = px_crds$y,
                               srad = srad)
            return(df)
          }) %>% do.call(rbind, .)
        parallel::stopCluster(cl)
        rm(srad_fls)
        srad <- srad %>% dplyr::arrange(Date)
        
        his_climate <- cbind(prec,
                             tmax %>% dplyr::select(tmax),
                             tmin %>% dplyr::select(tmin),
                             srad %>% dplyr::select(srad))
        rm(prec, tmax, tmin, srad)
        his_climate$srad <- his_climate$srad * 0.0864 # W/m-2 to MJ/m-2/day-1
      } else {
        his_climate <- cbind(prec,
                             tmax %>% dplyr::select(tmax),
                             tmin %>% dplyr::select(tmin))
        rm(prec, tmax, tmin)
      }
      
      cat(paste0('> Load future GCMs data\n'))
      # Future time series (GCMs)
      prec_fls <<- list.files(path=paste0(gcmFutDir,'/by-month'),pattern='^prec',full.names=T)
      cl       <- createCluster(30, export = list("prec_fls","px_crds","px"), lib = list("tidyverse","raster","fst"))
      prec <- prec_fls %>%
        parallel::parLapply(cl, ., function(x){
          date <- basename(x)
          date <- date %>% gsub('prec_','',.) %>% gsub('.nc','',.) %>% gsub('_','-',.)
          rsts <- raster::stack(x)
          prec <- raster::extract(x = rsts, y = px_crds) %>% as.numeric
          df   <- data.frame(Date = paste0(date,'-',c(paste0('0',1:9),10:length(prec))),
                             id   = px,
                             x    = px_crds$x,
                             y    = px_crds$y,
                             prec = prec)
          return(df)
        }) %>% do.call(rbind, .)
      parallel::stopCluster(cl)
      rm(prec_fls)
      prec <- prec %>% dplyr::arrange(Date)
      
      tmax_fls <<- list.files(path=paste0(gcmFutDir,'/by-month'),pattern='^tmax',full.names=T)
      cl       <- createCluster(30, export = list("tmax_fls","px_crds","px"), lib = list("tidyverse","raster","fst"))
      tmax <- tmax_fls %>%
        parallel::parLapply(cl, ., function(x){
          date <- basename(x)
          date <- date %>% gsub('tmax_','',.) %>% gsub('.nc','',.) %>% gsub('_','-',.)
          rsts <- raster::stack(x)
          tmax <- raster::extract(x = rsts, y = px_crds) %>% as.numeric
          df   <- data.frame(Date = paste0(date,'-',c(paste0('0',1:9),10:length(tmax))),
                             id   = px,
                             x    = px_crds$x,
                             y    = px_crds$y,
                             tmax = tmax)
          return(df)
        }) %>% do.call(rbind, .)
      parallel::stopCluster(cl)
      rm(tmax_fls)
      tmax <- tmax %>% dplyr::arrange(Date)
      
      tmin_fls <<- list.files(path=paste0(gcmFutDir,'/by-month'),pattern='^tmin',full.names=T)
      cl       <- createCluster(30, export = list("tmin_fls","px_crds","px"), lib = list("tidyverse","raster","fst"))
      tmin <- tmin_fls %>%
        parallel::parLapply(cl, ., function(x){
          date <- basename(x)
          date <- date %>% gsub('tmin_','',.) %>% gsub('.nc','',.) %>% gsub('_','-',.)
          rsts <- raster::stack(x)
          tmin <- raster::extract(x = rsts, y = px_crds) %>% as.numeric
          df   <- data.frame(Date = paste0(date,'-',c(paste0('0',1:9),10:length(tmin))),
                             id   = px,
                             x    = px_crds$x,
                             y    = px_crds$y,
                             tmin = tmin)
          return(df)
        }) %>% do.call(rbind, .)
      parallel::stopCluster(cl)
      rm(tmin_fls)
      tmin <- tmin %>% dplyr::arrange(Date)
      
      if(srad_avlb){
        srad_fls <<- list.files(path=paste0(gcmFutDir,'/by-month'),pattern='^rsds',full.names=T)
        cl       <- createCluster(30, export = list("srad_fls","px_crds","px"), lib = list("tidyverse","raster","fst"))
        srad <- srad_fls %>%
          parallel::parLapply(cl, ., function(x){
            date <- basename(x)
            date <- date %>% gsub('rsds_','',.) %>% gsub('.nc','',.) %>% gsub('_','-',.)
            rsts <- raster::stack(x)
            srad <- raster::extract(x = rsts, y = px_crds) %>% as.numeric
            df   <- data.frame(Date = paste0(date,'-',c(paste0('0',1:9),10:length(srad))),
                               id   = px,
                               x    = px_crds$x,
                               y    = px_crds$y,
                               srad = srad)
            return(df)
          }) %>% do.call(rbind, .)
        parallel::stopCluster(cl)
        rm(srad_fls)
        srad <- srad %>% dplyr::arrange(Date)
        
        fut_climate <- cbind(prec,
                             tmax %>% dplyr::select(tmax),
                             tmin %>% dplyr::select(tmin),
                             srad %>% dplyr::select(srad))
        rm(prec, tmax, tmin, srad)
        fut_climate$srad <- fut_climate$srad * 0.0864 # W/m-2 to MJ/m-2/day-1
      } else {
        fut_climate <- cbind(prec,
                             tmax %>% dplyr::select(tmax),
                             tmin %>% dplyr::select(tmin))
        rm(prec, tmax, tmin)
      }
      
      cat(paste0('> Applying Quantile mapping ...\n'))
      
      prec_fit <- qmap::fitQmap(obs=obs_climate$prec, mod=his_climate$prec, method="RQUANT", qstep=0.01, wet.day=TRUE, na.rm=TRUE)
      tmax_fit <- qmap::fitQmap(obs=obs_climate$tmax, mod=his_climate$tmax, method="RQUANT", qstep=0.01, wet.day=FALSE, na.rm=TRUE)
      tmin_fit <- qmap::fitQmap(obs=obs_climate$tmin, mod=his_climate$tmin, method="RQUANT", qstep=0.01, wet.day=FALSE, na.rm=TRUE)
      if(srad_avlb){
        srad_fit <- qmap::fitQmap(obs=obs_climate$srad, mod=his_climate$srad, method="RQUANT", qstep=0.01, wet.day=FALSE, na.rm=TRUE)
        qmap_results <- tibble::tibble(id   = px,
                                       x    = px_crds$x,
                                       y    = px_crds$y,
                                       prec = list(prec_fit),
                                       tmax = list(tmax_fit),
                                       tmin = list(tmin_fit),
                                       srad = list(srad_fit))
        if(!dir.exists(paste0(outDir,'/',gcm,'/',period,'/',rcp,'/qmap_fit'))){dir.create(paste0(outDir,'/',gcm,'/',period,'/',rcp,'/qmap_fit'))}
        fst::write_fst(qmap_results,paste0(outDir,'/',gcm,'/',period,'/',rcp,'/qmap_fit/',px,'.fst'))
        
        fut_clmt_bc <- data.frame(Date = fut_climate$Date,
                                  id   = fut_climate$id,
                                  x    = fut_climate$x,
                                  y    = fut_climate$y,
                                  prec = qmap::doQmap(x=fut_climate$prec, prec_fit, type="linear"),
                                  tmax = qmap::doQmap(x=fut_climate$tmax, tmax_fit, type="linear"),
                                  tmin = qmap::doQmap(x=fut_climate$tmin, tmin_fit, type="linear"),
                                  srad = qmap::doQmap(x=fut_climate$srad, srad_fit, type="linear"))
        
        his_clmt_bc <- data.frame(Date = his_climate$Date,
                                  id   = his_climate$id,
                                  x    = his_climate$x,
                                  y    = his_climate$y,
                                  prec = qmap::doQmap(x=his_climate$prec, prec_fit, type="linear"),
                                  tmax = qmap::doQmap(x=his_climate$tmax, tmax_fit, type="linear"),
                                  tmin = qmap::doQmap(x=his_climate$tmin, tmin_fit, type="linear"),
                                  srad = qmap::doQmap(x=his_climate$srad, srad_fit, type="linear"))
      } else {
        qmap_results <- tibble::tibble(id   = px,
                                       x    = px_crds$x,
                                       y    = px_crds$y,
                                       prec = list(prec_fit),
                                       tmax = list(tmax_fit),
                                       tmin = list(tmin_fit))
        
        fut_clmt_bc <- data.frame(Date = fut_climate$Date,
                                  id   = fut_climate$id,
                                  x    = fut_climate$x,
                                  y    = fut_climate$y,
                                  prec = qmap::doQmap(x=fut_climate$prec, prec_fit, type="linear"),
                                  tmax = qmap::doQmap(x=fut_climate$tmax, tmax_fit, type="linear"),
                                  tmin = qmap::doQmap(x=fut_climate$tmin, tmin_fit, type="linear"))
        
        his_clmt_bc <- data.frame(Date = his_climate$Date,
                                  id   = his_climate$id,
                                  x    = his_climate$x,
                                  y    = his_climate$y,
                                  prec = qmap::doQmap(x=his_climate$prec, prec_fit, type="linear"),
                                  tmax = qmap::doQmap(x=his_climate$tmax, tmax_fit, type="linear"),
                                  tmin = qmap::doQmap(x=his_climate$tmin, tmin_fit, type="linear"))
      }
      fst::write_fst(fut_clmt_bc,out)
      cat(paste0('Pixel: ',px,' future was correctly processed.\n'))
      if(!dir.exists(paste0(outDir,'/',gcm,'/1971_2000/',rcp))){dir.create(paste0(outDir,'/',gcm,'/1971_2000/',rcp))}
      fst::write_fst(his_clmt_bc,paste0(outDir,'/',gcm,'/1971_2000/',rcp,'/',px,'.fst'))
      cat(paste0('Pixel: ',px,' historical was correctly processed.\n'))
      
    } else {
      cat(paste0('Pixel: ',px,' is already processed.\n'))
    }
    
  })
    
}
BC_Qmap(country="Ethiopia",rcp="rcp85",gcm="ipsl_cm5a_mr",period="2021_2045",srad_avlb=T)

periodList <- '2021_2045' # c('2041_2065')
rcpList    <- 'rcp85'
gcmList    <- c("ipsl_cm5a_mr","miroc_esm_chem","ncc_noresm1_m")
ctryList   <- c('Ethiopia','Mali','Pakistan')
