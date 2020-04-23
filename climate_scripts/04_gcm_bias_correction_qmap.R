### Perform quantile-mapping bias correction of daily climate data
### H. Achicanoy - C. Navarro
### CIAT, 2020

# R options
options(warn = -1, scipen = 999)

# Load libraries
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(qmap, ncdf4, raster, tidyverse, compiler, vroom, gtools, fst))

# Paths
OSys <- Sys.info()[1]
root <<- switch(OSys,
                'Linux'   = '/home/jovyan/work/cglabs',
                'Windows' = '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles')

# Scripts
source(paste0(root,'/scripts/win_parallelization.R'))

# Quantile-mapping bias correction function for available pixels with solar radiation from NASA within a country
BC_Qmap <- function(country   = "Ethiopia",
                    county    = "Arsi",
                    rcp       = "rcp85",
                    gcm       = "ipsl_cm5a_mr",
                    period    = "2021_2045")
{
  
  bc_qmap <<- function(df_obs, df_his_gcm, df_fut_gcm){
    if('srad' %in% colnames(df_obs)){
      cat('> Fitting the Qmap function per variable\n')
      prec_fit <- qmap::fitQmap(obs=df_obs$prec, mod=df_his_gcm$prec, method="RQUANT", qstep=0.01, wet.day=TRUE, na.rm=TRUE)
      tmax_fit <- qmap::fitQmap(obs=df_obs$tmax, mod=df_his_gcm$tmax, method="RQUANT", qstep=0.01, wet.day=FALSE, na.rm=TRUE)
      tmin_fit <- qmap::fitQmap(obs=df_obs$tmin, mod=df_his_gcm$tmin, method="RQUANT", qstep=0.01, wet.day=FALSE, na.rm=TRUE)
      srad_fit <- qmap::fitQmap(obs=df_obs$srad, mod=df_his_gcm$srad, method="RQUANT", qstep=0.01, wet.day=FALSE, na.rm=TRUE)
      
      cat('> Doing bias correction per variable historical GCMs\n')
      bc_his_gcm <- df_his_gcm
      bc_his_gcm$prec <- qmap::doQmap(x=df_his_gcm$prec, prec_fit, type="linear")
      bc_his_gcm$tmax <- qmap::doQmap(x=df_his_gcm$tmax, tmax_fit, type="linear")
      bc_his_gcm$tmin <- qmap::doQmap(x=df_his_gcm$tmin, tmin_fit, type="linear")
      bc_his_gcm$srad <- qmap::doQmap(x=df_his_gcm$srad, srad_fit, type="linear")
      
      cat('> Doing bias correction per variable future GCMs\n')
      bc_fut_gcm <- df_fut_gcm
      bc_fut_gcm$prec <- qmap::doQmap(x=df_fut_gcm$prec, prec_fit, type="linear")
      bc_fut_gcm$tmax <- qmap::doQmap(x=df_fut_gcm$tmax, tmax_fit, type="linear")
      bc_fut_gcm$tmin <- qmap::doQmap(x=df_fut_gcm$tmin, tmin_fit, type="linear")
      bc_fut_gcm$srad <- qmap::doQmap(x=df_fut_gcm$srad, srad_fit, type="linear")
      
    } else {
      cat('> Fitting the Qmap function per variable\n')
      prec_fit <- qmap::fitQmap(obs=df_obs$prec, mod=df_his_gcm$prec, method="RQUANT", qstep=0.01, wet.day=TRUE, na.rm=TRUE)
      tmax_fit <- qmap::fitQmap(obs=df_obs$tmax, mod=df_his_gcm$tmax, method="RQUANT", qstep=0.01, wet.day=FALSE, na.rm=TRUE)
      tmin_fit <- qmap::fitQmap(obs=df_obs$tmin, mod=df_his_gcm$tmin, method="RQUANT", qstep=0.01, wet.day=FALSE, na.rm=TRUE)
      
      cat('> Doing bias correction per variable historical GCMs\n')
      bc_his_gcm <- df_his_gcm
      bc_his_gcm$prec <- qmap::doQmap(x=df_his_gcm$prec, prec_fit, type="linear")
      bc_his_gcm$tmax <- qmap::doQmap(x=df_his_gcm$tmax, tmax_fit, type="linear")
      bc_his_gcm$tmin <- qmap::doQmap(x=df_his_gcm$tmin, tmin_fit, type="linear")
      
      cat('> Doing bias correction per variable future GCMs\n')
      bc_fut_gcm <- df_fut_gcm
      bc_fut_gcm$prec <- qmap::doQmap(x=df_fut_gcm$prec, prec_fit, type="linear")
      bc_fut_gcm$tmax <- qmap::doQmap(x=df_fut_gcm$tmax, tmax_fit, type="linear")
      bc_fut_gcm$tmin <- qmap::doQmap(x=df_fut_gcm$tmin, tmin_fit, type="linear")
    }
    bc_data <- list(His = bc_his_gcm,
                    Fut = bc_fut_gcm)
    return(list(bc_data))
  }
  
  cat(paste0(' *** Performing Quantile-mapping bias correction for available pixels with SRAD (NASA) within ',county,' in the period ',period,', using: ',rcp,', GCM: ',gcm,'***\n'))
  cat(paste0('>>> Loading obs data\n'))
  
  obsDir <<- paste0(root,"/data/observational_data/",tolower(country))
  if(!file.exists(paste0(obsDir,'/',tolower(county),'.RDS'))){
    his_obs <<- readRDS(paste0(obsDir,'/',tolower(county),'_prec_temp.RDS'))
  } else {
    his_obs <<- readRDS(paste0(obsDir,'/',tolower(county),'.RDS'))
  }
  
  cat(paste0('>>> Loading historical GCM data\n'))
  hisGCMDir <<- paste0(root,"/data/gcm_0_05deg_lat_county/",tolower(country),"/",gcm,"/1971_2000")
  his_gcm <<- readRDS(paste0(hisGCMDir,'/',tolower(county),'.RDS'))
  
  cat(paste0('>>> Loading future GCM data\n'))
  futGCMDir <<- paste0(root,"/data/gcm_0_05deg_lat_county/",tolower(country),"/",gcm,"/",period)
  fut_gcm <<- readRDS(paste0(futGCMDir,'/',tolower(county),'.RDS'))
  
  his_gcm_bc <<- his_gcm
  fut_gcm_bc <<- fut_gcm
  
  cl <- createCluster(30, export = list("root","obsDir","his_obs","hisGCMDir","his_gcm","futGCMDir","fut_gcm","bc_qmap","his_gcm_bc","fut_gcm_bc"), lib = list("tidyverse","raster","qmap"))
  
  bc_data <- 1:nrow(his_obs) %>% parallel::parLapply(cl, ., function(i){
    bc_data <<- bc_qmap(df_obs    = his_obs$Climate[[i]],
                       df_his_gcm = his_gcm$Climate[[i]],
                       df_fut_gcm = fut_gcm$Climate[[i]])
    return(bc_data)
  })
  parallel::stopCluster(cl)
  his_gcm_bc$Climate <- bc_data %>% purrr::map(1) %>% purrr::map(1)
  fut_gcm_bc$Climate <- bc_data %>% purrr::map(1) %>% purrr::map(2)
  
  pDir <- paste0(root,'/data/bc_quantile_0_05deg_lat_county/',tolower(country),'/',gcm,'/1971_2000')
  if(!dir.exists(pDir)){dir.create(pDir, recursive = T)}
  fDir <- paste0(root,'/data/bc_quantile_0_05deg_lat_county/',tolower(country),'/',gcm,'/',period)
  if(!dir.exists(fDir)){dir.create(fDir, recursive = T)}
  
  if('srad' %in% colnames(his_obs$Climate[[1]])){
    outHis <- paste0(pDir,'/',tolower(county),'.RDS')
    if(!file.exists(outHis)){
      saveRDS(his_gcm_bc,outHis)
    }
    outFut <- paste0(fDir,'/',tolower(county),'.RDS')
    if(!file.exists(outFut)){
      saveRDS(fut_gcm_bc,outFut)
    }
  } else {
    outHis <- paste0(pDir,'/',tolower(county),'_prec_temp.RDS')
    if(!file.exists(outHis)){
      saveRDS(his_gcm_bc,outHis)
    }
    outFut <- paste0(fDir,'/',tolower(county),'_prec_temp.RDS')
    if(!file.exists(outFut)){
      saveRDS(fut_gcm_bc,outFut)
    }
  }
  
}
# Run once
periodList <- c('2021_2045','2041_2065')
rcpList    <- 'rcp85'
gcmList    <- c("ipsl_cm5a_mr","miroc_esm_chem","ncc_noresm1_m")
for(p in periodList){
  for(gcm in gcmList){
    BC_Qmap(country='Pakistan',county='Dadu',rcp='rcp85',gcm=gcm,period=p)
  }
}
