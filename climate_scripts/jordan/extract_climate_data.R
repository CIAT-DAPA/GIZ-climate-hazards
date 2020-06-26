# Extract climate data from .grd files
# By: H. Achicanoy
# Alliance Bioversity-CIAT, 2020

options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, fst))

cat('>>> Define root path\n')
root <- 'D:/JulianCIAT'

cat('>>> Load Jordan shapefile\n')
shp <- raster::shapefile('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data/shps/Jordan/JOR_adm0.shp')

cat('>>> Extract all climate data function by RCP')
extract_data <- function(rcp = 45){
  # RCP directory
  rcp_dir <- paste0(root,'/RCP',rcp,'/OUTPUTS')
  
  cat('>>> Define list of variables')
  varList <- c('PPT','Radiation','TMAX','TMIN')
  varName <- c('prec','srad','tmax','tmin')
  
  cat('>>> Run process by year\n')
  yrList <- 1980:2050
  
  yrList %>%
    purrr::map(.f = function(year){
      out <- paste0(root,'/hazard_inputs/rcp',rcp,'/yr_',year,'.fst')
      if(!file.exists(out)){
        variables_results <- 1:length(varList) %>%
          purrr::map(.f = function(i){
            # Prepare names and dates
            target <- list.files(path = paste0(rcp_dir,'/',varList[i]), pattern = '.gri$', full.names = F)
            target <- target[target %>% grep(year, x = .)]
            output <- target %>% gsub('.gri','.grd',.)
            dayyr  <- strsplit(target, split = '_') %>% purrr::map(4) %>% unlist() %>% gsub('.gri','',.) %>% as.numeric() - 1
            dates  <- as.Date(dayyr, origin = paste0(year,'-01-01'))
            
            # Prepare and write headers for .grd files
            header  <- readr::read_lines(paste0(root,'/ExampleHeader.grd'))
            headers <- list()
            for(j in 1:length(dates)){
              headers[[j]] <- header
              if(varList[i] == 'PPT'){
                headers[[j]][length(headers[[j]])] <- strsplit(headers[[j]][length(headers[[j]])], split = '=') %>% purrr::map(1) %>% unlist() %>% paste0(., '=prec_', dates[j])
              } else {
                if(varList[i] == 'Radiation'){
                  headers[[j]][length(headers[[j]])] <- strsplit(headers[[j]][length(headers[[j]])], split = '=') %>% purrr::map(1) %>% unlist() %>% paste0(., '=srad_', dates[j])
                } else {
                  if(varList[i] == 'TMAX'){
                    headers[[j]][length(headers[[j]])] <- strsplit(headers[[j]][length(headers[[j]])], split = '=') %>% purrr::map(1) %>% unlist() %>% paste0(., '=tmax_', dates[j])
                  } else {
                    if(varList[i] == 'TMIN'){
                      headers[[j]][length(headers[[j]])] <- strsplit(headers[[j]][length(headers[[j]])], split = '=') %>% purrr::map(1) %>% unlist() %>% paste0(., '=tmin_', dates[j])
                    }
                  }
                }
              }
            }; rm(j)
            1:length(headers) %>% purrr::map(.f = function(j){readr::write_lines(headers[[j]], paste0(rcp_dir,'/',varList[i],'/',output[j]))})
            
            # Create template raster
            r <- paste0(rcp_dir,'/',varList[i],'/',output[1]) %>%
              raster::raster(.) %>%
              raster::crop(., raster::extent(shp)) %>%
              raster::mask(., mask = shp)
            
            # Obtain data.frame
            tbl <- paste0(rcp_dir,'/',varList[i],'/',output) %>%
              raster::stack(.) %>%
              raster::crop(., raster::extent(shp)) %>%
              raster::mask(., mask = shp) %>%
              raster::as.data.frame(., xy = T, na.rm = T) %>%
              tidyr::pivot_longer(., cols = 3:ncol(.), names_to = 'Date', values_to = varName[i]) %>%
              dplyr::mutate(Date = Date %>% strsplit(., split = '_') %>% purrr::map(2) %>% unlist() %>% gsub('.','-',., fixed = T) %>% as.Date)
            tbl <- tbl %>%
              dplyr::mutate(id = raster::cellFromXY(object = r, xy = tbl %>% dplyr::select(x,y) %>% as.data.frame)) %>%
              dplyr::arrange(id, Date)
            return(tbl)
          })
        variables_results[[2]] <- variables_results[[2]][,3:5]
        variables_results[[3]] <- variables_results[[3]][,3:5]
        variables_results[[4]] <- variables_results[[4]][,3:5]
        results <- variables_results %>% purrr::reduce(dplyr::left_join, by = c('id','Date'))
        results$ISO3 <- 'JOR'
        results$Country <- 'Jordan'
        results <- results %>% dplyr::select(id, x, y, ISO3, Country, Date, prec, tmax, tmin, srad)
        results$srad <- results$srad * 0.0864
        fst::write_fst(results, out)
      } else{
        cat('File exists!\n')
      }
    })
  
}
extract_data(rcp = 45)
extract_data(rcp = 85)

rcp45 <- list.files('D:/JulianCIAT/hazard_inputs/rcp45/', pattern = '.fst', full.names = T)
tbl <- rcp45 %>% purrr::map(fst::read_fst)
tbl <- dplyr::bind_rows(tbl)
fst::write_fst(tbl, 'D:/JulianCIAT/hazard_inputs/tbl_rcp45.fst')

rcp85 <- list.files('D:/JulianCIAT/hazard_inputs/rcp85/', pattern = '.fst', full.names = T)
tbl <- rcp85 %>% purrr::map(fst::read_fst)
tbl <- dplyr::bind_rows(tbl)
fst::write_fst(tbl, 'D:/JulianCIAT/hazard_inputs/tbl_rcp85.fst')
