### Do climathology graph v2
### A. Esquivel, H. Achicanoy
### CIAT, 2020

# R options
options(warn = -1, scipen = 999)

# Load libraries
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, lubridate, fst))

# Paths
OSys <- Sys.info()[1]
root <<- switch(OSys,
                'Linux'   = '/home/jovyan/work/cglabs',
                'Windows' = '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles')

source(paste0(root, "/scripts/indices.R"))

do_climatology <- function(country = 'Zambia',
                           county  = 'Eastern',
                           seasons = list(s1 = c(11:12,1:4)))
{
  
  input1 <- paste0(root, "/data/observational_data/",tolower(country),"/",tolower(county),".fst")
  input2 <- paste0(root, "/data/observational_data/",tolower(country),"/",tolower(county),"_prec_temp.fst")
  if(!file.exists(input1)){
    site <- fst::read_fst(input2)
    site <- site %>%
      tidyr::nest(Climate = c('id','Date','prec','tmax','tmin')) %>%
      dplyr::rename(id = 'id1') %>%
      dplyr::select(id, everything(.))
  } else {
    site <- fst::read_fst(input1)
    site <- site %>%
      tidyr::nest(Climate = c('id','Date','prec','tmax','tmin','srad')) %>%
      dplyr::rename(id = 'id1') %>%
      dplyr::select(id, everything(.))
  }
  if(nrow(site) > 100){
    set.seed(1235)
    smpl <- sample(x = 1:nrow(site), size = 100, replace = F) %>% sort
  } else {
    smpl <- 1:nrow(site)
  }
  
  all_clmtlgy <- 1:length(smpl) %>%
    purrr::map(.f = function(i){
      clmtlgy <- site[smpl[i],] %>%
        dplyr::pull('Climate') %>%
        .[[1]] %>% 
        dplyr::mutate(Year  = lubridate::year(lubridate::as_date(Date)),
                      Month = lubridate::month(lubridate::as_date(Date))) %>%
        dplyr::group_by(Year, Month) %>%
        dplyr::summarise(Prec = sum(prec, na.rm = T),
                         Tmin = mean(tmin, na.rm = T),
                         Tmax = mean(tmax, na.rm = T)) %>%
        dplyr::group_by(Month) %>%
        dplyr::summarise(Prec = mean(Prec, na.rm = T),
                         Tmin = mean(Tmin, na.rm = T),
                         Tmax = mean(Tmax, na.rm = T))
      return(clmtlgy)
    }) %>%
    dplyr::bind_rows()
  
  avrgs <- all_clmtlgy %>%
    dplyr::group_by(Month) %>%
    dplyr::summarise(Prec = mean(Prec, na.rm = T),
                     Tmin = mean(Tmin, na.rm = T),
                     Tmax = mean(Tmax, na.rm = T))
  
  rlc <- mean(avrgs$Prec) / mean(avrgs$Tmax) * 2
  cols <- c("Tmin" = "blue", "Tmax" = "red")
  
  outDir <- paste0(root,'/results/',country,'/graphs/',tolower(county),'/')
  if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
  
  gg <- all_clmtlgy %>%
    dplyr::group_by(Month) %>%
    dplyr::summarise(Prec = mean(Prec, na.rm = T),
                     Tmin = mean(Tmin, na.rm = T),
                     Tmax = mean(Tmax, na.rm = T)) %>%
    ggplot2::ggplot(aes(x = Month, y = Prec)) +
    ggplot2::geom_bar(stat="identity", fill = 'lightblue') +
    ggplot2::xlab('Month') +
    ggplot2::ylab('Precipitation (mm)') +
    ggplot2::xlim(0,13)+
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = 1:12) +
    ggplot2::geom_line(aes(x = Month, y = Tmin*rlc, colour = 'blue'), size = 1.2) +
    ggplot2::geom_line(aes(x = Month, y = Tmax*rlc, colour = 'red'), size = 1.2) +
    #ggplot2::scale_color_identity(guide = 'legend') +
    ggplot2::theme(axis.text       = element_text(size = 17),
                   axis.title      = element_text(size = 20),
                   legend.text     = element_text(size = 17),
                   legend.title    = element_text(size = 20),
                   plot.title      = element_text(size = 20),
                   plot.subtitle   = element_text(size = 17),
                   strip.text.x    = element_text(size = 17),
                   legend.position = "top") +
    ggplot2::scale_y_continuous(sec.axis = sec_axis(~./rlc, name = 'Temperature ºC')) +
    ggplot2::scale_colour_discrete(name = "", labels = c("Tmin", "Tmax"))
  for(i in 1:length(seasons)){
    if(sum(diff(seasons[[i]]) < 0) > 0){
      gg <- gg +
        ggplot2::annotate("rect",
                          xmin  = seasons[[i]][1]-.5,
                          xmax  = seasons[[i]][which(diff(seasons[[i]]) < 0)]+.5,
                          ymin  = -Inf,
                          ymax  = Inf,
                          alpha =.3,
                          fill  = "forestgreen")
      gg <- gg +
        ggplot2::annotate("rect",
                          xmin  = seasons[[i]][which(diff(seasons[[i]]) < 0)+1]-.5,
                          xmax  = seasons[[i]][length(seasons[[i]])]+.5,
                          ymin  = -Inf,
                          ymax  = Inf,
                          alpha =.3,
                          fill  = "forestgreen")
    } else {
      gg <- gg +
        ggplot2::annotate("rect",
                          xmin  = seasons[[i]][1]-.5,
                          xmax  = seasons[[i]][length(seasons[[i]])]+.5,
                          ymin  = -Inf,
                          ymax  = Inf,
                          alpha =.3,
                          fill  = "forestgreen")
    }
  }
  ggplot2::ggsave(filename = paste0(root,'/results/',country,'/graphs/',tolower(county),'/',tolower(county),'_climatology_v3.png'), plot = gg, device = "png", width = 12, height = 6, units = "in")
  
}

do_climatology(country = 'Mozambique',
               county  = 'Manica',
               seasons = list(s1 = c(11:12,1:4)))
