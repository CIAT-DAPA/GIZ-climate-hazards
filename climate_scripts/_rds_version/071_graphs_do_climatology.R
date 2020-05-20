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

do_climatology <- function(country, county){
  
  input1 <- paste0(root, "/data/observational_data/",tolower(country),"/",tolower(county),".RDS")
  input2 <- paste0(root, "/data/observational_data/",tolower(country),"/",tolower(county),"_prec_temp.RDS")
  if(!file.exists(input1)){
    site <- readRDS(input2)
  } else {
    site <- readRDS(input1)
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
    ggplot2::scale_colour_discrete(name = "", labels = c("Tmin", "Tmax")) +
    ggplot2::annotate('rect', xmin = 5, xmax = 11, ymin = -Inf, ymax = Inf, alpha=.3, fill="forestgreen") +
    ggplot2::annotate('rect', xmin = 11, xmax = 12, ymin = -Inf, ymax = Inf, alpha=.3, fill="darkgoldenrod4") +
    ggplot2::annotate('rect', xmin = 1, xmax = 4.5, ymin = -Inf, ymax = Inf, alpha=.3, fill="darkgoldenrod4")
  ggplot2::ggsave(filename = paste0(root,'/results/',country,'/graphs/',tolower(county),'/',tolower(county),'_climatology_all_new.png'), plot = gg, device = "png", width = 12, height = 6, units = "in")
  
}

do_climatology(country = 'Pakistan', county = 'Mithi')
