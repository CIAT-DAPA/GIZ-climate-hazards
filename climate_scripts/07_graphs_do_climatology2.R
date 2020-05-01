library(tidyverse)
library(lubridate)
source("//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/scripts/indices.R")

do_climatology <- function(country, county, seasons = 2){
  
  input1 <- paste0("//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data/observational_data/",tolower(country),"/",tolower(county),".RDS")
  input2 <- paste0("//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data/observational_data/",tolower(country),"/",tolower(county),"_prec_temp.RDS")
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
  
  outDir <- paste0('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/',country,'/graphs/',tolower(county),'/')
  if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
  
  if(seasons == 1){
    seasonsInfo <- 1:length(smpl) %>%
      purrr::map(.f = function(i){
        info <- site[smpl[i],] %>%
          dplyr::pull('Climate') %>%
          .[[1]] %>% 
          dplyr::mutate(Year  = lubridate::year(lubridate::as_date(Date)),
                        Month = lubridate::month(lubridate::as_date(Date))) %>%
          dplyr::group_by(Year, add = T) %>%
          dplyr::group_split() %>%
          purrr::map(., .f = function(tbl){
            SummDays <- rsum.lapply(x = tbl$prec, n = 150)
            WetDays  <- SummDays[which.max(cumulative.r.sum(SummDays))]
            WetDays  <- WetDays %>% purrr::map(2) %>% unlist %>% data.frame()
            return(WetDays)
          }) %>%
          dplyr::bind_cols()
        return(info)
      }) %>%
      dplyr::bind_cols()
    tbl_summary <- data.frame(Season = 1,
                              DIni   = round(rowMeans(seasonsInfo, na.rm = T)[1]),
                              DEnd   = round(rowMeans(seasonsInfo, na.rm = T)[150])) %>%
      dplyr::mutate(Start = as.Date(DIni,'2000-01-01'),
                    End   = as.Date(DEnd,'2000-01-01'),
                    sMnth = Start %>% ymd() %>% { month(.) + day(.) / days_in_month(.) },
                    eMnth = End %>% ymd() %>% { month(.) + day(.) / days_in_month(.) })
  } else {
    if(seasons == 2){
      seasonsInfo <- 1:length(smpl) %>%
        purrr::map(.f = function(i){
          info <- site[smpl[i],] %>%
            dplyr::pull('Climate') %>%
            .[[1]] %>% 
            dplyr::mutate(Year     = lubridate::year(lubridate::as_date(Date)),
                          Month    = lubridate::month(lubridate::as_date(Date)),
                          Semester = ifelse(Month %in% 1:6, "1", "2")) %>%
            dplyr::group_by(Semester, add = T) %>%
            dplyr::group_split() %>%
            purrr::map(., .f = function(tbl){
              info2 <- tbl %>%
                dplyr::group_by(Year, add = T) %>%
                dplyr::group_split() %>%
                purrr::map(., .f = function(tbl2){
                  SummDays <- rsum.lapply(x = tbl2$prec, n = 150)
                  WetDays  <- SummDays[which.max(cumulative.r.sum(SummDays))]
                  WetDays  <- WetDays %>% purrr::map(2) %>% unlist %>% data.frame()
                  return(WetDays)
                }) %>%
                dplyr::bind_cols()
            })
          return(info)
        })
      seasonsInfo1 <- seasonsInfo %>% purrr::map(1) %>% dplyr::bind_cols()
      seasonsInfo2 <- seasonsInfo %>% purrr::map(2) %>% dplyr::bind_cols()
      tbl_summary <- data.frame(Season = c('1','2'),
                                DIni   = c(round(rowMeans(seasonsInfo1, na.rm = T)[1]),
                                           round(rowMeans(seasonsInfo2, na.rm = T)[1])+182),
                                DEnd   = c(round(rowMeans(seasonsInfo1, na.rm = T)[150]),
                                           round(rowMeans(seasonsInfo2, na.rm = T)[150])+182)) %>%
        dplyr::mutate(Start = as.Date(DIni,'2000-01-01'),
                      End   = as.Date(DEnd,'2000-01-01'),
                      sMnth = Start %>% ymd() %>% { month(.) + day(.) / days_in_month(.) },
                      eMnth = End %>% ymd() %>% { month(.) + day(.) / days_in_month(.) })
    }
  }
  
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
  if(exists('tbl_summary')){
    for(i in 1:nrow(tbl_summary)){
      gg <- gg +
        ggplot2::annotate("rect", xmin=tbl_summary$sMnth[i], xmax=tbl_summary$eMnth[i], ymin=-Inf, ymax=Inf, alpha=.3, fill="forestgreen")
    }
  }
  ggplot2::ggsave(filename = paste0('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/',country,'/graphs/',tolower(county),'/',tolower(county),'_climatology.png'), plot = gg, device = "png", width = 12, height = 6, units = "in")
  
}
# do_climatology(country = 'Pakistan', county = 'Kurram', seasons = 2)
# 
# counties <- c("Muzaffargarh",
#               "Rajan Pur",
#               "Jhang",
#               "Ghotki",
#               "Kashmore",
#               "Dadu",
#               "Mithi",
#               "Chitral",
#               "Dera Ismail Khan",
#               "South Waziristan",
#               "North Waziristan",
#               "Orakzai")
# for(cnt in counties){
  do_climatology(country = 'Pakistan', county = cnt, seasons = 2)
# }
