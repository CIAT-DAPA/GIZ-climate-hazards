library(tidyverse)
library(fst)

country <- 'Jordan'
county  <- 'jordan_rcp85'

df <- fst::read_fst('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/Jordan/past/jordan_rcp45_1980_2050_corrected.fst')
df1_ <- df %>%
  dplyr::select(year,season:ndws) %>%
  tidyr::pivot_longer(cols = 'CDD':'ndws', names_to = 'Indices', values_to = 'Value') %>%
  dplyr::group_split(Indices)

df2_ <- df %>%
  dplyr::select(year,gSeason:LGP) %>%
  tidyr::pivot_longer(cols = 'SLGP':'LGP', names_to = 'Indices', values_to = 'Value') %>%
  dplyr::group_split(Indices)

# Output folder
outDir <- paste0('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/',country,'/graphs/',tolower(county),'/time_series')
if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}

df1_ %>%
  purrr::map(.f = function(tbl){
    df_summ <- tbl %>%
      tidyr::drop_na() %>%
      dplyr::group_by(year,season) %>%
      dplyr::summarise(n      = n(),
                       mean   = mean(Value, na.rm = T),
                       sd     = sd(Value, na.rm = T)) %>%
      dplyr::mutate(sem       = sd/sqrt(n-1),
                    CI_lower  = mean + qt((1-0.95)/2, n - 1) * sem,
                    CI_upper  = mean - qt((1-0.95)/2, n - 1) * sem,
                    Serie     = ifelse(year %in% 1980:2019, 'Past','Fut'),
                    Year      = as.Date(ISOdate(year, 1, 1)),
                    season    = factor(season))
    if(length(unique(df_summ$season)) == 1){
      sem.labs <- 'S:1'
      names(sem.labs) <- 's1'
    } else {
      sem.labs <- c('S:1','S:2')
      names(sem.labs) <- c('1','2')
    }
    df_summ_ <- df_summ %>%
      dplyr::group_by(season) %>%
      dplyr::group_split(season)
    
    1:length(df_summ_) %>%
      purrr::map(.f = function(i){
        df_summ2 <- df_summ_[[i]]
        plt      <- df_summ2 %>%
          dplyr::filter(Serie == 'Past') %>%
          ggplot2::ggplot(aes(x = Year, y = mean, colour = season)) +
          ggplot2::geom_line() +
          ggplot2::scale_x_date(date_labels = "%Y", date_breaks = '10 years', limits = as.Date(c('1980-01-01', '2050-12-31'))) +
          ggplot2::ylim(min(df_summ2$mean)-5, max(df_summ2$mean)+5) +
          ggplot2::geom_line(data = df_summ2 %>% dplyr::filter(Serie == 'Fut'), aes(x = Year, y = mean, colour = season)) +
          ggplot2::geom_ribbon(data = df_summ2 %>% dplyr::filter(Serie == 'Fut'), aes(ymin = CI_lower, ymax = CI_upper, fill = season), color = "grey70", alpha = 0.4) +
          ggplot2::ylab('Average') +
          ggplot2::labs(title    = tbl$Indices %>% unique,
                        subtitle = paste0(country,", ",county),
                        caption  = "Data source: Alliance Bioversity-CIAT") +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text       = element_text(size = 17),
                         axis.title      = element_text(size = 20),
                         legend.text     = element_text(size = 17),
                         legend.title    = element_text(size = 20),
                         plot.title      = element_text(size = 20),
                         plot.subtitle   = element_text(size = 17),
                         strip.text.x    = element_text(size = 17),
                         legend.position = 'none') +
          ggplot2::facet_wrap(~season, labeller = labeller(season = sem.labs)) +
          ggplot2::ggsave(filename = paste0(outDir,'/ts_',tbl$Indices %>% unique,'_season_',i,'.png'), device = "png", width = 12, height = 6, units = "in")
      })
    
    return(cat('Graphs done\n'))
  })

df2_ %>%
  purrr::map(.f = function(tbl){
    tbl <- tbl %>%
      tidyr::drop_na()
    df_summ <- tbl %>%
      dplyr::group_by(year,gSeason) %>%
      dplyr::summarise(n      = dplyr::n(),
                       mean   = mean(Value, na.rm = T),
                       sd     = sd(Value, na.rm = T)) %>%
      dplyr::mutate(sem       = sd/sqrt(n-1),
                    CI_lower  = mean + qt((1-0.95)/2, n - 1) * sem,
                    CI_upper  = mean - qt((1-0.95)/2, n - 1) * sem,
                    Serie     = ifelse(year %in% 1980:2019, 'Past','Fut'),
                    Year      = as.Date(ISOdate(year, 1, 1)),
                    gSeason  = as.factor(gSeason))
    if(length(unique(df_summ$gSeason)) == 1){
      sem.labs <- 'S:1'
      names(sem.labs) <- '1'
    } else {
      sem.labs <- c('S:1','S:2')
      names(sem.labs) <- c('1','2')
    }
    df_summ_ <- df_summ %>%
      dplyr::group_by(gSeason) %>%
      dplyr::group_split(gSeason)
    
    1:length(df_summ_) %>%
      purrr::map(.f = function(i){
        df_summ2 <- df_summ_[[i]]
        plt      <- df_summ2 %>%
          dplyr::filter(Serie == 'Past') %>%
          ggplot2::ggplot(aes(x = Year, y = mean, colour = gSeason)) +
          ggplot2::geom_line() +
          ggplot2::scale_x_date(date_labels = "%Y", date_breaks = '10 years', limits = as.Date(c('1980-01-01', '2050-12-31'))) +
          ggplot2::ylim(min(df_summ2$mean)-5, max(df_summ2$mean)+5) +
          ggplot2::geom_line(data = df_summ2 %>% dplyr::filter(Serie == 'Fut'), aes(x = Year, y = mean, colour = gSeason)) +
          ggplot2::geom_ribbon(data = df_summ2 %>% dplyr::filter(Serie == 'Fut'), aes(ymin = CI_lower, ymax = CI_upper, fill = gSeason), color = "grey70", alpha = 0.4) +
          ggplot2::ylab('Average') +
          ggplot2::labs(title    = tbl$Indices %>% unique,
                        subtitle = paste0(country,", ",county),
                        caption  = "Data source: Alliance Bioversity-CIAT") +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text       = element_text(size = 17),
                         axis.title      = element_text(size = 20),
                         legend.text     = element_text(size = 17),
                         legend.title    = element_text(size = 20),
                         plot.title      = element_text(size = 20),
                         plot.subtitle   = element_text(size = 17),
                         strip.text.x    = element_text(size = 17),
                         legend.position = 'none') +
          ggplot2::facet_wrap(~gSeason, labeller = labeller(gSeason = sem.labs)) +
          ggplot2::ggsave(filename = paste0(outDir,'/ts_',tbl$Indices %>% unique,'_season_',i,'.png'), device = "png", width = 12, height = 6, units = "in")
      })
    
    return(cat('Graphs done\n'))
  })
