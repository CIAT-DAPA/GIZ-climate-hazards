library(tidyverse)
library(fst)

country <- 'Pakistan'
county  <- 'Mithi'

past    <- fst::fst(paste0("//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/",country,"/past/",county,"_1985_2015_all_new.fst")) %>% data.frame
futDir  <- paste0('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/',country,'/future')
fut_fls <- list.files(futDir, pattern = paste0('^',county,'_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_all_new.fst'), recursive = T)
fut_fls <- paste0(futDir,'/',fut_fls)
future  <- fut_fls %>%
  purrr::map(.f = function(x){df <- fst::fst(x) %>% data.frame; return(df)}) %>%
  do.call(rbind, .)

df  <- rbind(past, future)
df$season[df$season == 'Kafir'] <- 'Khariff'
df1_ <- df %>%
  dplyr::select(year,season,CDD:IRR) %>%
  tidyr::pivot_longer(cols = 'CDD':'IRR', names_to = 'Indices', values_to = 'Value') %>%
  dplyr::group_split(Indices)

# Output folder
outDir <- paste0('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/',country,'/graphs/',tolower(county),'/time_series')
if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}

df1_ %>%
  purrr::map(.f = function(tbl){
    df_summ <- tbl %>%
      dplyr::group_by(year,season) %>%
      dplyr::summarise(n      = n(),
                       mean   = mean(Value, na.rm = T),
                       sd     = sd(Value, na.rm = T)) %>%
      dplyr::mutate(sem       = sd/sqrt(n-1),
                    CI_lower  = mean + qt((1-0.95)/2, n - 1) * sem,
                    CI_upper  = mean - qt((1-0.95)/2, n - 1) * sem,
                    Serie     = ifelse(year %in% 1985:2015, 'Past','Fut'),
                    Year      = as.Date(ISOdate(year, 1, 1)),
                    season  = as.factor(season))
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
          ggplot2::scale_x_date(date_labels = "%Y", date_breaks = '10 years', limits = as.Date(c('1985-01-01', '2065-01-01'))) +
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
          ggplot2::facet_wrap(~season) +
          ggplot2::ggsave(filename = paste0(outDir,'/ts_',tbl$Indices %>% unique,'_season_',i,'_all_new.png'), device = "png", width = 12, height = 6, units = "in")
      })
    
    return(cat('Graphs done\n'))
  })
