library(tidyverse)
library(fst)

country <- 'Pakistan' # Ethiopia
county  <- 'Muzaffargarh' # Arsi

past    <- fst(paste0("//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/",country,"/past/",county,"_1985_2015.fst")) %>% data.frame
futDir  <- paste0('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/',country,'/future')
fut_fls <- list.files(futDir, pattern = county, recursive = T)
fut_fls <- paste0(futDir,'/',fut_fls)
future  <- fut_fls %>%
  purrr::map(.f = function(x){df <- fst(x) %>% data.frame; return(df)}) %>%
  do.call(rbind, .)

df  <- rbind(past, future)
df_ <- df %>%
  tidyr::pivot_longer(cols = 'CDD':'NT35', names_to = 'Indices', values_to = 'Value') %>%
  dplyr::group_split(Indices)

df_ %>%
  purrr::map(.f = function(tbl){
    df_summ <- tbl %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(n      = n(),
                       median = median(Value),
                       mad    = mad(Value)) %>%
      dplyr::mutate(sem       = mad/sqrt(n-1),
                    CI_lower  = median + qt((1-0.95)/2, n - 1) * sem,
                    CI_upper  = median - qt((1-0.95)/2, n - 1) * sem,
                    Serie     = ifelse(year %in% 1985:2015, 'Past','Fut'),
                    Year      = as.Date(ISOdate(year, 1, 1)))
    plt <- df_summ %>%
      dplyr::filter(Serie == 'Past') %>%
      ggplot2::ggplot(aes(x = Year, y = median)) +
      ggplot2::geom_line() +
      ggplot2::scale_x_date(date_labels = "%Y", date_breaks = '5 years', limits = as.Date(c('1985-01-01', '2065-01-01'))) +
      ggplot2::ylim(min(df_summ$median)-5, max(df_summ$median)+5) +
      ggplot2::geom_line(data = df_summ %>% dplyr::filter(Serie == 'Fut'), aes(x = Year, y = median), colour = 'red') +
      ggplot2::geom_ribbon(data = df_summ %>% dplyr::filter(Serie == 'Fut'), aes(ymin = CI_lower, ymax = CI_upper), fill = "red", color = "grey70", alpha = 0.4) +
      ggplot2::ylab('Median') +
      ggplot2::labs(title    = tbl$Indices %>% unique,
                    subtitle = paste0(country,", ",county),
                    caption  = "Data source: Alliance Bioversity-CIAT") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text = element_text(size = 17),
                     axis.title = element_text(size = 20),
                     legend.text = element_text(size = 17),
                     legend.title = element_text(size = 20),
                     plot.title = element_text(size = 20),
                     plot.subtitle = element_text(size = 17)) +
      ggplot2::ggsave(filename = paste0('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/',country,'/graphs/ts_',tbl$Indices %>% unique,'.png'), device = "png", width = 12, height = 6, units = "in")
    return(plt)
  })
