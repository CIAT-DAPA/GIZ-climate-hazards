library(tidyverse)
library(fst)

country <- 'Ethiopia' # Ethiopia
county  <- 'Arsi' # Arsi

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
      dplyr::group_by(year,semester) %>%
      dplyr::summarise(n      = n(),
                       mean   = mean(Value),
                       sd     = sd(Value)) %>%
      dplyr::mutate(sem       = sd/sqrt(n-1),
                    CI_lower  = mean + qt((1-0.95)/2, n - 1) * sem,
                    CI_upper  = mean - qt((1-0.95)/2, n - 1) * sem,
                    Serie     = ifelse(year %in% 1985:2015, 'Past','Fut'),
                    Year      = as.Date(ISOdate(year, 1, 1)),
                    semester  = as.factor(semester))
    if(length(unique(df_summ$semester)) == 1){
      sem.labs <- 'S:1'
      names(sem.labs) <- '1'
    } else {
      sem.labs <- c('S:1','S:2')
      names(sem.labs) <- c('1','2')
    }
    df_summ_ <- df_summ %>%
      dplyr::group_by(semester) %>%
      dplyr::group_split(semester)
    
    1:length(df_summ_) %>%
      purrr::map(.f = function(i){
        df_summ2 <- df_summ_[[i]]
        plt      <- df_summ2 %>%
          dplyr::filter(Serie == 'Past') %>%
          ggplot2::ggplot(aes(x = Year, y = mean, colour = semester)) +
          ggplot2::geom_line() +
          ggplot2::scale_x_date(date_labels = "%Y", date_breaks = '10 years', limits = as.Date(c('1985-01-01', '2065-01-01'))) +
          ggplot2::ylim(min(df_summ2$mean)-5, max(df_summ2$mean)+5) +
          ggplot2::geom_line(data = df_summ2 %>% dplyr::filter(Serie == 'Fut'), aes(x = Year, y = mean, colour = semester)) +
          ggplot2::geom_ribbon(data = df_summ2 %>% dplyr::filter(Serie == 'Fut'), aes(ymin = CI_lower, ymax = CI_upper, fill = semester), color = "grey70", alpha = 0.4) +
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
          ggplot2::facet_wrap(~semester, labeller = labeller(semester = sem.labs)) +
          ggplot2::ggsave(filename = paste0('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/',country,'/graphs/',tolower(county),'/ts_',tbl$Indices %>% unique,'_season_',i,'.png'), device = "png", width = 12, height = 6, units = "in")
      })
    
    return(cat('Graphs done\n'))
  })
