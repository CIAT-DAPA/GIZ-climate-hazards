library(tidyverse)
library(fst)
library(raster)

# Load AEZ Jordan shapefile
shp <- raster::shapefile("//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data/shps/Jordan/Jordan_wgs84_pol.shp")
vly <- shp[shp@data$Climatic_r == 'Jordan Valley',]
hgl <- shp[shp@data$Climatic_r == 'Highlands',]
dst <- shp[shp@data$Climatic_r == 'Desert',]

country <- 'Jordan'
county  <- 'jordan_rcp45'

df <- fst::read_fst(paste0('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/',country,'/past/',county,'_1980_2050_corrected.fst'))

pnt <- df %>% dplyr::select('x','y') %>% sp::SpatialPoints(coords = .)
raster::crs(pnt) <- raster::crs(shp)
vly_mch <- sp::over(pnt, vly) %>% data.frame %>% dplyr::select('Climatic_r') %>% complete.cases() %>% which()
df_vly <- df[vly_mch,]
df_vly$AEZ <- NA
df_vly$AEZ <- 'Jordan Valley'
hgl_mch <- sp::over(pnt, hgl) %>% data.frame %>% dplyr::select('Climatic_r') %>% complete.cases() %>% which()
df_hgl <- df[hgl_mch,]
df_hgl$AEZ <- NA
df_hgl$AEZ <- 'Highlands'
dst_mch <- sp::over(pnt, dst) %>% data.frame %>% dplyr::select('Climatic_r') %>% complete.cases() %>% which()
df_dst <- df[dst_mch,]
df_dst$AEZ <- NA
df_dst$AEZ <- 'Desert'

df <- dplyr::bind_rows(df_vly, df_hgl, df_dst); rm(df_vly, df_hgl, df_dst, shp, vly, vly_mch, hgl, hgl_mch, dst, dst_mch)

df1_ <- df %>%
  dplyr::select(year,AEZ,season:ndws) %>%
  tidyr::pivot_longer(cols = 'CDD':'ndws', names_to = 'Indices', values_to = 'Value') %>%
  dplyr::group_split(Indices)

df2_ <- df %>%
  dplyr::select(year,AEZ,gSeason:LGP) %>%
  tidyr::pivot_longer(cols = 'SLGP':'LGP', names_to = 'Indices', values_to = 'Value') %>%
  dplyr::group_split(Indices)

# Output folder
outDir <- paste0('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/',country,'/graphs/',tolower(county),'/time_series')
if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}

df1_ %>%
  purrr::map(.f = function(tbl){
    df_summ <- tbl %>%
      tidyr::drop_na() %>%
      dplyr::group_by(year,AEZ,season) %>%
      dplyr::summarise(n      = n(),
                       mean   = mean(Value, na.rm = T),
                       sd     = sd(Value, na.rm = T)) %>%
      dplyr::mutate(sem       = sd/sqrt(n-1),
                    CI_lower  = mean + qt((1-0.95)/2, n - 1) * sem,
                    CI_upper  = mean - qt((1-0.95)/2, n - 1) * sem,
                    Serie     = ifelse(year %in% 1980:2019, 'Past','Fut'),
                    Year      = as.Date(ISOdate(year, 1, 1)),
                    season    = factor(season))
    df_summ_ <- df_summ %>%
      dplyr::group_by(season) %>%
      dplyr::group_split(season)
    
    1:length(df_summ_) %>%
      purrr::map(.f = function(i){
        df_summ2 <- df_summ_[[i]]
        plt      <- df_summ2 %>%
          dplyr::filter(Serie == 'Past') %>%
          ggplot2::ggplot(aes(x = Year, y = mean, colour = AEZ)) +
          ggplot2::geom_line() +
          ggplot2::scale_x_date(date_labels = "%Y", date_breaks = '10 years', limits = as.Date(c('1980-01-01', '2050-12-31'))) +
          ggplot2::ylim(min(df_summ2$mean)-5, max(df_summ2$mean)+5) +
          ggplot2::geom_line(data = df_summ2 %>% dplyr::filter(Serie == 'Fut'), aes(x = Year, y = mean, colour = AEZ)) +
          ggplot2::geom_ribbon(data = df_summ2 %>% dplyr::filter(Serie == 'Fut'), aes(ymin = CI_lower, ymax = CI_upper, fill = AEZ), color = "grey70", alpha = 0.4) +
          ggplot2::ylab('Average') +
          ggplot2::labs(title    = tbl$Indices %>% unique,
                        subtitle = paste0(country,", ",county),
                        caption  = "Data source: Alliance Bioversity-CIAT") +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text       = element_text(size = 15),
                         axis.title      = element_text(size = 20),
                         legend.text     = element_text(size = 17),
                         legend.title    = element_text(size = 20),
                         plot.title      = element_text(size = 20),
                         plot.subtitle   = element_text(size = 17),
                         strip.text.x    = element_text(size = 15),
                         legend.position = 'none') +
          ggplot2::facet_wrap(~AEZ) +
          ggplot2::ggsave(filename = paste0(outDir,'/ts_',tbl$Indices %>% unique,'_season_',i,'_AEZ.png'), device = "png", width = 14, height = 6, units = "in")
      })
    
    return(cat('Graphs done\n'))
  })

df2_ %>%
  purrr::map(.f = function(tbl){
    tbl <- tbl %>%
      tidyr::drop_na()
    df_summ <- tbl %>%
      dplyr::group_by(year,AEZ,gSeason) %>%
      dplyr::summarise(n      = dplyr::n(),
                       mean   = mean(Value, na.rm = T),
                       sd     = sd(Value, na.rm = T)) %>%
      dplyr::mutate(sem       = sd/sqrt(n-1),
                    CI_lower  = mean + qt((1-0.95)/2, n - 1) * sem,
                    CI_upper  = mean - qt((1-0.95)/2, n - 1) * sem,
                    Serie     = ifelse(year %in% 1980:2019, 'Past','Fut'),
                    Year      = as.Date(ISOdate(year, 1, 1)),
                    gSeason  = as.factor(gSeason))
    df_summ_ <- df_summ %>%
      dplyr::group_by(gSeason) %>%
      dplyr::group_split(gSeason)
    
    1:length(df_summ_) %>%
      purrr::map(.f = function(i){
        df_summ2 <- df_summ_[[i]]
        plt      <- df_summ2 %>%
          dplyr::filter(Serie == 'Past') %>%
          ggplot2::ggplot(aes(x = Year, y = mean, colour = AEZ)) +
          ggplot2::geom_line() +
          ggplot2::scale_x_date(date_labels = "%Y", date_breaks = '10 years', limits = as.Date(c('1980-01-01', '2050-12-31'))) +
          ggplot2::ylim(min(df_summ2$mean)-5, max(df_summ2$mean)+5) +
          ggplot2::geom_line(data = df_summ2 %>% dplyr::filter(Serie == 'Fut'), aes(x = Year, y = mean, colour = AEZ)) +
          ggplot2::geom_ribbon(data = df_summ2 %>% dplyr::filter(Serie == 'Fut'), aes(ymin = CI_lower, ymax = CI_upper, fill = AEZ), color = "grey70", alpha = 0.4) +
          ggplot2::ylab('Average') +
          ggplot2::labs(title    = tbl$Indices %>% unique,
                        subtitle = paste0(country,", ",county),
                        caption  = "Data source: Alliance Bioversity-CIAT") +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text       = element_text(size = 15),
                         axis.title      = element_text(size = 20),
                         legend.text     = element_text(size = 17),
                         legend.title    = element_text(size = 20),
                         plot.title      = element_text(size = 20),
                         plot.subtitle   = element_text(size = 17),
                         strip.text.x    = element_text(size = 15),
                         legend.position = 'none') +
          ggplot2::facet_wrap(~AEZ) +
          ggplot2::ggsave(filename = paste0(outDir,'/ts_',tbl$Indices %>% unique,'_season_',i,'_AEZ.png'), device = "png", width = 14, height = 6, units = "in")
      })
    
    return(cat('Graphs done\n'))
  })
