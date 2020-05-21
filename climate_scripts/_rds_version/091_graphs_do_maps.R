rm(list = ls())
gc(reset = TRUE)

# =--------------------
# Packages 
library(tidyverse)
library(raster)
library(ncdf4)
library(sf)
library(future)
library(furrr)
library(lubridate)
library(glue)
library(cowsay)
library(fst)
library(ggspatial)

options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, vroom, sp, compiler))
# =--------------------

# =----------------------------------
# Identificacion de pixel para ETH
# =----------------------------------
country <- 'Pakistan'
count_i  <- c('Mithi')
adm_lvl <- 3
iso3c <- 'PAK'
Type <- 'A'

for(i in 1:length(count_i)){
  county  <- count_i[i]
  #----------------------
  C_shp <- county
  co <- tolower(country)
  country1 <- co
  # =--------------------
  
  # Ruta Principal para guardados: 
  OSys <- Sys.info()[1]
  root <<- switch(OSys,
                  'Linux'   = '/home/jovyan/work/cglabs',
                  'Windows' = '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles')
  
  # Load county shapefile
  country1 <- raster::shapefile(paste0(root,'/data/shps/',country,'/',iso3c,'_adm',adm_lvl,'.shp'))
  shp <- raster::shapefile(paste0(root,'/data/shps/',country,'/',iso3c,'_adm',adm_lvl,'.shp'))
  glue::glue('shp <- shp[shp@data$NAME_{adm_lvl} == county,]') %>%
    as.character %>%
    parse(text = .) %>%
    eval(expr = ., envir = .GlobalEnv)
  plot(shp)
  
  # Load id coords
  crd <- vroom('//dapadfs/workspace_cluster_8/climateriskprofiles/data/id_all_country.csv')
  crd <- crd %>%
    dplyr::filter(Country == country)
  pnt <- crd %>% dplyr::select(x,y) %>% sp::SpatialPoints(coords = .)
  crs(pnt) <- crs(shp)
  # Filter coordinates that are present in the county
  pnt <- sp::over(pnt, shp) %>% data.frame %>% dplyr::select(ISO) %>% complete.cases() %>% which()
  crd <- crd[pnt,]
  crd <<- crd
  
  # =------------------- Functions...
  do_climate_maps <- function(data_split){
    
    ISO3 <- unique(data_split$ISO3); county <- unique(data_split$county)
    # time <- unique(data_split$time);
    season <-  unique(data_split$season) ;  Country <- unique(data_split$Country)
    
    median_data <- data_split %>%
      group_by(id, county, Country, x, y, ISO3, season, time) %>%
      dplyr::select(-year) %>%
      summarise_all(mean) %>% ungroup() %>%
      mutate(NT35 = round(x = NT35, digits = 0))
    
    limits <- dplyr::select(median_data, CDD, P5D, P95, NT35, ndws, IRR) %>% summarise_all(.funs = c('min', 'max'))
    
    
    median_data <- median_data %>% dplyr::filter(time == 'future') %>%
      rename('CDD_f'='CDD','P5D_f'='P5D','P95_f'='P95','NT35_f'='NT35','ndws_f'='ndws','IRR_f'='IRR') %>%
      dplyr::select(-time) %>%
      dplyr::inner_join(dplyr::filter(median_data , time == 'past') %>% dplyr::select(-time), .) %>%
      dplyr::mutate(CDD_c = CDD_f-CDD,
                    P5D_c = P5D_f-P5D,
                    P95_c = P95_f-P95,
                    NT35_c = NT35_f-NT35,
                    ndws_c = ndws_f-ndws,
                    IRR_c = IRR_f-IRR)
    
    # Aqui se hace solo la figura base...
    shp_sf <- shp  %>% sf::st_as_sf()
    country <- country1 %>% sf::st_as_sf()
    xlims <- sf::st_bbox(shp_sf)[c(1, 3)]
    ylims <- sf::st_bbox(shp_sf)[c(2, 4)]
    
    b <- ggplot() +
      geom_sf(data = shp_sf, fill = 'red', color = gray(.1)) +
      geom_sf(data = country, fill = NA, color = gray(.5)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.x = element_blank(), axis.text.y = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 7, face = "bold"))
    
    ggsave(glue::glue('{path}{Country}/graphs/{tolower(county)}/maps/{county}_all_new.png'), width = 8, height = 5.5)
    
    #===---------------------------------------------------------
    #=------------------------------------------------------------
    
    # Primero dejaré hechos los de presente... luego repito los de futuro...
    # Esta función va a quedar super manual.
    index_a <- 'CDD'
    
    a <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = CDD)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_a}\n(days)'), title = 'Historic', x = 'Longitude', y = 'Latitude') +
      scale_fill_viridis_c(limits = c(limits$CDD_min, limits$CDD_max)) +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_a}_past_S{season}_all_new.png') , width = 8, height = 5.5)
    
    # =- Lo mismo para futuro...
    
    a1 <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = CDD_f)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_a}\n(days)'), title = 'Future',x = 'Longitude', y = 'Latitude') +
      scale_fill_viridis_c(limits = c(limits$CDD_min, limits$CDD_max)) +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_a}_future_S{season}_all_new.png') , width = 8, height = 5.5)
    
    a_d <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = CDD_c)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_a}\n(days)'), title = 'Change', x = 'Longitude', y = 'Latitude') +
      # scale_fill_distiller(palette = "RdYlGn") +
      scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099') +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_a}_S{season}_all_new.png') , width = 8, height = 5.5)
    
    png(filename = glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_a}_{season}_all_new.png') , width = 1580, height = 720)
    print(    gridExtra::grid.arrange(a, a1, a_d, ncol=3,  
                                      top = glue::glue('{Country}, {county}\nS:{season}',
                                                       bottom =   "Data source: Alliance Bioversity-CIAT")))
    dev.off()
    
    #===---------------------------------------------------------
    #=------------------------------------------------------------
    
    #   # =---------------------------------------
    #   # Siguiente indice...
    index_c <- 'P5D'
    
    c <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = P5D)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_c}\n(mm)'), title = 'Historic',x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits = c(limits$P5D_min, limits$P5D_max)) +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_c}_past_S{season}_all_new.png') , width = 8, height = 5.5)
    
    # =- Futuro.
    c1 <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = P5D_f)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_c}\n(mm)'), title = 'Future', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits = c(limits$P5D_min, limits$P5D_max)) +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_c}_future_S{season}_all_new.png') , width = 8, height = 5.5)
    
    c_d <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = P5D_c)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_c}\n(mm)'), title = 'Change', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099') +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_c}_S{season}_all_new.png') , width = 8, height = 5.5)
    
    
    png(filename = glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_c}_{season}_all_new.png') , width = 1580, height = 720)
    print(    gridExtra::grid.arrange(c, c1, c_d, ncol=3,  
                                      top = glue::glue('{Country}, {county}\nS:{season}',
                                                       bottom =   "Data source: Alliance Bioversity-CIAT")))
    dev.off()
    
    #===---------------------------------------------------------
    #=------------------------------------------------------------
    # =------------
    index_d <- 'P95'
    
    d <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = P95)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_d}\n(mm)'), title = 'Historic', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits = c(limits$P95_min, limits$P95_max)) +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_d}_past_S{season}_all_new.png') , width = 8, height = 5.5)
    
    # =----
    
    d1 <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = P95_f)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_d}\n(mm)'), title = 'Future', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits = c(limits$P95_min, limits$P95_max)) +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_d}_future_S{season}_all_new.png') , width = 8, height = 5.5)
    
    d_d <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = P95_c)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_d}\n(mm)'), title = 'Change', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099') +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_d}_S{season}_all_new.png') , width = 8, height = 5.5)
    
    png(filename = glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_d}_{season}_all_new.png') , width = 1580, height = 720)
    print(    gridExtra::grid.arrange(d, d1, d_d, ncol=3,  
                                      top = glue::glue('{Country}, {county}\nS:{season}',
                                                       bottom =   "Data source: Alliance Bioversity-CIAT")))
    dev.off()
    
    #===---------------------------------------------------------
    #=------------------------------------------------------------
    
    # =------------
    index_e <- 'NT35'
    
    e <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = NT35)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_e}\n(days)'), title = 'Historic', x = 'Longitude', y = 'Latitude') +
      scale_fill_viridis_c(limits = c(limits$NT35_min, limits$NT35_max)) +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_e}_past_S{season}_all_new.png') , width = 8, height = 5.5)
    
    # =------------
    
    e1 <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = NT35_f)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_e}\n(days)'), title = 'Future', x = 'Longitude', y = 'Latitude') +
      scale_fill_viridis_c(limits = c(limits$NT35_min, limits$NT35_max)) +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_e}_future_S{season}_all_new.png') , width = 8, height = 5.5)
    
    e_d <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = NT35)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_e}\n(days)'), title = 'Change', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099') +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_e}_S{season}_all_new.png') , width = 8, height = 5.5)
    
    png(filename = glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_e}_{season}_all_new.png') , width = 1580, height = 720)
    print(    gridExtra::grid.arrange(e, e1, e_d, ncol=3,  
                                      top = glue::glue('{Country}, {county}\nS:{season}',
                                                       bottom =   "Data source: Alliance Bioversity-CIAT")))
    dev.off()
    
    #===---------------------------------------------------------
    #=------------------------------------------------------------
    
    index_f <- 'ndws'
    
    f <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = ndws)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_f}\n(days)'), title = 'Historic', x = 'Longitude', y = 'Latitude') +
      scale_fill_viridis_c(limits = c(limits$ndws_min, limits$ndws_max)) +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_f}_past_S{season}_all_new.png') , width = 8, height = 5.5)
    
    # =- Lo mismo para futuro...
    
    f1 <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = ndws_f)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_f}\n(days)'), title = 'Future',x = 'Longitude', y = 'Latitude') +
      scale_fill_viridis_c(limits = c(limits$ndws_min, limits$ndws_max)) +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_f}_future_S{season}_all_new.png') ,width = 8, height = 5.5)
    
    f_d <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = ndws_c)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_f}\n(days)'), title = 'Change', x = 'Longitude', y = 'Latitude') +
      # scale_fill_distiller(palette = "RdYlGn") +
      scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099') +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_f}_S{season}_all_new.png') , width = 8, height = 5.5)
    
    png(filename = glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_f}_{season}_all_new.png') , width = 1580, height = 720)
    print(    gridExtra::grid.arrange(f, f1, f_d, ncol=3,  
                                      top = glue::glue('{Country}, {county}\nS:{season}',
                                                       bottom =   "Data source: Alliance Bioversity-CIAT")))
    dev.off()
    
    #===---------------------------------------------------------
    #=------------------------------------------------------------
    
    index_g <- 'IRR'
    
    g <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = IRR)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_g}\n(mm)'), title = 'Historic', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits = c(limits$IRR_min, limits$IRR_max)) +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_g}_past_S{season}_all_new.png'), width = 8, height = 5.5)
    
    # =----
    
    g1 <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = IRR_f)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_g}\n(mm)'), title = 'Future', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits = c(limits$IRR_min, limits$IRR_max)) +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_g}_future_S{season}_all_new.png'), width = 8, height = 5.5)
    
    g_d <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = IRR_c)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_g}\n(mm)'), title = 'Change', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099') +
      ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                        pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                        style = north_arrow_fancy_orienteering) +
      theme_bw()
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_g}_S{season}_all_new.png') , width = 8, height = 5.5)
    
    png(filename = glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_g}_{season}_all_new.png') , width = 1580, height = 720)
    print(    gridExtra::grid.arrange(g, g1, g_d, ncol=3,  
                                      top = glue::glue('{Country}, {county}\nS:{season}',
                                                       bottom =   "Data source: Alliance Bioversity-CIAT")))
    dev.off()
    
    return(median_data)}
  
  # =---------------------
  # Para correr Mali hay que cambiar los resultados por fst y hacer que quedé de forma correcta. 
  # Observed data for each country... 
  observacional_data <- readRDS(glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/data/observational_data/{co}/{county}.RDS'))
  # =--------------------
  
  path <- '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/'
  
  dir.create(glue::glue('{path}{country}/graphs/{tolower(county)}/maps'),recursive = TRUE) 
  
  if(Type == 'B'){
    past <-  fst::fst(glue::glue('{path}{country}/past/{county}_1985_2015_prec_temp.fst')) %>%
      as_tibble() %>% mutate(time = 'past') %>% unique()
    
    futDir  <- paste0('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/',country,'/future')
    fut_fls <- list.files(futDir, pattern = paste0('^',county,'_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_prec_temp.fst'), recursive = T)
    fut_fls <- paste0(futDir,'/',fut_fls)  
    
    future  <- fut_fls %>%
      purrr::map(.f = function(x){df <- fst::fst(x) %>% as_tibble() %>% mutate(time = 'future'); return(df)}) %>%
      bind_rows() %>% unique()
  } else{
    past <- fst::fst(glue::glue('{path}{country}/past/{county}_1985_2015_all_new.fst')) %>%
      tibble::as_tibble() %>% dplyr::mutate(time = 'past') %>%
      unique()
    
    futDir  <- paste0('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/',country,'/future')
    fut_fls <- list.files(futDir, pattern = paste0('^',county,'_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_all_new.fst'), recursive = T)
    fut_fls <- paste0(futDir,'/',fut_fls)  
    
    future  <- fut_fls %>%
      purrr::map(.f = function(x){df <- fst::fst(x) %>% as_tibble() %>% mutate(time = 'future'); return(df)}) %>%
      dplyr::bind_rows() %>%
      unique()
  }
  
  data_all <- dplyr::bind_rows(past, future)
  data_all$season[data_all$season == 'Kafir'] <- 'Khariff'
  
  data_all <- data_all %>% 
    dplyr::full_join(., observacional_data) %>% 
    dplyr::mutate(county = county, Country = Country) %>% 
    dplyr::select(id, ISO3, county, Country, x, y, time, season, year, CDD, P5D, P95, NT35, ndws, IRR) %>% 
    dplyr::group_split(season)
  
  data_all %>% purrr::walk(.f = do_climate_maps)
  
}
