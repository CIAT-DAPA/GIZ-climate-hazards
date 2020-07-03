# Calculate all agro-climatic indices
# A. Esquivel and H. Achicanoy
# CIAT, 2020

options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, ncdf4, sf, future, furrr, lubridate, glue, vroom, sp, fst, compiler))

OSys <- Sys.info()[1]
root <<- switch(OSys,
                'Linux'   = '/home/jovyan/work/cglabs',
                'Windows' = '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles')

# RCP 8.5
ind1 <- fst::read_fst(paste0(root,'/results/Jordan/past/jordan_rcp85_1980_2050_corrected_original.fst'))
ind2 <- fst::read_fst(paste0(root,'/results/Jordan/past/jordan_rcp85_1980_2050_gSeasons.fst'))
ind2 <- ind2 %>%
  dplyr::filter(gSeason %in% 1:2 & SLGP > 200)
ind2$gSeason <- 1

ind1 <- ind1 %>%
  dplyr::select(-gSeason,-SLGP,-LGP)
ind1 <- ind1 %>%
  unique()

all <- dplyr::left_join(x = ind1, y = ind2 %>% dplyr::select(id,year,gSeason,SLGP,LGP))

fst::write_fst(all, paste0(root,'/results/Jordan/past/jordan_rcp85_1980_2050_corrected.fst'))

# RCP 4.5
ind1 <- fst::read_fst(paste0(root,'/results/Jordan/past/jordan_rcp45_1980_2050_corrected_original.fst'))
ind2 <- fst::read_fst(paste0(root,'/results/Jordan/past/jordan_rcp45_1980_2050_gSeasons.fst'))
ind2 <- ind2 %>%
  dplyr::filter(gSeason %in% 1:2 & SLGP > 200)
ind2$gSeason <- 1

ind1 <- ind1 %>%
  dplyr::select(-gSeason,-SLGP,-LGP)
ind1 <- ind1 %>%
  unique()

all <- dplyr::left_join(x = ind1, y = ind2 %>% dplyr::select(id,year,gSeason,SLGP,LGP))

fst::write_fst(all, paste0(root,'/results/Jordan/past/jordan_rcp45_1980_2050_corrected.fst'))

