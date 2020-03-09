# Read crop location files for GIZ
library(readxl)
m <- read_excel("data/master_list.xlsx", sheet = 1) 

# unique value chains
vc <- cleanVC(m$`Value Chain`)

# filter by country
md <- m[m$Country == "Mali",]
crp <- cleanVC(md$`Value Chain`)

# Aggregate crop locations from different sources
# source 1: GBIF

# source 2: Genesys
# devtools::install_git('https://gitlab.croptrust.org/genesys-pgr/genesysr')

# USGS cropland data
# install.packages("genesysr")
  
  
