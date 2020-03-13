source("0_functions.R")

# directory
indir <- "C:\\Users\\anibi\\Documents\\ciat\\giz\\data\\"
indir <- "/share/spatial02/users/anighosh/ciat/giz/data/"
dir.create(indir, FALSE, TRUE)

# Read crop location files for GIZ
library(readxl)
m <- read_excel("data/master_list.xlsx", sheet = 1) 

# unique value chains
vc <- cleanVC(m$`Value Chain`)
vc
# filter by country
md <- m[m$Country == "Mali",]
crp <- cleanVC(md$`Value Chain`)

# Aggregate crop locations from different sources

# Source 1.1: GBIF
# install.packages("rgbif")
library(rgbif)

head(name_lookup(query = 'Rice', rank="species", return = 'data'))
key <- name_suggest(q='Oryza', rank='species')
occdf <- occ_get(key=qq[2], return='data')

# https://www.gbif.org/species/2703455
 
# Source 1.2: Crop wild relative global occurance database
# https://www.cwrdiversity.org/checklist/cwr-occurrences.php
# urls might not work in future and will need new request
u1 <- "http://api.gbif.org/v1/occurrence/download/request/0013840-200221144449610.zip"
download.file(u1, paste0(indir,"crop_wild_relative_complete.zip"), mode = "wb")

u2 <- "http://api.gbif.org/v1/occurrence/download/request/0013846-200221144449610.zip"
download.file(u2, paste0(indir,"crop_wild_relative_simple.zip"), mode = "wb")
unzip(paste0(indir,"crop_wild_relative_simple.zip"), exdir= paste0(indir,"cwr_simple"))

library(data.table)
cwr <- fread(paste0(indir, "cwr_simple/0013846-200221144449610.csv"))

# Source 2: Genesys
# devtools::install_git('https://gitlab.croptrust.org/genesys-pgr/genesysr')
# install.packages("genesysr")
# https://www.genesys-pgr.org/a/map/v2xqd5LzYpP/@0,0,3z


# use a constant set of backgorund point


# Source 3: USGS cropland data
for (i in 1:2){
  uurl <- paste0("https://api.croplands.org/data/download?page=",i,"&page_size=1000000")
  download.file(uurl, paste0(indir,"usgs30m_cropland_page_",i,".csv"), mode = "wb")
}

ff <- list.files(indir, "usgs30m_cropland_page_", full.names = TRUE)
d <- lapply(ff, read.csv, stringsAsFactors = FALSE)
d <- do.call(rbind, d)
write.csv(d, paste0(indir,"usgs30m_cropland.csv"), row.names = FALSE)
# for mali, we want rice and potato, or at least start with rice first from the USGS data?   
  
  
