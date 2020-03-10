# Raster data preparation
indir <- "C:\\Users\\anibi\\Documents\\ciat\\giz\\data\\"
dir.create(indir, FALSE, TRUE)

# current climate data @2.5m from worldclim 2.1
wurl <- "http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_bio.zip"
bczip <- file.path(indir, basename(wurl))

if(!file.exists(bczip)){
  download.file(wurl, destfile = bczip, mode = "wb")
  unzip(bczip, exdir = indir)
}

# TODO: future climate data @2.5m

# TODO: SRTM derived slope
# https://www.worldpop.org/geodata/listing?id=57

# TODO: distance to irrigation sources
# inland water: https://www.worldpop.org/project/categories?id=14
# major riverways: https://www.worldpop.org/geodata/listing?id=34

# soil
# ISRIC datahub
# global: https://files.isric.org/soilgrids/data/recent/
# africa: https://files.isric.org/public/afsis250m/

# accessibility layer
# https://www.nature.com/articles/s41597-019-0265-5#Sec10
aurl <- "https://ndownloader.figshare.com/articles/7638134/versions/3"
acczip <- file.path(indir, "accessibility.zip")

if(!file.exists(acczip)){
  download.file(aurl, destfile = acczip, mode = "wb")
  unzip(acczip, exdir = indir)
}

# population density
# worldpop global layers (2000-2020): https://www.worldpop.org/geodata/listing?id=64
# UN adjusted countries (2000-2020): https://www.worldpop.org/geodata/listing?id=69
# ftp://ftp.worldpop.org.uk/GIS/Population/Global_2000_2020/

