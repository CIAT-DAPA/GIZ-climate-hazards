indir <- "C:\\Users\\anibi\\Documents\\ciat\\giz\\data\\"
indir <- "/share/spatial02/users/anighosh/ciat/giz/data/"

cdir <- "/share/spatial03/worldclim/wc2.1/tif/2.5m"

# models
library(sf)
library(stringr)
rc <- st_read(paste0(indir,"genesys_rice.kml"))

rc <- st_read("data/genesys_rice.kml")

# are all of them rice?
desc <- str_match(rc$Description, "<p><i>(.*?)</i>")
desc <- gsub("(<p>|<i>|</i>)","",desc)
unique(desc)

rc <- st_read("G:\\My Drive\\work\\ciat\\giz_profiles\\data\\genesys_rice.kml")
desc <- str_match(rc$Description, "<p><i>(.*?)</i>")
desc <- gsub("(<p>|<i>|</i>)","",desc)
unique(desc)
table(desc)

# worldclim current layers
library(raster)
library(dismo)
wcc <- list.files(cdir, "wc2.1_2.5m_bio_*", full.names = TRUE) 
wcr <- stack(wcc)

# extract values
presvals <- extract(wcr, rc)
# random set of background points
set.seed(0)
backgr <- randomPoints(wcr, 5000)
absvals <- extract(wcr, backgr)
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))

# save for future use
saveRDS(sdmdata, paste0(indir, "sdm.Rds"))
saveRDS(presvals, paste0(indir, "pvals.Rds"))


# model definition
# gbm.x ~ predictor variables
# gbm.y ~ response variable
fitted.model <- gbm.step(data=sdmdata, gbm.x = 2:20, gbm.y = 1,
                            family = "bernoulli", tree.complexity = 5,
                            learning.rate = 0.01, bag.fraction = 0.5)
summary(fitted.model)
