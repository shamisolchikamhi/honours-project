# packages ----------------------------------------------------------------
require(geoR)
require(caret)
require(raster)
require(sp)
require(rgdal)
require(pracma)
require(rgeos)
require(data.table)
require(ggplot2)
require(ranger)
require(spcosa)
require(readr)
require(gstat)
require(reshape2)
library(sf)
library(randomForest)
library(FRK)
library(tmap)
library(rgdal)

# setting the working directory -------------------------------------------
wrkn_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wrkn_dir)


#ABG 1  ------------------------------------------------------------
#xmin, ymin, xmax, ymax
# bounding_box = c(14.417656249999986, -13.411444700107198,
#                   30.413749999999986,  -0.7914537537262674)
bounding_box = c(16.35018, -12.70012,
                 31.03487,  -4.779805)


files <- list.files(path = paste0(wrkn_dir, '/biomass sample/'), 
                    recursive = FALSE, pattern = "\\.tif$")
abg <- stack(paste0(wrkn_dir, '/biomass sample/', files))
files <- list.files(path = paste0(wrkn_dir, '/january/'), 
                    recursive = FALSE, pattern = "\\.tif$")
s <- stack(paste0(wrkn_dir, '/january/', files))

abg <- projectRaster(abg$Mg,crs = crs(s[[1]]))
# xmax(abg) <- 30
# xmin(abg) <- 10
# ymin(abg) <- -8
# ymax(abg) <- 6

# comented ----------------------------------------------------------------


# Convert raster to SpatialPointsDataFrame
# df <- rasterToPoints(abg, spatial=TRUE)
# proj4string(df)
# to lat/lon
# llprj <-  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
# llpts <- spTransform(df, CRS(llprj))
# abg1 <- raster(llpts)



# covariates --------------------------------------------------------------
#january files 

# processed_stacks <- NULL
#stack 
pstacks <- stack(abg)
for(i in 1:length(names(s))){
  print(i)
  print(names(s)[i])
  rl <- s[[i]]
  print('extract')
  rl <-crop(rl,extent(abg))
  print('croped')
  rl <-resample(rl,abg,method = 'bilinear')
  print('resampled')
  pstacks <- addLayer(pstacks, rl)
  print('processed')
}

df <- as.data.frame(pstacks,xy=T, na.rm=T)
plot(pstacks[[1]])

# Plot the data  ----------------------------------------------------------

valuetable <- na.omit(as.data.frame(df))
coordinates(valuetable) <- ~x+y

plot(headvaluetable, add=T)


abg <- stack(paste0(wrkn_dir, '/imageToDriveExample (1).tif'))
plot(abg[[1]])

df <- as.data.frame(abg)

