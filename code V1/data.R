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
# require(spcosa)
require(readr)
require(gstat)
require(reshape2)
library(sf)
library(randomForest)
library(FRK)
library(tmap)
library(rgdal)
library(stars)
library(motif)
library(dplyr)
library(readr)
library(blockCV)
require(caret)
library(ncdf4)
library(rstudioapi)
library(maptools)

# setting the working directory -------------------------------------------
wrkn_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wrkn_dir)


#agb 1  ------------------------------------------------------------

files <- list.files(path = paste0(wrkn_dir, '/biomass sample/'), 
                    recursive = FALSE, pattern = "\\.tif$")
agb <- stack(paste0(wrkn_dir, '/biomass sample/', files))

agb <- agb$Mg 
# feb data (4 variables)
files <- list.files(path = paste0(wrkn_dir, '/feb/'), 
                    recursive = FALSE, pattern = "\\.tif$")
s <- stack(paste0(wrkn_dir, '/feb/', files))

# other covariates (from google)
# cannot be stacked as tehy have different extents
clay <- stack(paste0(wrkn_dir, '/data/Congo_clay_content.tif'))
cloud <- stack(paste0(wrkn_dir, '/data/cloud_cover.tif'))
sand <- stack(paste0(wrkn_dir, '/data/Congo_sand_content.tif'))
soil_carbon <- stack(paste0(wrkn_dir, '/data/Congo_soil_organic_carbon.tif'))
elevation  <- stack(paste0(wrkn_dir, '/data/Congo_elevation.tif'))

#project to use lat and lon 
agb <- projectRaster( agb$Mg,crs = crs(s[[1]]))
names(agb) <- 'ABG'

# covariates --------------------------------------------------------------
#january data 

#stack 
pstacks <- stack(agb)
for(i in 1:length(names(s))){
  print(i)
  print(names(s)[i])
  rl <- s[[i]]
  print('extract')
  rl <-crop(rl,extent(agb))
  print('croped')
  rl <-resample(rl,agb,method = 'bilinear')
  print('resampled')
  pstacks <- addLayer(pstacks, rl)
  print('processed')
}

#clay 
clay <- clay$mean_0_20
clay <-  projectRaster(clay,crs = crs(s[[1]]))
clay <- crop(clay$mean_0_20,extent(agb))
clay <-resample(clay,agb,method = 'bilinear')
names(clay) <- 'clay_content'
# clay <- crop(clay$clay_content,extent(pstacks$AGB))
pstacks <- addLayer(pstacks, clay)

#sand
sand <- sand$mean_0_20
sand <-  projectRaster(sand,crs = crs(s[[1]]))
sand <- crop(sand,extent(agb))
sand <-resample(sand,agb,method = 'bilinear')
names(sand) <- 'sand_content'
# sand <- crop(sand$sand_content,extent(pstacks$AGB))
pstacks <- addLayer(pstacks, sand)

#soil_carbon
soil_carbon <- soil_carbon$mean_0_20
soil_carbon <-  projectRaster(soil_carbon,crs = crs(s[[1]]))
soil_carbon <- crop(soil_carbon,extent(agb))
soil_carbon <-resample(soil_carbon,agb,method = 'bilinear')
names(soil_carbon) <- 'soil_carbon_content'
# soil_carbon <- crop(soil_carbon$soil_carbon_content,extent(pstacks$AGB))
pstacks <- addLayer(pstacks, soil_carbon)

#elevation 
elevation <- elevation$elevation
elevation <-  projectRaster(elevation,crs = crs(s[[1]]))
elevation <- crop(elevation,extent(agb))
elevation <-resample(elevation,agb,method = 'bilinear')
names(elevation) <- 'elevation'
# elevation <- crop(elevation$elevation,extent(pstacks$AGB))
pstacks <- addLayer(pstacks, elevation)

#cloud (not working)
cloud <- cloud$cloud_cover
cloud <-  projectRaster(cloud,crs = crs(s[[1]]))
cloud <- crop(cloud,extent(pstacks$AGB))
cloud <-resample(cloud,agb,method = 'bilinear')
names(cloud) <- 'cloud_cover'
# cloud <- crop(cloud$cloud_cover,extent(pstacks$AGB))
pstacks <- addLayer(pstacks, cloud)


# save the raster layers
pstacks[is.na(pstacks[])] <- 0
for(i in 0:length(names(pstacks))+1){
  fname <-  paste0(wrkn_dir, '/processed_data/',names(pstacks[[i]]),'.tif')
  writeRaster(pstacks[[i]], filename=fname, 
              options="INTERLEAVE=BAND", overwrite=TRUE)
}
df <- as.data.frame(pstacks,xy=T, na.rm=T)
save(df, file = paste0(wrkn_dir, '/processed_data/dataset.RData'))
rm(list=ls())




