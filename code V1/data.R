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

files <- list.files(path = paste0(wrkn_dir, '/biomass sample/'), 
                    recursive = FALSE, pattern = "\\.tif$")
abg <- stack(paste0(wrkn_dir, '/biomass sample/', files))

files <- list.files(path = paste0(wrkn_dir, '/feb/'), 
                    recursive = FALSE, pattern = "\\.tif$")
s <- stack(paste0(wrkn_dir, '/feb/', files))

#project to use lat and lon 
abg <- projectRaster(abg$Mg,crs = crs(s[[1]]))

# covariates --------------------------------------------------------------
#january data 

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

# save the raster layers
for(i in 1:length(names(s))+1){
  fname <-  paste0(wrkn_dir, '/processed_data/',names(pstacks[[i]]),'.tif')
  writeRaster(pstacks[[i]], filename=fname, 
            options="INTERLEAVE=BAND", overwrite=TRUE)
}
df <- as.data.frame(pstacks,xy=T, na.rm=T)
save(df, file = paste0(wrkn_dir, '/processed_data/dataset.RData'))
ls(rm())

# Plot the data & measure correlation ----------------------------------------------------------


plot(pstacks[[1]])
df <- rasterToPoints(pstacks, spatial=TRUE)

'functions to measure the sapatial autocrealtion, plot variograms 
and general data exploration'
explore_sample <- function(){}

# sampling ---------------------------------------------------------

set.seed(123) 
#random sampling 
randomsample <- function(sp_data, sp_size){
  r_sample <- spsample(sp_data,n=sp_size, type ='random')
  valuetable <- na.omit(as.data.frame(r_sample))
  coordinates(valuetable) <- ~x+y
  plot(abg[[1]])
  plot(valuetable, add=T)
  return(r_sample)
}
# different sample sizes 
randomsample(df,400)

# claster sampling 
# 

# Machine learning models -------------------------------------------------

# cross validation --------------------------------------------------------

k_fold_cv <- function(){
  # spliting the data 
  # plot to show training and test set
  #calculate errors
}

spatial_k_fold_cv <- function(){
  # spliting the data 
  # plot to show training and test set
  # calculate errors
}

