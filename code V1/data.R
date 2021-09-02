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
for(i in 0:length(names(s))+1){
  fname <-  paste0(wrkn_dir, '/processed_data/',names(pstacks[[i]]),'.tif')
  writeRaster(pstacks[[i]], filename=fname, 
            options="INTERLEAVE=BAND", overwrite=TRUE)
}
df <- as.data.frame(pstacks,xy=T, na.rm=T)
save(df, file = paste0(wrkn_dir, '/processed_data/dataset.RData'))
rm(list=ls())

# Plot the data & measure correlation ----------------------------------------------------------
# upload data
wrkn_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wrkn_dir)
files <- list.files(path = paste0(wrkn_dir, '/processed_data/'), 
                    recursive = FALSE, pattern = "\\.tif$")
pstacks <- stack(paste0(wrkn_dir, '/processed_data/', files))
df <- rasterToPoints(pstacks, spatial=TRUE)

# plot 
paletteGoogleEE=c('#FFFFFF', '#CE7E45', '#DF923D', '#F1B555', '#FCD163', '#99B718', '#74A901',
                  '#66A000', '#529400', '#3E8601', '#207401', '#056201', '#004C00', '#023B01',
                  '#012E01', '#011D01', '#011301')
drc <-tm_shape(pstacks$Mg)+tm_raster(palette = 
                                 paletteGoogleEE) + tm_layout(
                                   legend.outside = TRUE,
                                   legend.outside.position = c("left", "bottom"), 
                                   legend.outside.size = 0.2,
                                   bg.color = "white")


drc

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
  ?plot(valuetable, add=T)
  return(r_sample)
}
# different sample sizes 
randomsample(df,100)

# claster sampling 
# 

# Machine learning models -------------------------------------------------
s_sizes <- c(500, 1000,1500)

for (samp in s_sizes){
  pt.SRS_500 <- as.data.frame(randomsample(df,samp))
  valuetable <- na.omit(as.data.frame(pt.SRS_500))
  coordinates(valuetable) <- ~x+y
  train_set <- sample(1:nrow(pt.SRS_500), 0.8*nrow(pt.SRS_500))
  df_train <- pt.SRS_500[train_set,]
  train <- df_train$ABG1
  df_test <- df[-train_set,]
  test <- df_test$ABG1
  
  # random Forrest tree 
  
  rf_tree <- randomForest(ABG1 ~ ., data = df_train, 
                          ntree = 500)
  
  pred <- predict(rf_tree, newdata = df_train, type ='response')
  rmse <- sqrt( mean((test - pred)^2))
  print(rmse)
}
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

