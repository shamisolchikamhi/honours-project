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
library(stars)
library(motif)
library(tmap)
library(dplyr)
library(readr)
library(blockCV)
require(caret)


# Plot the data & measure correlation ----------------------------------------------------------
# upload data
wrkn_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wrkn_dir)
files <- list.files(path = paste0(wrkn_dir, '/processed_data/'), 
                    recursive = FALSE, pattern = "\\.tif$")
pstacks <- stack(paste0(wrkn_dir, '/processed_data/', files))
pstacks[is.na(pstacks[])] <- 0
df <- rasterToPoints(pstacks, spatial=TRUE)
df <- as.data.frame(df)
df <- df[rowSums(df[,1:9])>0,]
# plot 
rasterPlot <- function(r_layer){
  paletteGoogleEE=c('#FFFFFF', '#CE7E45', '#DF923D', '#F1B555', '#FCD163', '#99B718', '#74A901',
                    '#66A000', '#529400', '#3E8601', '#207401', '#056201', '#004C00', '#023B01',
                    '#012E01', '#011D01', '#011301')
  plt <-tm_shape(r_layer)+tm_raster(palette = 
                                      paletteGoogleEE) + tm_layout(
                                        legend.outside = TRUE,
                                        legend.outside.position = c("left", "bottom"), 
                                        legend.outside.size = 0.2,
                                        bg.color = "white")
  return(plt)
}
rasterPlot(pstacks$ABG)

'functions to measure the sapatial autocrealtion, plot variograms 
and general data exploration'
explore_sample <- function(){}

# sampling ---------------------------------------------------------


#random sampling 
randomsample <- function(sp_data, sp_size, plt = FALSE){
  # r_sample <- spsample(sp_data,n=sp_size, type ='clustered')
  r_sample <- sample(nrow(sp_data), sp_size )
  r_sample <- df[r_sample,]
  # if(plt){
  # valuetable <- na.omit(as.data.frame(r_sample))
  # coordinates(valuetable) <- ~x+y
  # plot(pstacks[[1]])
  # plot(valuetable, add=T, pch = 20, cex = 0.8)}
  
  return(r_sample)
  
}



# cluster sampling --------------------------------------------------------

xx <- rasterToPoints(pstacks$ABG, spatial=TRUE)

biomass_signature = lsp_signature(st_as_stars(xx),type = "incove",window = 900)
#getDistMethods()
#st_dimensions
biomass_dist = lsp_to_dist(biomass_signature)
bio_hclust = hclust(biomass_dist, method = "ward.D2") #expects a distance matrix as the first argument and a linkage method as the second one
plot(bio_hclust)
clusters = cutree(eco_hclust, k = 9) #Still to decide on the number of clusters
biomass_grid_sf = lsp_add_clusters(biomass_signature,clusters)
tm_clu = tm_shape(biomass_grid_sf) +
  tm_polygons("clust", style = "cat", palette = "Set2", title = "Cluster:") +
  tm_layout(legend.position = c("LEFT", "BOTTOM"))
tm_clu

#If most clusters form continuous regions,we can merge areas of the same clusters into larger polygons.
biomass_grid_sf2 = biomass_grid_sf %>%
  dplyr::group_by(clust) %>%
  dplyr::summarize()

#Crea


# Machine learning models -------------------------------------------------
set.seed(123)
rf_hyperpaarmeter <- function(data ){
  # d <- randomsample(df, 1000)
  rf <- randomForest(ABG ~ clay_content + elevation + precipitation + 
                       sand_content + soil_carbon_content + sola_radiation + 
                       temperature + vapour + x +y, data = d, 
                     ntree = 1000, mtry = 8,nodesize =2, maxnodes= 100,
                     importance = T)
  #number of trees
  plot( rf$mse, type = 'l', xlab = 'Number of trees', ylab = 'OOB MSE', 
        col = 'navy', lwd = 2, ylim = c(0, max(rf$mse)))
  
  # variabe importance 
  varImpPlot(rf, type = 2) 
  rf_varimp <- randomForest::importance(rf, type=2)
  rf_varimp <- rf_varimp[order(rf_varimp, decreasing=FALSE),]
  par(mar = c(4,8,4,6))
  barplot(rf_varimp, horiz = T, col = 'navy', las = 1,
          xlab = 'Mean decrease in Gini index', cex.lab = 0.8, cex.axis = 0.8,
          main = 'Variable importance', cex.main = 0.8, cex.names = 0.8)
  
  return(which(rf$mse == min(rf$mse))) #number of trees
  
}
num_of_trees = rf_hyperpaarmeter(df)

rf_model <- function(data, ntree = 100, mtry = 8,nodesize =2, maxnodes= 100){
  rf <- randomForest(ABG ~ clay_content + elevation + precipitation + 
                       sand_content + soil_carbon_content + sola_radiation + 
                       temperature + vapour + x +y, data = d, 
                     ntree = ntree, mtry = mtry,nodesize =nodesize, maxnodes= maxnodes,
                     importance = T)
  return(rf)
}

# cross validation --------------------------------------------------------

# k_fold_cv <- function(){
#   # spliting the data 
#   # plot to show training and test set
#   #calculate errors
# }
# 
# spatial_k_fold_cv <- function(){
#   # spliting the data 
#   # plot to show training and test set
#   # calculate errors
# }
# 

# block_cv <- function(){

#                     BLOCKING STRATEGIES
# spatial blocking  
  y <- raster::extract(d, ABG, D = TRUE)
  b_sizse <- spatialAutoRange(pstacks, doParallel = TRUE, showPlots = TRUE, degMetre = 111325, maxpixels = 1e+05, plotVariograms = TRUE, progress = TRUE)
  block_range <- b_sizse$range
  spat_bl <- spatialBlock(speciesData = st_as_sf(rasterToPoints(pstacks$ABG, spatial = TRUE)),
                            #as(pstacks$ABG,'spatial'),
                          species = NULL,
                          rasterLayer = abg,
                          theRange = block_range, # size of the blocks
                          k = 10, #number of folds
                          selection = "random",
                          iteration = 100, # find evenly dispersed folds
                          biomod2Format = TRUE,
                          xOffset = 0, # shift the blocks horizontally
                          yOffset = 0,
                          seed = 123)
  
  
#Buffering: Generates spatially separated train and test folds by considering buffers of the specified distance around each observation point.
bf1 <- buffering(speciesdata = st_as_sf(rasterToPoints(pstacks$ABG, spatial = TRUE)),
                 species = "Species",
                 theRAnge = block_range,
                 spDataType = "PB",
                 addBG = TRUE,
                 progress = TRUE)

#Environmental Blocking
# Environmental blocking for cross-validation. This function uses clustering methods to specify sets
# of similar environmental conditions based on the input covariates. Species data corresponding to
# any of these groups or clusters are assigned to a fold. This function does the clustering in raster
# space and species data. Clustering is done using kmeans for both approaches. This function works
# on single or multiple raster files; multiple rasters need to be in a raster brick or stack format.

envBlock(
  rasterLayer = pstacks,
  speciesData = st_as_sf(rasterToPoints(pstacks$ABG, spatial = TRUE)),
  species = NULL, #null beacue response variable is continuous.
  k = 5,
  standardization = "normal",
  rasterBlock = TRUE,
  sampleNumber = 10000,
  biomod2Format = TRUE,
  numLimit = 0,
  verbose = TRUE
)


  ##                      Visualization Tools:
  foldExplorer(blocks = spat_bl,
               rasterLayer = pstacks,
               speciesData = st_as_sf(rasterToPoints(pstacks$ABG, spatial = TRUE)))

  # explore the block size
  rangeExplorer(rasterLayer = pstacks)




# cv ----------------------------------------------------------------------


s_sizes <- c(500, 1000)
k = 10
for (samp in s_sizes){
  
  errors <- NULL
  d = randomsample(df,samp)
  
  flds <- createFolds(1:length((d$ABG)), k = k, list = TRUE, returnTrain = FALSE)
  for(i in 1:k){
    
    ms_errors = NULL
    test_set <- d[flds[[i]],]
    t <- c(1:k)[-i]
    train_set <- d[unlist(flds[t]),]
    
    valuetable_test <- na.omit(as.data.frame(test_set))
    valuetable_train <- na.omit(as.data.frame(train_set))
    coordinates(valuetable_test) <- ~x+y
    coordinates(valuetable_train) <- ~x+y
    
    plot(pstacks[[1]])
    plot(valuetable_train, add=T, col ='red', pch = 20, cex= 0.9)
    plot(valuetable_test, add=T, col ='blue',pch = 20, cex= 0.9)
    
    
    rf <- rf_model(train_set)
    agb_pred <- predict(rf, newdata = test_set[,!(colnames(test_set) == "ABG")])
    mse <- mean((test_set$ABG - agb_pred)^2)
    ms_errors[i] <- mse
  }
}
# fitting models ----------------------------------------------------------