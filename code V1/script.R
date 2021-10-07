# packages ----------------------------------------------------------------
pacman::p_load( geoR,caret,raster,sp,rgdal,pracma,
                rgeos,data.table,ggplot2,ranger,spcosa,
                readr,gstat,reshape2,sf,randomForest,FRK,tmap,
                rgdal,stars,motif,tmap,dplyr,readr,blockCV,caret,tree)

# Data loading and variables ----------------------------------------------

wrkn_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wrkn_dir)
files <- list.files(path = paste0(wrkn_dir, '/processed_data/'), 
                    recursive = FALSE, pattern = "\\.tif$")
pstacks <- stack(paste0(wrkn_dir, '/processed_data/', files))
pstacks[is.na(pstacks[])] <- 0
df <- rasterToPoints(pstacks, spatial=TRUE)
df <- as.data.frame(df)
df <- df[rowSums(df[,1:9])>0,]
df <- df[df$ABG>0,]
sp_df <- df
coordinates(sp_df) <- ~x+y 

paletteGoogleEE=c('#FFFFFF', '#CE7E45', '#DF923D', '#F1B555', '#FCD163', '#99B718', '#74A901',
                  '#66A000', '#529400', '#3E8601', '#207401', '#056201', '#004C00', '#023B01',
                  '#012E01', '#011D01', '#011301')

# functions ---------------------------------------------------------------

# raster layer plot -------------------------------------------------------
#raster plot of the Biomass raster layer 
rasterPlot <- function(r_layer, colour_palette){
  plt <-tm_shape(r_layer)+tm_raster(palette = 
                                      colour_palette) + tm_layout(
                                        legend.outside = TRUE,
                                        legend.outside.position = c("left", "bottom"), 
                                        legend.outside.size = 0.2,
                                        bg.color = "white") 
  return(plt)
}

# rasterPlot(pstacks$ABG,paletteGoogleEE)


# data sampling  ----------------------------------------------------------

#random sampling 
randomsample <- function(sp_data, sp_size, plt = FALSE){
  # r_sample <- spsample(sp_data,n=sp_size, type ='clustered')
  r_rows <- sample(nrow(sp_data), sp_size )
  r_sample <- sp_data[r_rows,]
  if(plt){
    valuetable <- as.data.frame(r_sample)
    coordinates(valuetable) <- ~x+y
    plot(pstacks[[1]], col = paletteGoogleEE)
    plot(valuetable, add=T, pch = 20, cex = 0.6)}
  
  return(r_sample)
}

# claster samples
clusterSample <- function(data,esp_value=0.15,min_pts = 10, plt = F){
  d.scale <- scale(data)
  db <- fpc::dbscan(d.scale[,c('x','y','ABG')], eps = esp_value, MinPts =min_pts)
  data$clusters <- db$cluster
  if (plt){
    # Plot DBSCAN results 
    d <- data
    coordinates(d) <- ~x+y
    plot(pstacks[[1]], col = paletteGoogleEE)
    plot(d,col=as.factor(d$clusters),add=T,pch = 20, cex = 0.6, legend =T)
  }
  return(data)
}



# Block CV (spatial cv) ---------------------------------------------------

blocks <- function(rast_layer, dat,sample_size){

  

  
  spat_bl <- spatialBlock(speciesData = dat,
                          species = NULL,
                          rasterLayer = pstacks$ABG,
                          theRange = block_range, # size of the blocks
                          # k = 20, #number of folds
                          selection = "systematic",
                          # iteration = 100, # find evenly dispersed folds
                          biomod2Format = TRUE,
                          xOffset = 0, # shift the blocks horizontally
                          yOffset = 0,
                          seed = 123,showBlocks = F)
  return(spat_bl)
  
}

sp_CV  <- function(data, sample_size){
  rs <- data
  coordinates(rs) <- ~ x + y
  gridded(rs) <- TRUE
  rs <- raster(rs)
  rs <-resample(rs,pstacks$ABG,method = 'bilinear')
  
  st_data <- st_as_sf(data, coords = c("x", "y"),
                  crs = '+proj=longlat +datum=WGS84 +no_defs',na.fail = F) 
  
  b_sizse <- spatialAutoRange(rs, doParallel = TRUE, showPlots = T, 
                              degMetre = 111325, maxpixels = 10000, 
                              plotVariograms = TRUE, progress = TRUE,
                              sampleNumber = sample_size)
  block_range <- b_sizse$range
  
  spat_bl <- spatialBlock(speciesData = st_data,
                          species = NULL,
                          rasterLayer = pstacks$ABG,
                          theRange = block_range, # size of the blocks
                          # k = length(n_folds$blocks$layer), #number of folds
                          k = 10,
                          selection = "systematic",
                          iteration = 1000, # find evenly dispersed folds
                          biomod2Format = F,
                          xOffset = 0, # shift the blocks horizontally
                          yOffset = 0,
                          seed = 123,showBlocks = T)
  
  return(list(spat_bl= spat_bl, rs = rs))
  
}

# conventional cv ---------------------------------------------------------

con_cv <- function(data,k){
  flds <- createFolds(1:nrow(d), k = k, list = T, returnTrain = F)
  
  return(flds)
}

# machine learning models -------------------------------------------------

# random forest model 
# hyperparameter tuning 

rf_hyperpaarmeter <- function(data ){
  # d <- randomsample(df, 1000)
  rf <- randomForest(ABG ~ clay_content + elevation + precipitation + 
                       sand_content + soil_carbon_content + sola_radiation + 
                       temperature + vapour + x +y, data = data, 
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
num_of_trees = rf_hyperpaarmeter(randomsample(df, 2000))
#the RF model 
rf_model <- function(data, ntree = 100, mtry = 8,nodesize =2, maxnodes= 100){
  rf <- randomForest(ABG ~ clay_content + elevation + precipitation + 
                       sand_content + soil_carbon_content + sola_radiation + 
                       temperature + vapour + x +y, data = d, 
                     ntree = ntree, mtry = mtry,nodesize =nodesize, maxnodes= maxnodes,
                     importance = T)
  return(rf)
}


# decision trees
dt_dat = randomsample(df, 2000)
stopcrit <- tree.control(nobs=nrow(dt_dat), mincut = 5, 
                         minsize = 10 , mindev = 0.003)
tree1 <- tree(ABG ~ ., data = dt_dat, control = stopcrit)
# Cost complexity pruning through CV
cv_bigtree <- cv.tree(tree1,K = 10)

plot(cv_bigtree$size, cv_bigtree$dev, type = 'b', pch = 16,
     xlab = 'Number of terminal nodes', ylab = 'CV error')
axis(side = 1, at = 1:max(cv_bigtree$size))


dt_model <- function(data,  mincut = 5, minsize = 10 , 
                     mindev = 0.003, treecut= 15){
  
  stopcrit <- tree.control(nobs=nrow(data), mincut = mincut, 
                           minsize = minsize , mindev = mindev)
  tree1 <- tree(ABG ~ ., data = data,  control = stopcrit)
  pruned_tree <- prune.tree(tree1, best = treecut)
  plot(pruned_tree)
  text(pruned_tree, cex=0.8)
  return(pruned_tree)
}


# CV output storage -------------------------------------------------------

## prepare output table
rand_samp_cv <- data.frame(model = NA, sample_size = NA,ME = NA,RMSE = NA,r2 = NA,MEC = NA)
rand_samp_scv <- data.frame(model = NA, sample_size = NA,ME = NA,RMSE = NA,r2 = NA,MEC = NA)
clus_samp_cv <- data.frame(model = NA, nclusters = NA,ME = NA,RMSE = NA,r2 = NA,MEC = NA)
clus_samp_scv <- data.frame(model = NA, nclusters = NA,ME = NA,RMSE = NA,r2 = NA,MEC = NA)


# evaluation metrics ------------------------------------------------------
errors <- function(sample_obs, pred){
  
  # mean absolute error
  ME <- mean(abs(pred - sample_obs), na.rm = TRUE)
  
  # root mean square error
  RMSE <-   sqrt(mean((pred - sample_obs)^2, na.rm = TRUE))
  
  # Pearson's correlation squared
  r2 <-  (cor(pred, sample_obs, method = 'spearman', 
              use = 'pairwise.complete.obs')^2)
  
  # coefficient of determination - 'MEC' - 'AVE'
  SSE <- sum((pred - sample_obs) ^ 2, na.rm = T)
  SST <- sum((sample_obs - mean(sample_obs, na.rm = T)) ^ 2, na.rm = T)
  MEC <- (1 - SSE/SST)
  
  return(data.frame(ME = ME, RMSE = RMSE, r2 = r2, MEC = MEC))
}



# random sampling  --------------------------------------------------------

s_sizes <- c(1000,2000)
random_sampling_data_cv = NULL
random_sampling_data_scv = NULL
k =10
for(samp in s_sizes){
  d = randomsample(df,samp,plt = TRUE)
  row.names(d) <- 1:nrow(d)
  
  # coventional cross validation
  flds <- con_cv(d,k)
  
  for(i in 1:k){
    test_set <- d[flds[[i]],]
    t <- c(1:k)[-i]
    train_set <- d[unlist(flds[t]),]
    
    valuetable_test <- na.omit(as.data.frame(test_set))
    valuetable_train <- na.omit(as.data.frame(train_set))
    coordinates(valuetable_test) <- ~x+y
    coordinates(valuetable_train) <- ~x+y
    
    plot(pstacks[[1]], col = paletteGoogleEE)
    plot(valuetable_train, add=T, col ='red', pch = 20, cex= 0.7)
    plot(valuetable_test, add=T, col ='blue',pch = 20, cex= 0.7)
    
    # random forest 
  
    rf <- rf_model(train_set)
    abg_pred <- predict(rf, newdata = test_set[,!(colnames(test_set) == "ABG")])
    
    # store 
    
    pred_test <- data.frame(x=test_set$x,y=test_set$y ,ABG =test_set$ABG, 
                            predictions = abg_pred, fold = i, 
                            model = 'random forest', samp_size =samp )
    if(i == 1){
      rf.k_fold_pred_test = pred_test
    }
    if(i > 1){
      rf.k_fold_pred_test = rbind(rf.k_fold_pred_test, pred_test)
    }
    
    
    # decision tree 
    
    dt <- dt_model(train_set)
    abg_pred <- predict(dt, newdata = test_set[,!(colnames(test_set) == "ABG")])
    
    # store 
    
    pred_test <- data.frame(x=test_set$x,y=test_set$y ,ABG =test_set$ABG, 
                            predictions = abg_pred, fold = i, model = 'decision tree' 
                            ,samp_size =samp)
    if(i == 1){
      dt.k_fold_pred_test = pred_test
    }
    if(i > 1){
      dt.k_fold_pred_test = rbind(dt.k_fold_pred_test, pred_test)
    }
    
    
  }
  
  sample_size = samp
  ME <-  errors(rf.k_fold_pred_test$ABG, rf.k_fold_pred_test$predictions)$ME
  RMSE <- errors(rf.k_fold_pred_test$ABG, rf.k_fold_pred_test$predictions)$RMSE
  r2 <-   errors(rf.k_fold_pred_test$ABG, rf.k_fold_pred_test$predictions)$r2
  MEC <-   errors(rf.k_fold_pred_test$ABG, rf.k_fold_pred_test$predictions)$MEC
  
  rand_samp_cv[nrow(rand_samp_cv) + 1,] <- c(model= 'random forest',sample_size= sample_size, ME,
                                             RMSE, r2, MEC)
  
  ME <-  errors(dt.k_fold_pred_test$ABG, dt.k_fold_pred_test$predictions)$ME
  RMSE <- errors(dt.k_fold_pred_test$ABG, dt.k_fold_pred_test$predictions)$RMSE
  r2 <-   errors(dt.k_fold_pred_test$ABG, dt.k_fold_pred_test$predictions)$r2
  MEC <-   errors(dt.k_fold_pred_test$ABG, dt.k_fold_pred_test$predictions)$MEC
  
  rand_samp_cv[nrow(rand_samp_cv) + 1,] <- c(model= 'decision tree',sample_size= sample_size, ME,
                                             RMSE, r2, MEC)
    
  
  
  # Spatial cross validation using spatila cross validation 
  
  
  spat_cv <- sp_CV(d,samp)
  flds <- spat_cv$spat_bl$folds
  for(i in 1:k){
    train_set <- d[flds[[i]][[1]],]
    test_set <- d[flds[[i]][[2]],]
    geom_plot <- spat_cv$spat_bl$plots + 
      geom_sf(data = st_as_sf(train_set, coords = c("x", "y"),
                              crs = '+proj=longlat +datum=WGS84 +no_defs',na.fail = F), alpha = 0.3, col= 'blue') +
      geom_sf(data = st_as_sf(test_set, coords = c("x", "y"),
                              crs = '+proj=longlat +datum=WGS84 +no_defs',na.fail = F), alpha = 0.3, col= 'red')
    print(geom_plot)
    
    
    # random forest tree
    
    rf <- rf_model(train_set)
    abg_pred <- predict(rf, newdata = test_set[,!(colnames(test_set) == "ABG")])
    
    # store 
    
    pred_test <- data.frame(x=test_set$x,y=test_set$y ,ABG =test_set$ABG, 
                            predictions = abg_pred, fold = i, 
                            model = 'random forest', samp_size =samp )
    if(i == 1){
      rf.k_fold_pred_test = pred_test
    }
    if(i > 1){
      rf.k_fold_pred_test = rbind(rf.k_fold_pred_test, pred_test)
    }
    
    
    # decision tree 
    
    dt <- dt_model(train_set)
    abg_pred <- predict(dt, newdata = test_set[,!(colnames(test_set) == "ABG")])
    
    # store 
    
    pred_test <- data.frame(x=test_set$x,y=test_set$y ,ABG =test_set$ABG, 
                            predictions = abg_pred, fold = i, model = 'decision tree' 
                            ,samp_size =samp)
    if(i == 1){
      dt.k_fold_pred_test = pred_test
    }
    if(i > 1){
      dt.k_fold_pred_test = rbind(dt.k_fold_pred_test, pred_test)
    }
    
  }
  
  sample_size = samp
  ME <-  errors(rf.k_fold_pred_test$ABG, rf.k_fold_pred_test$predictions)$ME
  RMSE <- errors(rf.k_fold_pred_test$ABG, rf.k_fold_pred_test$predictions)$RMSE
  r2 <-   errors(rf.k_fold_pred_test$ABG, rf.k_fold_pred_test$predictions)$r2
  MEC <-   errors(rf.k_fold_pred_test$ABG, rf.k_fold_pred_test$predictions)$MEC
  
  rand_samp_scv[nrow(rand_samp_scv) + 1,] <- c(model= 'random forest',sample_size= sample_size, ME,
                                             RMSE, r2, MEC)
  
  ME <-  errors(dt.k_fold_pred_test$ABG, dt.k_fold_pred_test$predictions)$ME
  RMSE <- errors(dt.k_fold_pred_test$ABG, dt.k_fold_pred_test$predictions)$RMSE
  r2 <-   errors(dt.k_fold_pred_test$ABG, dt.k_fold_pred_test$predictions)$r2
  MEC <-   errors(dt.k_fold_pred_test$ABG, dt.k_fold_pred_test$predictions)$MEC
  
  rand_samp_scv[nrow(rand_samp_scv) + 1,] <- c(model= 'decision tree',sample_size= sample_size, ME,
                                             RMSE, r2, MEC)
    
}

# cluster sampling --------------------------------------------------------

cl_data = clusterSample(randomsample(df,10000,plt = T),esp_value = 0.15,
                        min_pts = 10, plt = T)
clusters = unique(cl_data$clusters)
ncl <- c(20,40,60)
for (cl in ncl){
  
  d = cl_data[cl_data$clusters %in% sample(clusters,cl,replace = F),]
  
  valuetable <- d
  coordinates(valuetable) <- ~x+y
  plot(pstacks[[1]], col = paletteGoogleEE)
  plot(valuetable,col=as.factor(valuetable$clusters),add=T,pch = 20, cex = 0.8, legend =T)
  
  flds <- createFolds(1:length((d$ABG)), k = k, list = TRUE, returnTrain = FALSE)
  for(i in 1:k){
    test_set <- d[flds[[i]],]
    t <- c(1:k)[-i]
    train_set <- d[unlist(flds[t]),]
    
    valuetable_test <- na.omit(as.data.frame(test_set))
    valuetable_train <- na.omit(as.data.frame(train_set))
    coordinates(valuetable_test) <- ~x+y
    coordinates(valuetable_train) <- ~x+y
    
    plot(pstacks[[1]], col = paletteGoogleEE)
    plot(valuetable_train, add=T, col ='red', pch = 20, cex= 0.7)
    plot(valuetable_test, add=T, col ='blue',pch = 20, cex= 0.7)
    
    
    rf <- rf_model(train_set)
    abg_pred <- predict(rf, newdata = test_set[,!(colnames(test_set) == "ABG")])
    
    # store 
    
    pred_test <- data.frame(x=test_set$x,y=test_set$y ,ABG =test_set$ABG, 
                            predictions = abg_pred, fold = i )
    if(i == 1){
      k_fold_pred_test = pred_test
    }
    if(i > 1){
      k_fold_pred_test = rbind(k_fold_pred_test, pred_test)
    }
  }
  # plot
  plot(k_fold_pred_test$ABG, k_fold_pred_test$predictions, pch=19, 
       xlab="Observed ABG (in Mg/ha)", ylab="Predicted ABG (in Mg/ha)", 
       main= paste("Cluster K-fold CV (cluster size:",cl,')'))
  errors(k_fold_pred_test$ABG, k_fold_pred_test$predictions)
  abline(0,1)
  
  nclusters = cl
  ME <-  errors(k_fold_pred_test$ABG, k_fold_pred_test$predictions)$ME
  RMSE <- errors(k_fold_pred_test$ABG, k_fold_pred_test$predictions)$RMSE
  r2 <-   errors(k_fold_pred_test$ABG, k_fold_pred_test$predictions)$r2
  MEC <-   errors(k_fold_pred_test$ABG, k_fold_pred_test$predictions)$MEC
  
  clus_samp_cv[nrow(rand_samp_cv) + 1,] <- c(nclusters = cl, ME, 
                                             RMSE, r2, MEC)
  
}

    



































