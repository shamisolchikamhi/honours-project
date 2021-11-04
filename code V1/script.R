# packages ----------------------------------------------------------------
pacman::p_load( geoR,caret,raster,sp,rgdal,pracma,
                rgeos,data.table,ggplot2,ranger,
                readr,gstat,reshape2,sf,randomForest,FRK,tmap,
                rgdal,stars,motif,tmap,dplyr,readr,blockCV,caret,tree)
# spcosa


# Data loading and variables ----------------------------------------------

wrkn_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wrkn_dir)
files <- list.files(path = paste0(wrkn_dir, '/processed_data/'), 
                    recursive = FALSE, pattern = "\\.tif$")
pstacks <- stack(paste0(wrkn_dir, '/processed_data/', files))

pstacks <- crop(pstacks, extent(13, 16.5 , -1, 1.5))
# pstacks <- crop(pstacks, extent(13, 15 , -2, 1.5 ))
pstacks[is.na(pstacks[])] <- 0
df <- rasterToPoints(pstacks, spatial=T)
df <- as.data.frame(df)
df <- df[rowSums(df[,1:9])>0,]
df <- df[df$ABG>0,]
sp_df <- df
coordinates(sp_df) <- ~x+y 

paletteGoogleEE=c('#FFFFFF', '#CE7E45', '#DF923D', '#F1B555', '#FCD163', '#99B718', '#74A901',
                  '#66A000', '#529400', '#3E8601', '#207401', '#056201', '#004C00', '#023B01',
                  '#012E01', '#011D01', '#011301')
brks <- c(0,1,seq(50, 450, by=25) ) 

# functions ---------------------------------------------------------------

# raster layer plot -------------------------------------------------------
#raster plot of the Biomass raster layer 
rasterPlot <- function(r_layer, colour_palette){
  plt <-tm_shape(r_layer)+tm_raster(palette = 
                                      colour_palette,breaks =brks) + tm_layout(
                                        legend.outside = TRUE,
                                        legend.outside.position = c("left", "bottom"), 
                                        legend.outside.size = 0.2,
                                        bg.color = "white") 
  return(plt)
}

rasterPlot(pstacks$ABG,paletteGoogleEE)


# data sampling  ----------------------------------------------------------

#random sampling 
randomsample <- function(sp_data, sp_size, plt = FALSE){
  rm(.Random.seed, envir=globalenv())
  r_rows <- sample(nrow(sp_data), sp_size, replace = F )
  r_sample <- sp_data[r_rows,]
  if(plt){
    valuetable <- as.data.frame(r_sample)
    coordinates(valuetable) <- ~x+y
    plot(pstacks[[1]], col = paletteGoogleEE,
       xlim = c(extent(pstacks[[1]])[1], extent(pstacks[[1]])[2]),
       ylim = c(extent(pstacks[[1]])[3], extent(pstacks[[1]])[4]),
         axes = F, box = F, legend=T) 
    plot(valuetable, add=T, pch = 20, cex = 0.3, col = 'red')}
  
  return(r_sample)
}

# claster samples
clusterSample <- function(data,esp_value=0.15,min_pts = 10, plt = F){
  d.scale <- scale(data)
  db <- fpc::dbscan(d.scale[,c('x','y','ABG')], eps = esp_value, 
                    MinPts =min_pts,method = 'hybrid')
  data$clusters <- db$cluster
  if (plt){
    # Plot DBSCAN results 
    d <- data
    coordinates(d) <- ~x+y
    plot(pstacks[[1]], col = paletteGoogleEE,
         xlim = c(extent(pstacks[[1]])[1], extent(pstacks[[1]])[2]), 
         ylim = c(extent(pstacks[[1]])[3], extent(pstacks[[1]])[4]), 
         axes = FALSE, box = FALSE, legend=T)
    plot(d,col=as.factor(d$clusters),add=T,pch = 20, cex = 0.3, legend =T)
  }
  return(data)
}


# Block CV (spatial cv) ---------------------------------------------------

blocks <- function(dat,sz, bk_range){
  
  spat_bl <- spatialBlock(speciesData = dat,
                          species = NULL,
                          rasterLayer = pstacks$ABG,
                          theRange = bk_range, # size of the blocks
                          k = 2, #number of folds
                          selection = "systematic",
                          # iteration = 100, # find evenly dispersed folds
                          biomod2Format = TRUE,
                          xOffset = 0, # shift the blocks horizontally
                          yOffset = 0,
                          seed = 123,showBlocks = F)
  return(spat_bl)
  
}

sp_CV  <- function(data, sample_size=500){
  
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
  
  n_folds <- blocks(dat = st_data, sz = sample_size, bk_range=block_range)
  
  if(length(n_folds$blocks$layer)>= 10){
    
    f = 10
    
  }else{
    f = length(n_folds$blocks$layer)
  }
  
  print(f)
  spat_bl <- spatialBlock(speciesData = st_data,
                          species = NULL,
                          rasterLayer = pstacks$ABG,
                          theRange = block_range, # size of the blocks
                          # k = length(n_folds$blocks$layer), #number of folds
                          k = f,
                          selection = "systematic",
                          iteration = 100, # find evenly dispersed folds
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
                       temperature + vapour, data = data, 
                     ntree = 500, mtry = 6,nodesize = 10, maxnodes= 14,
                     importance = T)
  #number of trees
  plot( rf$mse, type = 'l', xlab = 'Number of trees', ylab = 'OOB MSE', 
        col = 'navy', lwd = 2, ylim = c(min(rf$mse)-10, max(rf$mse)))
  
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
num_of_trees = rf_hyperpaarmeter(randomsample(df, 5000, plt = T))
#the RF model 
rf_model <- function(data, ntree = 150, mtry = 6,nodesize =10, maxnodes= 14){
  rf <- randomForest(ABG ~ clay_content + elevation + precipitation + 
                       sand_content + soil_carbon_content + sola_radiation + 
                       temperature + vapour, data = d, 
                     ntree = ntree, mtry = mtry,nodesize =nodesize, maxnodes= maxnodes,
                     importance = T)
  return(rf)
}


# decision trees
dt_dat = randomsample(df, 5000)
stopcrit <- tree.control(nobs=nrow(dt_dat), 
                         minsize = 2 , mindev = 0.003)
tree1 <- tree(ABG ~ ., data = dt_dat, control = stopcrit)
# Cost complexity pruning through CV
cv_bigtree <- cv.tree(tree1,K = 10)

plot(cv_bigtree$size, cv_bigtree$dev, type = 'b', pch = 16,
     xlab = 'Number of terminal nodes', ylab = 'CV error')
axis(side = 1, at = 1:max(cv_bigtree$size))

dt_model <- function(data, minsize = 2, 
                     mindev = 0.003, treecut= 14){
  
  stopcrit <- tree.control(nobs=nrow(data), minsize = minsize , mindev = mindev)
  tree1 <- tree(ABG ~ clay_content + elevation + precipitation + 
                  sand_content + soil_carbon_content + sola_radiation + 
                  temperature + vapour, data = data,  control = stopcrit)
  pruned_tree <- prune.tree(tree1, best = treecut)
  # pruned_tree <- tree1
  plot(pruned_tree)
  text(pruned_tree, cex=0.8)
  return(pruned_tree)
}


# CV output storage -------------------------------------------------------

## prepare output table
rand_samp_cv <- data.frame(model = NA, sample_size = 0,range = 0,sill =0,ME = 0,RMSE = 0,r2 = 0,MEC = 0)
rand_samp_scv <- data.frame(model = NA, sample_size = 0,range = 0,sill =0,ME = 0,RMSE = 0,r2 = 0,MEC = 0)
clus_samp_cv <- data.frame(model = NA, sample_size = 0,range = 0,sill =0,ME = 0,RMSE = 0,r2 = 0,MEC = 0)
clus_samp_scv <- data.frame(model = NA, sample_size = 0,range = 0,sill =0,ME = 0,RMSE = 0,r2 = 0,MEC = 0)

rf.k_fold_pred <- data.frame()
dt.k_fold_pred <- data.frame()
rf.c_fold_pred <- data.frame()
dt.c_fold_pred <- data.frame()

rf.k_fold_pred2 <- data.frame()
dt.k_fold_pred2 <- data.frame()
rf.c_fold_pred2 <- data.frame()
dt.c_fold_pred2 <- data.frame()

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

s_sizes <- c(500,1000,2000, 5000,10000, 20000)
# s_sizes <-  c(10000,15000,20000,25000)
random_sampling_data_cv = NULL
random_sampling_data_scv = NULL
k =10
niter <- c(1:50)
for(i in niter){
  print(paste0('iteration number ********************', i))
for(samp in s_sizes){
  d = randomsample(df,samp,plt = TRUE)
  row.names(d) <- 1:nrow(d)
  
  rst <- d
  coordinates(rst) <- ~ x + y
  gridded(rst) <- TRUE
  rst <- raster(rst)
  
  spat_auto <- spatialAutoRange(rst, doParallel = TRUE, showPlots = T,
                   degMetre = 111325, maxpixels = 10000,
                   plotVariograms = TRUE, progress = TRUE,
                   sampleNumber = samp)
  print('auto cor done')
  # coventional cross validation
  flds <- con_cv(d,k)
  print('con cv')
  for(i in 1:k){
    test_set <- d[flds[[i]],]
    t <- c(1:k)[-i]
    train_set <- d[unlist(flds[t]),]
    
    valuetable_test <- na.omit(as.data.frame(test_set))
    valuetable_train <- na.omit(as.data.frame(train_set))
    coordinates(valuetable_test) <- ~x+y
    coordinates(valuetable_train) <- ~x+y
    
    plot(pstacks[[1]], col = paletteGoogleEE, axes = FALSE, box = FALSE, legend=T)
    plot(valuetable_train, add=T, col ='red', pch = 20, cex= 0.4)
    plot(valuetable_test, add=T, col ='yellow',pch = 20, cex= 0.4)
    
    # random forest 
    
    rf <- rf_model(train_set)
    abg_pred <- predict(rf, newdata = test_set[,!(colnames(test_set) == "ABG")])
    
    # store 
    
    pred_test <- data.frame(x=test_set$x,y=test_set$y ,ABG =test_set$ABG, 
                            predictions = abg_pred, fold = i, 
                            model = 'random forest', samp_size =samp,range = spat_auto$variograms$var_model$range[2], 
                            sill= spat_auto$variograms$var_model$psill[2])
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
                            ,samp_size =samp,range = spat_auto$variograms$var_model$range[2], 
                            sill= spat_auto$variograms$var_model$psill[2])
    if(i ==1){
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
  
  rand_samp_cv[nrow(rand_samp_cv) + 1,] <- c(model= 'random forest',sample_size= sample_size,range = spat_auto$variograms$var_model$range[2], 
                                             sill= spat_auto$variograms$var_model$psill[2], ME, RMSE, r2, MEC)
  
  ME <-  errors(dt.k_fold_pred_test$ABG, dt.k_fold_pred_test$predictions)$ME
  RMSE <- errors(dt.k_fold_pred_test$ABG, dt.k_fold_pred_test$predictions)$RMSE
  r2 <-   errors(dt.k_fold_pred_test$ABG, dt.k_fold_pred_test$predictions)$r2
  MEC <-   errors(dt.k_fold_pred_test$ABG, dt.k_fold_pred_test$predictions)$MEC
  
  rand_samp_cv[nrow(rand_samp_cv) + 1,] <- c(model= 'decision tree',sample_size= sample_size,range = spat_auto$variograms$var_model$range[2], 
                                             sill= spat_auto$variograms$var_model$psill[2], ME,RMSE, r2, MEC)
  
  
  print('spat cv')
  # Spatial cross validation using spatila cross validation 
  try(spat_cv <- sp_CV(d,samp))
  flds <- spat_cv$spat_bl$folds
  for(i in 1:min(k, length(flds))){
    train_set <- d[flds[[i]][[1]],]
    test_set <- d[flds[[i]][[2]],]
    geom_plot <- spat_cv$spat_bl$plots + geom_density(color = paletteGoogleEE) +
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
                            model = 'random forest', samp_size =samp,range = spat_auto$variograms$var_model$range[2], 
                            sill= spat_auto$variograms$var_model$psill[2])
    if(i == 1){
      rf.c_fold_pred_test = pred_test
    }
    if(i > 1){
      rf.c_fold_pred_test = rbind(rf.c_fold_pred_test, pred_test)
    }
    
    # decision tree 
    
    dt <- dt_model(train_set)
    abg_pred <- predict(dt, newdata = test_set[,!(colnames(test_set) == "ABG")])
    
    # store 
    
    pred_test <- data.frame(x=test_set$x,y=test_set$y ,ABG =test_set$ABG, 
                            predictions = abg_pred, fold = i, model = 'decision tree' 
                            ,samp_size =samp,range = spat_auto$variograms$var_model$range[2], 
                            sill= spat_auto$variograms$var_model$psill[2])
    if(i == 1){
      dt.c_fold_pred_test = pred_test
    }
    if(i > 1){
      dt.c_fold_pred_test = rbind(dt.c_fold_pred_test, pred_test)
    }
  
  }
  
  sample_size = samp
  ME <-  errors(rf.c_fold_pred_test$ABG, rf.c_fold_pred_test$predictions)$ME
  RMSE <- errors(rf.c_fold_pred_test$ABG, rf.c_fold_pred_test$predictions)$RMSE
  r2 <-   errors(rf.c_fold_pred_test$ABG, rf.c_fold_pred_test$predictions)$r2
  MEC <-   errors(rf.c_fold_pred_test$ABG, rf.c_fold_pred_test$predictions)$MEC
  
  rand_samp_scv[nrow(rand_samp_scv) + 1,] <- c(model= 'random forest',sample_size= sample_size,range = spat_auto$variograms$var_model$range[2], 
                                               sill= spat_auto$variograms$var_model$psill[2], ME,
                                               RMSE, r2, MEC)
  
  ME <-  errors(dt.c_fold_pred_test$ABG, dt.c_fold_pred_test$predictions)$ME
  RMSE <- errors(dt.c_fold_pred_test$ABG, dt.c_fold_pred_test$predictions)$RMSE
  r2 <-   errors(dt.c_fold_pred_test$ABG, dt.c_fold_pred_test$predictions)$r2
  MEC <-   errors(dt.c_fold_pred_test$ABG, dt.c_fold_pred_test$predictions)$MEC
  
  rand_samp_scv[nrow(rand_samp_scv) + 1,] <- c(model= 'decision tree',sample_size= sample_size,range = spat_auto$variograms$var_model$range[2], 
                                               sill= spat_auto$variograms$var_model$psill[2], ME,
                                               RMSE, r2, MEC)
  
  write.csv(rand_samp_scv, 'rand_samp_scv5.csv', row.names = F, append = T)
  write.csv(rand_samp_cv, 'rand_samp_cv5.csv', row.names = F, append = T)
  
  rf.k_fold_pred <- rbind(rf.k_fold_pred, rf.k_fold_pred_test)
  write.csv(rf.k_fold_pred, 'rf.k_fold_pred.csv',row.names = F, append = T )
  
  dt.k_fold_pred <- rbind(dt.k_fold_pred, dt.k_fold_pred_test)
  write.csv(dt.k_fold_pred, 'dt.k_fold_pred.csv',row.names = F, append = T )
  
  rf.c_fold_pred <- rbind(rf.c_fold_pred, rf.c_fold_pred_test)
  write.csv(rf.c_fold_pred, 'rf.c_fold_pred.csv',row.names = F, append = T )
  dt.c_fold_pred <- rbind(dt.c_fold_pred, dt.c_fold_pred_test)
  write.csv(dt.c_fold_pred, 'dt.c_fold_pred.csv',row.names = F, append = T )
  
} 


}


# cluster sampling --------------------------------------------------------

for(ii in niter){
  
  for(sp in c(10000, 15000, 20000,1000,2000,5000)){
    cl_data = clusterSample(randomsample(df, sp,plt = T),esp_value = 0.1,
                            min_pts = sp/1000, plt = T)

    clust = unique(cl_data$clusters)
for (cl in clust){
  
  d = cl_data[cl_data$clusters==cl,]
  samp <- nrow(d)
  samp <- round(samp, digits = 3)
  row.names(d) <- 1:nrow(d)
  
  if(nrow(d) >= 400){
  
  rst <- d
  coordinates(rst) <- ~ x + y
  gridded(rst) <- TRUE
  rst <- raster(rst)
  
  spat_auto <- spatialAutoRange(rst, doParallel = TRUE, showPlots = T,
                                degMetre = 111325, maxpixels = 10000,
                                plotVariograms = TRUE, progress = TRUE,
                                sampleNumber = samp)
  print('auto cor done')
  # coventional cross validation
 
  
  valuetable <- d
  coordinates(valuetable) <- ~x+y
  plot(pstacks[[1]], col = paletteGoogleEE)
  plot(valuetable,col=as.factor(valuetable$clusters),add=T,pch = 20, cex = 0.8, legend =T)
  
  flds <- con_cv(d,k)
  print('con cv')
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
                            model = 'random forest', samp_size =samp,range = spat_auto$variograms$var_model$range[2], 
                            sill= spat_auto$variograms$var_model$psill[2])
    if(i == 1){
      rf.k_fold_pred2_test = pred_test
    }
    if(i > 1){
      rf.k_fold_pred2_test = rbind(rf.k_fold_pred2_test, pred_test)
    }
    
    # decision tree 
    
    dt <- dt_model(train_set)
    abg_pred <- predict(dt, newdata = test_set[,!(colnames(test_set) == "ABG")])
    
    # store 
    
    pred_test <- data.frame(x=test_set$x,y=test_set$y ,ABG =test_set$ABG, 
                            predictions = abg_pred, fold = i, model = 'decision tree' 
                            ,samp_size =samp,range = spat_auto$variograms$var_model$range[2], 
                            sill= spat_auto$variograms$var_model$psill[2])
    if(i == 1){
      dt.k_fold_pred2_test = pred_test
    }
    if(i > 1){
      dt.k_fold_pred2_test = rbind(dt.k_fold_pred2_test, pred_test)
    }
    
  }
  sample_size = samp
  ME <-  errors(rf.k_fold_pred2_test$ABG, rf.k_fold_pred2_test$predictions)$ME
  RMSE <- errors(rf.k_fold_pred2_test$ABG, rf.k_fold_pred2_test$predictions)$RMSE
  r2 <-   errors(rf.k_fold_pred2_test$ABG, rf.k_fold_pred2_test$predictions)$r2
  MEC <-   errors(rf.k_fold_pred2_test$ABG, rf.k_fold_pred2_test$predictions)$MEC
  
  clus_samp_cv[nrow(clus_samp_cv) + 1,] <- c(model= 'random forest',sample_size= sample_size,range = spat_auto$variograms$var_model$range[2], 
                                             sill= spat_auto$variograms$var_model$psill[2], ME,
                                             RMSE, r2, MEC)
  
  ME <-  errors(dt.k_fold_pred2_test$ABG, dt.k_fold_pred2_test$predictions)$ME
  RMSE <- errors(dt.k_fold_pred2_test$ABG, dt.k_fold_pred2_test$predictions)$RMSE
  r2 <-   errors(dt.k_fold_pred2_test$ABG, dt.k_fold_pred2_test$predictions)$r2
  MEC <-   errors(dt.k_fold_pred2_test$ABG, dt.k_fold_pred2_test$predictions)$MEC
  
  clus_samp_cv[nrow(clus_samp_cv) + 1,] <- c(model= 'decision tree',sample_size= sample_size,range = spat_auto$variograms$var_model$range[2], 
                                             sill= spat_auto$variograms$var_model$psill[2], ME,
                                             RMSE, r2, MEC)
  
  
  print('spat cv')
  
  # Spatial cross validation using spatila cross validation 
  
  
  try(spat_cv <- sp_CV(d,samp))
  flds <- spat_cv$spat_bl$folds
  for(i in 1:min(k, length(flds))){
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
                            model = 'random forest', samp_size =samp,range = spat_auto$variograms$var_model$range[2], 
                            sill= spat_auto$variograms$var_model$psill[2])
    if(i == 1){
      rf.c_fold_pred2_test = pred_test
    }
    if(i > 1){
      rf.c_fold_pred2_test = rbind(rf.c_fold_pred2_test, pred_test)
    }
    
    # decision tree 
    
    dt <- dt_model(train_set)
    abg_pred <- predict(dt, newdata = test_set[,!(colnames(test_set) == "ABG")])
    
    # store 
    
    pred_test <- data.frame(x=test_set$x,y=test_set$y ,ABG =test_set$ABG, 
                            predictions = abg_pred, fold = i, model = 'decision tree' 
                            ,samp_size =samp,range = spat_auto$variograms$var_model$range[2], 
                            sill= spat_auto$variograms$var_model$psill[2])
    if(i == 1){
      dt.c_fold_pred2_test = pred_test
    }
    if(i > 1){
      dt.c_fold_pred2_test = rbind(dt.c_fold_pred2_test, pred_test)
    }
    
  }
  
  sample_size = samp
  ME <-  errors(rf.c_fold_pred2_test$ABG, rf.c_fold_pred2_test$predictions)$ME
  RMSE <- errors(rf.c_fold_pred2_test$ABG, rf.c_fold_pred2_test$predictions)$RMSE
  r2 <-   errors(rf.c_fold_pred2_test$ABG, rf.c_fold_pred2_test$predictions)$r2
  MEC <-   errors(rf.c_fold_pred2_test$ABG, rf.c_fold_pred2_test$predictions)$MEC
  
  clus_samp_scv[nrow(clus_samp_scv) + 1,] <- c(model= 'random forest',sample_size= sample_size,range = spat_auto$variograms$var_model$range[2], 
                                               sill= spat_auto$variograms$var_model$psill[2], ME,
                                               RMSE, r2, MEC)
  
  ME <-  errors(dt.c_fold_pred2_test$ABG, dt.c_fold_pred2_test$predictions)$ME
  RMSE <- errors(dt.c_fold_pred2_test$ABG, dt.c_fold_pred2_test$predictions)$RMSE
  r2 <-   errors(dt.c_fold_pred2_test$ABG, dt.c_fold_pred2_test$predictions)$r2
  MEC <-   errors(dt.c_fold_pred2_test$ABG, dt.c_fold_pred2_test$predictions)$MEC
  
  clus_samp_scv[nrow(clus_samp_scv) + 1,] <- c(model= 'decision tree',sample_size= sample_size,range = spat_auto$variograms$var_model$range[2], 
                                               sill= spat_auto$variograms$var_model$psill[2], ME,
                                               RMSE, r2, MEC)
  
  write.csv(clus_samp_scv, 'clus_samp_scv5.csv',row.names = F, append = T)
  write.csv(clus_samp_cv, 'clus_samp_cv5.csv',row.names = F, append = T)
  
  dt.c_fold_pred2 <- rbind(dt.c_fold_pred2, dt.c_fold_pred2_test)
  write.csv(dt.c_fold_pred2, 'dt.c_fold_pred2.csv',row.names = F, append = T )
  rf.c_fold_pred2 <- rbind(rf.c_fold_pred2, rf.c_fold_pred2_test)
  write.csv(rf.c_fold_pred2, 'rf.c_fold_pred2.csv',row.names = F, append = T )
  
  rf.k_fold_pred2 <- rbind(rf.k_fold_pred2, rf.k_fold_pred2_test)
  write.csv(rf.k_fold_pred2, 'rf.k_fold_pred2.csv',row.names = F, append = T )
  dt.k_fold_pred2 <- rbind(dt.k_fold_pred2, dt.k_fold_pred2_test)
  write.csv(dt.k_fold_pred2, 'dt.k_fold_pred2.csv',row.names = F, append = T )
  
}
  
}
}

}



# Results -----------------------------------------------------------------

###### cluster data #####
clus_samp_cv <- read_csv("clus_samp_cv5.csv")
clus_samp_scv <- read_csv("clus_samp_scv5.csv")
dt_k_fold_pred2 <- read_csv("dt.k_fold_pred2.csv")
rf_k_fold_pred2 <- read_csv("rf.k_fold_pred2.csv")
dt_c_fold_pred2 <- read_csv("dt.c_fold_pred2.csv")
rf_c_fold_pred2 <- read_csv("rf.c_fold_pred2.csv")

clus_samp_cv <- clus_samp_cv[!is.na(clus_samp_cv$model), ]
clus_samp_cv$cv_type <- 'cv'
clus_samp_scv <- clus_samp_scv[!is.na(clus_samp_scv$model), ]
clus_samp_scv$cv_type <- 'scv'
clust_sample <- rbind(as.data.frame(clus_samp_scv), as.data.frame(clus_samp_cv))
# clust_sample <- rbind(clust_sample,clust_sample2)
clust_sample$size[clust_sample$sample_size <= 1500] <- 'small'    
clust_sample$size[clust_sample$sample_size <= 5000 & clust_sample$sample_size > 1500] <- 'medium'   
clust_sample$size[clust_sample$sample_size > 10000 & clust_sample$sample_size > 5000] <- 'large'   
# clust_sample$size[clust_sample$sample_size > 10000] <- 'x-large' 


#####random sample data ####

rand_samp_cv <- read_csv("rand_samp_cv5.csv")
rand_samp_scv <- read_csv("rand_samp_scv5.csv")
dt_k_fold_pred <- read_csv("dt.k_fold_pred.csv")
rf_k_fold_pred <- read_csv("rf.k_fold_pred.csv")
dt_c_fold_pred <- read_csv("dt.c_fold_pred.csv")
rf_c_fold_pred <- read_csv("rf.c_fold_pred.csv")

rand_samp_cv <- rand_samp_cv[!is.na(rand_samp_cv$model), ]
rand_samp_cv$cv_type <- 'cv'
rand_samp_scv <- rand_samp_scv[!is.na(rand_samp_scv$model), ]
rand_samp_scv$cv_type <- 'scv'
rand_sample <- rbind(as.data.frame(rand_samp_scv), as.data.frame(rand_samp_cv))


# cluster sample  ---------------------------------------------------------

ggplot(data =clust_sample, aes(x= factor(size, levels = c('small','medium', 'large', 'x-large')), y=RMSE)) +
  geom_boxplot(aes(fill=paste0(model,'-',cv_type))) + ggtitle('Cluster sample - sample sizes') + 
  xlab('sample size') + guides(fill=guide_legend("model-cv_type"))

ggplot(data =clust_sample[clust_sample$MEC > 0,], aes(x= factor(size, levels = c('small','medium', 'large', 'x-large')), y= MEC)) + 
  geom_boxplot(aes(fill=paste0(model,'-',cv_type))) +ggtitle('Cluster sample - sample sizes') + 
  xlab('sample size') + ylab('propotion explained') + guides(fill=guide_legend("model-cv_type"))

ggplot(data =clust_sample, aes(x= factor(size, levels = c('small','medium', 'large', 'x-large')), y=r2)) +
  geom_boxplot(aes(fill=paste0(model,'-',cv_type))) + ggtitle('Cluster sample - sample sizes') + 
  xlab('sample size') + guides(fill=guide_legend("model-cv_type"))
ggplot(data =clust_sample, aes(x= factor(size, levels = c('small','medium', 'large', 'x-large')), y=ME)) +
  geom_boxplot(aes(fill=paste0(model,'-',cv_type))) + ggtitle('Cluster sample - sample sizes') + 
  xlab('sample size') + guides(fill=guide_legend("model-cv_type"))




# semi-variance 
ggplot(data = clust_sample[clust_sample$sill <= 5000,], aes(x= factor(size, levels = c('small','medium', 'large', 'x-large')), y=sill)) + 
  geom_boxplot(fill = 'grey') +xlab('sample size') + ylab('semi_variance') + ggtitle('Cluster samples')

ggplot(data = clust_sample[clust_sample$sill <= 2000,], aes(x=sill, y=RMSE)) + 
  geom_line(aes(colour=as.factor(paste0(model,'-',cv_type)))) +xlab('semi_variance') +
  ylab('RMSE') + guides(col=guide_legend("model-cv_type")) + ggtitle('Cluster samples')

ggplot(data = clust_sample[clust_sample$sill <= 2000 & clust_sample$MEC >0,], aes(x=sill, y=MEC)) + 
  geom_line(aes(colour=as.factor(paste0(model,'-',cv_type)))) +xlab('semi_variance') +
  ylab('propotion explained') + guides(col=guide_legend("model-cv_type")) + ggtitle('Cluster samples')

ggplot(data = clust_sample[clust_sample$sill <= 2000,], aes(x=sill, y=r2)) + 
  geom_line(aes(colour=as.factor(paste0(model,'-',cv_type)))) +xlab('semi_variance') +
  guides(col=guide_legend("model-cv_type")) + ggtitle('Cluster samples')

ggplot(data = clust_sample[clust_sample$sill <= 2000,], aes(x=sill, y=ME)) + 
  geom_line(aes(colour=as.factor(paste0(model,'-',cv_type)))) +xlab('semi_variance') +
  guides(col=guide_legend("model-cv_type")) + ggtitle('Cluster samples')



# random sample  ----------------------------------------------------------

ggplot(data =rand_sample[rand_sample$RMSE <100,], aes(x=as.factor(sample_size), y=RMSE)) + geom_boxplot(aes(fill=paste0(model,'-',cv_type))) +
  ggtitle('Random sample') + xlab('sample size') + guides(fill=guide_legend("model-cv_type"))
ggplot(data =rand_sample[rand_sample$MEC >= 0,], aes(x=as.factor(sample_size), y=MEC)) + geom_boxplot(aes(fill=paste0(model,'-',cv_type))) +
  ggtitle('Random sample') + xlab('sample size') + guides(fill=guide_legend("model-cv_type")) + ylab('propotion explained') 

ggplot(data =rand_sample[rand_sample$RMSE <100,], aes(x=as.factor(sample_size), y=RMSE)) + geom_boxplot(aes(fill=paste0(model,'-',cv_type))) +
  ggtitle('Random sample') + xlab('sample size') + guides(fill=guide_legend("model-cv_type"))
ggplot(data =rand_sample[rand_sample$RMSE <100,], aes(x=as.factor(sample_size), y=RMSE)) + geom_boxplot(aes(fill=paste0(model,'-',cv_type))) +
  ggtitle('Random sample') + xlab('sample size') + guides(fill=guide_legend("model-cv_type"))

ggplot(data =rand_sample[rand_sample$RMSE <100,], aes(x=as.factor(sample_size), y=r2)) + geom_boxplot(aes(fill=paste0(model,'-',cv_type))) +
  ggtitle('Random sample') + xlab('sample size') + guides(fill=guide_legend("model-cv_type"))
ggplot(data =rand_sample[rand_sample$RMSE <100,], aes(x=as.factor(sample_size), y=ME)) + geom_boxplot(aes(fill=paste0(model,'-',cv_type))) +
  ggtitle('Random sample') + xlab('sample size') + guides(fill=guide_legend("model-cv_type"))


#semivariance 
ggplot(data = rand_sample[rand_sample$sill <= 10000,], aes(x= as.factor(sample_size), y=sill)) + 
  geom_boxplot(fill = 'grey') +xlab('sample size') +
  ylab('semi_variance') + ggtitle('Random samples')

# Random sample 
ggplot(data = rand_sample[rand_sample$sill <= 10000,], aes(x=sill, y=RMSE)) + 
  geom_line(aes(colour=as.factor(paste0(model,'-',cv_type)))) +xlab('semi_variance') +
  ylab('RMSE')  + ggtitle('Random Sample') +  guides(col=guide_legend("model-cv_type"))

ggplot(data = rand_sample[rand_sample$sill <= 10000 & rand_sample$MEC >=0,], aes(x=sill, y=MEC)) + 
  geom_line(aes(colour=as.factor(paste0(model,'-',cv_type)))) +xlab('semi_variance') +
  ggtitle('Random Sample') +  guides(col=guide_legend("model-cv_type")) +  ylab('propotion explained')

ggplot(data = rand_sample[rand_sample$sill <= 10000,], aes(x=sill, y=r2)) + 
  geom_line(aes(colour=as.factor(paste0(model,'-',cv_type)))) +xlab('semi_variance') +
  ggtitle('Random Sample') +  guides(col=guide_legend("model-cv_type"))

ggplot(data = rand_sample[rand_sample$sill <= 10000,], aes(x=sill, y=ME)) + 
  geom_line(aes(colour=as.factor(paste0(model,'-',cv_type)))) +xlab('semi_variance') +
  ggtitle('Random Sample') +  guides(col=guide_legend("model-cv_type"))