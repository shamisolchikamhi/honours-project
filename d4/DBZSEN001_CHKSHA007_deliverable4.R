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




# importing tif files -----------------------------------------------------


root_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_dir)

# make a raster stack 
files <- list.files(path = paste0(root_dir, '/data2/'), 
                    recursive = FALSE, pattern = "\\.tif$")
s <- stack(paste0(root_dir, '/data2/', files))
df <- as.data.frame(s, xy=T, na.rm=T)
head(df)


test <- df
gridded(test) = ~x+y
spplot(test[1])


# Random sampling-------------------------------------------------------

s_sizes <- c(500, 1000,1500)

for (samp in s_sizes){
pt.SRS_500 <- df[sample(nrow(df), samp, replace = F), ]
valuetable <- na.omit(as.data.frame(pt.SRS_500))
coordinates(valuetable) <- ~x+y
plot(s[[1]])
plot(valuetable, add=T)

set.seed(1998)
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

