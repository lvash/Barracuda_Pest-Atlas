### Rasters for Ag Booth
### LVA
### 29 December 2022

################################
# Libraries
#require(rgdal)
require(terra)
library(tidyverse)
require(randomForest)
################################

################################
# Climatic Data

bioLayersCurrent <- rast("bioLayersCurrent.tif") 
bioLayersFuture <- rast("bioLayersFuture.tif") 

################################
# Point Data

### Species: Drosophila suzukii; spotted wing drosophila (SWD)
SWDdata <- read.csv("GBIF_Occurrence/SWD_Drosophila_suzukii/occurrence_D_suzukii.csv", header=T)
SWDpoints <- SWDdata %>%
  dplyr::select(Long=decimalLongitude, Lat=decimalLatitude)

### Extracted data
## remove layers with > 2% NAs
bioData <- terra::extract(bioLayersCurrent, SWDpoints[, c("Long", "Lat")]) 
badVars <- names(which(colSums(is.na(bioCurrentData)) > 0.02*nrow(SWDpoints))) 
bioLayersCurrent <- terra::subset(bioLayersCurrent, which(!names(bioLayersCurrent) %in% badVars))

####### extract data
bioCurrentData <- terra::extract(bioLayersCurrent, SWDpoints[, c("Long", "Lat")]) 

bioCurrentData$Long <- SWDpoints$Long
bioCurrentData$Lat <- SWDpoints$Lat
bioCurrentData <- bioCurrentData %>%
  dplyr::select(-ID) %>%
  dplyr::select(Long, Lat, everything())

## Future
## remove columns with > 5% NAs
bioLayersFuture <- terra::subset(bioLayersFuture, which(!names(bioLayersFuture) %in% badVars))

bioFutureData <- terra::extract(bioLayersFuture, SWDpoints[, c("Long", "Lat")]) 

#which(colSums(is.na(bioFutureData)) > 11)
bioFutureData$Long <- SWDpoints$Long
bioFutureData$Lat <- SWDpoints$Lat
bioFutureData <- bioFutureData %>%
  dplyr::select(-ID) %>%
  dplyr::select(Long, Lat, everything())


################################
# Current Model
######
set.seed(123)
backg <- terra::spatSample(bioLayersCurrent, 10000, "random", na.rm=TRUE, xy=T)
names(backg)[1:2]<-c("Long","Lat")

#### exploratory
# fullCurrentData<-rbind(bioCurrentData, backg)
# occ <- factor(c(rep("1", nrow(bioCurrentData)), rep("0", nrow(backg))))
# fullCurrentData <- cbind(occ, fullCurrentData)
# fullCurrentData %>%
#   ggplot(aes(Long, Lat)) +
#   geom_point(size = 0.5, alpha = 0.4, aes(colour = occ)) +
#   labs(color = NULL)


## Training and testing
set.seed(123)
group <- kfold(bioCurrentData, 5)
pres_train <- bioCurrentData[group != 1, ]
pres_test <- bioCurrentData[group == 1, ]
group <- kfold(backg, 5)
backg_train <- backg[group != 1, ]
backg_test <- backg[group == 1, ]

train <- rbind(pres_train, backg_train)
pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
envtrain <- raster::extract(bioLayersCurrent, train[,1:2])
envtrain <- data.frame(cbind(Extant=pb_train, envtrain[,2:ncol(envtrain)]) )

## Test
testpres <- data.frame(Extant = c(rep(1, nrow(pres_test))), raster::extract(bioLayersCurrent, pres_test[,1:2])[,2:ncol(envtrain)])
testbackg <- data.frame(Extant = c(rep(0, nrow(backg_test))), raster::extract(bioLayersCurrent, backg_test[,1:2])[,2:ncol(envtrain)])
envtest <- data.frame(rbind(testpres, testbackg))

fullData <- rbind(envtrain, envtest)

# Random Forest (regression, not classification)
library(randomForest)
library(caret)

tr <- trainControl(method = "cv", number = 5)
rf_train <- train(Extant ~ ., data=fullData, method="rf", trControl= tr, na.action=na.exclude) #mtry=2
# options here: http://topepo.github.io/caret/train-models-by-tag.html#Boosting

rf <- randomForest(Extant ~ ., data=envtrain, mtry=rf_train$bestTune$mtry, na.action = na.exclude, importance=T)

# Evaluate variable importance (use caution - does not indicate biological significance)
varImpPlot(rf) 

# bio19, bio18, bio15, bio13, bio4
# bio19=precip of coldest (18=warmest) quarter
# bio15 precipitation seasonality
# bio13 precipitation of wettest month
# bio4 temp seasonality

erf <- dismo::evaluate(p=testpres, a=testbackg, model=rf)
erf

# full model using all data (alternatively you could use the training data for the final model)
rfFull <- randomForest(Extant ~ ., data=fullData, mtry=2, na.action = na.omit, importance=T)
# varImpPlot(rfFull) #bio18,15,13 #nodepurity: bio15,npp
prFull <- predict(bioLayersCurrent, rfFull, type="response")

plot(prFull)

rfDFCurrent<-terra::as.data.frame(prFull, xy=T)

prFuture <- predict(bioLayersFuture, rfFull, type="response")
plot(prFuture)
rfDFFuture<-terra::as.data.frame(prFuture, xy=T)

# Save Rasters
#prFull
#prFuture
terra::writeRaster(prFull, "Pest Atlas Pages/Maps/SWD_Current_RF-SDM.tif", overwrite=TRUE)
terra::writeRaster(prFuture, "Pest Atlas Pages/Maps/SWD_Future_RF-SDM.tif", overwrite=TRUE)

# Data Frames
# rfDFCurrent
# rfDFFuture
#write.csv(rfDFCurrent, "SWD_Current_RF-SDM.csv")
#write.csv(rfDFFuture, "SWD_Future_RF-SDM.csv")


###############################
################################### STILL TESTING CODE
################## spatial subsetting using blockCV package
# set.seed(623)
# spatDat <- fullCurrentData %>%
#   select(Long, Lat, Presence)
# spatDat <- sf::st_as_sf(spatDat, coords = c("Long", "Lat"), crs = 4326)
# 
# scv <- cv_cluster(x = spatDat,
#                   column = "Presence", 
#                   k = 11)
# 
# range <- cv_spatial_autocor(
#   x = spatDat, # species data
#   column = "Presence", # column storing presence-absence records (0s and 1s)
#   plot = F
# )
# 
# range$range # 1113506
# 
# blockCV::cv_plot(cv = scv, 
#                  x = spatDat)
# 
# ################################### 
# 
# scv1 <- cv_spatial(
#   x = spatDat,
#   column = "Presence", # the response column (binary or multi-class)
#   r = bioLayersCurrent,
#   k = 9, # number of folds
#   size = 1113506, # size of the blocks in metres
#   selection = "random", # random blocks-to-fold
#   iteration = 50, # find evenly dispersed folds
#   progress = FALSE, # turn off progress bar
#   biomod2 = F, # also create folds for biomod2
#   raster_colors = terrain.colors(10, rev = TRUE) # options from cv_plot for a better colour contrast
# ) 
# 
# scv2 <- cv_nndm(
#   x = spatDat,
#   column = "Presence",
#   r = bioLayersCurrent,
#   #size = 1113506, # range of spatial autocorrelation
#   #num_sample = 1000, # number of samples of prediction points
#   sampling = "regular", # sampling methods; it can be random as well
#   min_train = 0.1, # minimum portion to keep in each train fold
#   plot = TRUE
# )
# 
# folds <- scv2$folds_list
# 
# # create a data.frame to store the prediction of each fold (record)
# test_table <- pa_data
# test_table$preds <- NA
# 
# for(k in seq_len(length(folds))){
#   # extracting the training and testing indices
#   # this way works with folds_list list (but not folds_ids)
#   trainSet <- unlist(folds[[k]][1]) # training set indices; first element
#   testSet <- unlist(folds[[k]][2]) # testing set indices; second element
#   rf <- randomForest(occ ~ ., model_data[trainSet, ], ntree = 500) # model fitting on training set
#   test_table$preds[testSet] <- predict(rf, model_data[testSet, ], type = "prob")[,2] # predict the test set
# }
# 
# # calculate Area Under the ROC and PR curves and plot the result
# precrec_obj <- evalmod(scores = test_table$preds, labels = test_table$occ)
# auc(precrec_obj)

#############################################


