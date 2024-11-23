### Cydia pomonella (Codling Moth)
### LVA
### 31 December 2022

################################
# Libraries
require(raster)
require(rgdal)
require(terra)
library(tidyverse)
require(dismo)
require(randomForest)
################################

### Species: Cydia pomonella; Codling Moth

CpCMdata <- read.csv("GBIF_Occurrence/Codling_moth_Cydia_pomonella/occurrence_C_pomonella.csv", header=T)
CpCMpoints <- CpCMdata %>%
  dplyr::select(Long=decimalLongitude, Lat=decimalLatitude)

# Climatic Data

bioLayersCurrent <- rast("bioLayersCurrent.tif") 
bioLayersFuture <- rast("bioLayersFuture.tif") 

### Making names the same for models
library(stringr)
txt <- names(bioLayersCurrent)
newNames<-txt |> 
  str_match("CHELSA_(\\w*)_\\d+.*")
names(bioLayersCurrent) <- newNames[,2]
names(bioLayersCurrent)

txt <- names(bioLayersFuture)
newNames<-txt |> 
  str_match("CHELSA_(\\w*)_\\d+.*")
names(bioLayersFuture) <- newNames[,2]
names(bioLayersFuture)

### Extracted data
# remove layers with lots of NAs
#which(colSums(is.na(bioCurrentData)) > 10)
## remove columns with > 10 NAs
bioLayersCurrent <- terra::subset(bioLayersCurrent, c(1:19, 22:24, 31:39, 41:45))

####### extract data
bioCurrentData <- terra::extract(bioLayersCurrent, CpCMpoints[, c("Long", "Lat")]) 

bioCurrentData$Long <- CpCMpoints$Long
bioCurrentData$Lat <- CpCMpoints$Lat
bioCurrentData <- bioCurrentData %>%
  select(-ID) %>%
  select(Long, Lat, everything())

## Future
## remove columns with > 10 NAs
bioLayersFuture <- terra::subset(bioLayersFuture, c(1:19, 22:24, 31:39, 41:45))

bioFutureData <- terra::extract(bioLayersFuture, CpCMpoints[, c("Long", "Lat")]) 
# remove layers with lots of NAs
#which(colSums(is.na(bioFutureData)) > 11)
bioFutureData$Long <- CpCMpoints$Long
bioFutureData$Lat <- CpCMpoints$Lat
bioFutureData <- bioFutureData %>%
  select(-ID) %>%
  select(Long, Lat, everything())

################################
# Current Models
set.seed(123)
bg <- terra::spatSample(bioLayersCurrent, 3000, "random", na.rm=TRUE, xy=T)
#bg$Presence<- rep(0, dim(bg)[1])
names(bg)[1:2]<-c("Long","Lat")

fullCurrentData<-rbind(bioCurrentData, bg)
fullCurrentData$Presence <- c(rep(1, nrow(CpCMdata)), rep(0,nrow(bg)))

#exploratory
fullCurrentData %>%
  ggplot(aes(Long, Lat), col=Presence) +
  geom_point(size = 0.5, alpha = 0.4) +
  labs(color = NULL)

# ## Training and testing (random subset)
set.seed(123)
group <- kfold(bioCurrentData, 5)
pres_train <- bioCurrentData[group != 1, ]
pres_test <- bioCurrentData[group == 1, ]
group <- kfold(bg, 5)
backg_train <- bg[group != 1, ]
backg_test <- bg[group == 1, ]

train <- rbind(pres_train, backg_train)
pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
envtrain <- raster::extract(bioLayersCurrent, train[,1:2])
envtrain <- data.frame(cbind(Extant=pb_train, envtrain[,2:ncol(envtrain)]) )
#head(envtrain)

## Test
testpres <- data.frame(Extant = c(rep(1, nrow(pres_test))), raster::extract(bioLayersCurrent, pres_test[,1:2])[,2:ncol(envtrain)])
testbackg <- data.frame(Extant = c(rep(0, nrow(backg_test))), raster::extract(bioLayersCurrent, backg_test[,1:2])[,2:ncol(envtrain)])

envtest <- data.frame(rbind(testpres, testbackg))

fullData <- rbind(envtrain, envtest)


library(randomForest)
library(caret)

tr <- trainControl(method = "cv", number = 5)
train(Extant ~ ., data=fullData, method="rf", trControl= tr, na.action=na.exclude) #mtry=2

rf3 <- randomForest(Extant ~ ., data=envtrain, mtry=2, na.action = na.exclude, importance=T)
#Evaluate variable importance
varImpPlot(rf3) # bio15, bio16, bio18
#incnodepurity: gdd10, gsp

erf3 <- dismo::evaluate(p=testpres, a=testbackg, model=rf3)
erf3 #AUC=0.9569; cor: 0.73; max TPR+TNR at : 0.089

# full model
rfFull <- randomForest(Extant ~ ., data=fullData, mtry=2, na.action = na.exclude, importance=T)
varImpPlot(rfFull)
prFull <- predict(bioLayersCurrent, rfFull, type="response")

plot(prFull)

#rfDFCurrent<-terra::as.data.frame(prFull, xy=T)

prFuture <- predict(bioLayersFuture, rfFull, type="response")
plot(prFuture)
#rfDFFuture<-terra::as.data.frame(prFuture, xy=T)

# Rasters
#prFull
#prFuture
#terra::writeRaster(prFull, "CpCM_Current_RF-SDM.tif")
#terra::writeRaster(prFuture, "CpCM_Future_RF-SDM.tif")

# Data Frames
# rfDFCurrent
# rfDFFuture
#write.csv(rfDFCurrent, "CpCM_Current_RF-SDM.csv")
#write.csv(rfDFFuture, "CpCM_Future_RF-SDM.csv")


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


## Training and testing (spatial subset)
library(blockCV)
### info here: https://methodsblog.com/2018/11/29/blockcv-english/
### vignette: https://github.com/rvalavi/blockCV/blob/master/vignettes/tutorial_1.Rmd
### use cv_spatial: https://rdrr.io/github/rvalavi/blockCV/man/cv_spatial.html
# make an sf object from data.frame
#pa_data <- sf::st_as_sf(points, coords = c("x", "y"), crs = 7845)
#sb2 <- cv_spatial(x = pa_data,
# column = "occ",
# rows_cols = c(8, 10),
# k = 5,
# hexagon = FALSE,
# selection = "systematic")

## OR spatialRF https://cran.r-project.org/web/packages/spatialRF/spatialRF.pdf
library(spatialRF)
