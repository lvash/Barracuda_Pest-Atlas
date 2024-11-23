## Lycorma delicatula
## Spotted lanternfly
## LVA
## 2023 25 May

################################
# Libraries
require(raster)
require(rgdal)
require(terra)
library(tidyverse)
require(dismo)
require(randomForest)
library(caret)
################################

### Species: Lycorma delicatula Spotted lanternfly

LdSLdata <- read.csv("GBIF_Occurrence/Lanternfly_Lycorma_delicatula/occurrence_L-delicatula.csv", sep=",", header=T)
LdSLpoints <- LdSLdata %>%
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
#which(colSums(is.na(bioLayersCurrent)) > 10)
## remove columns with > 10 NAs
bioLayersCurrent <- terra::subset(bioLayersCurrent, c(1:19, 22:24, 31:39, 41:45))

####### extract data
bioCurrentData <- terra::extract(bioLayersCurrent, LdSLpoints[, c("Long", "Lat")]) 

bioCurrentData$Long <- LdSLpoints$Long
bioCurrentData$Lat <- LdSLpoints$Lat
bioCurrentData <- bioCurrentData %>%
  select(-ID) %>%
  select(Long, Lat, everything())

## Future
## remove columns with > 10 NAs
bioLayersFuture <- terra::subset(bioLayersFuture, c(1:19, 22:24, 31:39, 41:45))

bioFutureData <- terra::extract(bioLayersFuture, LdSLpoints[, c("Long", "Lat")]) 
# remove layers with lots of NAs
#which(colSums(is.na(bioFutureData)) > 11)
bioFutureData$Long <- LdSLpoints$Long
bioFutureData$Lat <- LdSLpoints$Lat
bioFutureData <- bioFutureData %>%
  select(-ID) %>%
  select(Long, Lat, everything())

################################
# Current Models

bg <- terra::spatSample(bioLayersCurrent, 10000, "random", na.rm=TRUE, xy=T)
#bg$Presence<- rep(0, dim(bg)[1])
names(bg)[1:2]<-c("Long","Lat")

fullCurrentData<-rbind(bioCurrentData, bg)

fullCurrentData$Presence <- c(rep(1, nrow(LdSLdata)), rep(0,nrow(bg)))

#exploratory
fullCurrentData %>%
  ggplot(aes(Long, Lat), col=Presence) +
  geom_point(size = 0.5, alpha = 0.4) +
  labs(color = NULL)


## Training and testing
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

#tr <- trainControl(method = "cv", number = 5)

## warning, code takes a while to run
#params <- caret::train(Extant ~ ., data=fullData, method="rf", trControl= tr, na.action=na.exclude) #mtry=2
#params ##

## all variables
rf3 <- randomForest(Extant ~ ., data=envtrain, na.action = na.exclude, importance=T)
#Evaluate variable importance
varImpPlot(rf3) #%incMSE: bio15, bio9, bio3 #incnodepurity: kg5

erf3 <- dismo::evaluate(p=testpres, a=testbackg, model=rf3)
erf3 #AUC=0.998; cor: 0.976; max TPR+TNR at : 0.473

## important variables
rf4 <- randomForest(Extant ~ bio15 + bio9 + bio3 + kg5, data=envtrain, na.action = na.exclude, importance=T)
varImpPlot(rf4) #
erf4 <- dismo::evaluate(p=testpres, a=testbackg, model=rf4)
erf4 #AUC= 0.9948; cor = 0.96; max TPR+TNR: 0.53


# full model
rfFull <- randomForest(Extant ~ bio15 + bio9 + bio3 + kg5, data=fullData, na.action = na.exclude, importance=T) #, mtry=2
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
terra::writeRaster(prFull, "LdSL_Current_RF-SDM.tif")
terra::writeRaster(prFuture, "LdSL_Future_RF-SDM.tif")
# Data Frames
# rfDFCurrent
# rfDFFuture
#write.csv(rfDFCurrent, "LdSL_Current_RF-SDM.csv")
#write.csv(rfDFFuture, "LdSL_Future_RF-SDM.csv")



###############################
################################### STILL TESTING CODE
################## spatial subsetting using blockCV package
set.seed(623)
spatDat <- fullCurrentData %>%
  select(Long, Lat, Presence)
spatDat <- sf::st_as_sf(spatDat, coords = c("Long", "Lat"), crs = 4326)

scv <- cv_cluster(x = spatDat,
                  column = "Presence", 
                  k = 11)

range <- cv_spatial_autocor(
  x = spatDat, # species data
  column = "Presence", # column storing presence-absence records (0s and 1s)
  plot = F
)

range$range # 1113506

blockCV::cv_plot(cv = scv, 
                 x = spatDat)

################################### 

scv1 <- cv_spatial(
  x = spatDat,
  column = "Presence", # the response column (binary or multi-class)
  r = bioLayersCurrent,
  k = 9, # number of folds
  size = 1113506, # size of the blocks in metres
  selection = "random", # random blocks-to-fold
  iteration = 50, # find evenly dispersed folds
  progress = FALSE, # turn off progress bar
  biomod2 = F, # also create folds for biomod2
  raster_colors = terrain.colors(10, rev = TRUE) # options from cv_plot for a better colour contrast
) 

scv2 <- cv_nndm(
  x = spatDat,
  column = "Presence",
  r = bioLayersCurrent,
  #size = 1113506, # range of spatial autocorrelation
  #num_sample = 1000, # number of samples of prediction points
  sampling = "regular", # sampling methods; it can be random as well
  min_train = 0.1, # minimum portion to keep in each train fold
  plot = TRUE
)

folds <- scv2$folds_list

# create a data.frame to store the prediction of each fold (record)
test_table <- pa_data
test_table$preds <- NA

for(k in seq_len(length(folds))){
  # extracting the training and testing indices
  # this way works with folds_list list (but not folds_ids)
  trainSet <- unlist(folds[[k]][1]) # training set indices; first element
  testSet <- unlist(folds[[k]][2]) # testing set indices; second element
  rf <- randomForest(occ ~ ., model_data[trainSet, ], ntree = 500) # model fitting on training set
  test_table$preds[testSet] <- predict(rf, model_data[testSet, ], type = "prob")[,2] # predict the test set
}

# calculate Area Under the ROC and PR curves and plot the result
precrec_obj <- evalmod(scores = test_table$preds, labels = test_table$occ)
auc(precrec_obj)

#############################################

