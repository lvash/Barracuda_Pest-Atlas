### Processing Bioclimatic Layers
### LVA
### created 29 December 2022
### edited 27 August 2024

################################
# Libraries
require(terra)
library(stringr)
library(geodata)

################################

### Files were originally downloaded from the CHELSA Bioclim + website: https://chelsa-climate.org/exchelsa-extended-bioclim/ Download > climatologies > 1981-2010 > bio: 

### All files were selected under 2041-2070 SSP 585 (only GFDL-ESM4 earth model) bio folder and the corresponding files were downloaded for 1981-2010 bio

### Advice for future: average all 5 earth models instead of choosing 1. Average predictions if you want to include uncertainty (std dev, etc) or rasters themselves if not interested in reporting uncertainty

################################

extent <- ext(-100, -60, 25, 50)

##### Current CHELSA Bioclim+ 
files <- list.files(path="CHELSA/climatologies/ssp585/1981-2010/", pattern='\\.tif$', full.names = TRUE)
bioLayersCurrent <- rast(files)
bioLayersCurrent <- crop(bioLayersCurrent, extent)

# geodata's world function to download our base map
world_map <- world(resolution = 3,
                   path = ".")
w <- crop(world_map, ext(-100, -60, 25, 50))
bioLayersCurrent <- terra::mask(bioLayersCurrent, w) #mask to exclude ocean values

## Renaming to base names 
txt <- names(bioLayersCurrent)
newNames<- txt |> 
  str_match("CHELSA_(\\w*)_\\d+.*")
names(bioLayersCurrent) <- newNames[,2]
names(bioLayersCurrent)

writeRaster(bioLayersCurrent, "bioLayersCurrent.tif", overwrite=TRUE)

##### Future
files2 <- list.files("CHELSA/climatologies/ssp585/2041-2070/bio/", pattern='\\.tif$', full.names = TRUE)
bioLayersFuture <- rast(files2)
bioLayersFuture <- crop(bioLayersFuture, extent)
bioLayersFuture <- terra::mask(bioLayersFuture, w) # mask to exclude ocean values

## Renaming to base names, so predictors between raster stacks match 
txt <- names(bioLayersFuture)
newNames<- txt |> 
  str_match("CHELSA_(\\w*)_\\d+.*")
names(bioLayersFuture) <- newNames[,2]
names(bioLayersFuture)
terra::writeRaster(bioLayersFuture, "bioLayersFuture.tif", overwrite = TRUE)