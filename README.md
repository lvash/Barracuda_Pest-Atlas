# Barracuda Pest Atlas

We predicted the range shift of 5 agricultural pest species: 
* Codling Moth (*Cydia pomonella*)  
* Brown marmorated stinkbug (*Halyomorpha halys*)  
* Spotted lantern fly (*Lycorma delicatula*)  
* Cabbage White Butterfly (*Pieris rapae*)  
* Tufted Apple Bud Moth (*Platynota idaeusalis*)  
* Spotted wing *Drosophila* (*D. suzukii*)  

### Files and folders within this repository:

The **ODMap csv file** contains specifics about the project, predictors, and modeling process and follows the format of [ODMAP (Overview, Data, Model, Assessment and Prediction) protocol](https://odmap.wsl.ch/). 

Occurrence data for each species were obtained via iNaturalist GBIF Research-grade observations and are located in the **GBIF_Occurrence folder.**
Climatic predictors for current (1981-2010) and future (2041-2070) decades were obtained from [CHELSA BioClim + ](https://chelsa-climate.org/exchelsa-extended-bioclim/ and cropped and processed in the Raster-Processing.R script within the **Model Scripts folder.** Each species has its own script within that folder which uses occurrence data to train/test/predict and plots the current and future distributions using a Random Forest model. 

Predicted Distribution Maps were saved within the **Pest Atlas Pages folder** and used to create RMarkdown and html pages for each species within the same folder. The pages were converted to Rpubs and stored in the **rsconnect folder** within the Pest Atlas Pages folder.

The home page of the Pest Atlas can be found [here](https://rpubs.com/lvash/Pest_Atlas_Home).
