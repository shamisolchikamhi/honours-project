
### honours-project

This repository conatins our final year honours project in Statistics;
`Assessing the Predictive Performance of StatisticalMethods and Machine Learning Algorithms Using Spatial Cross Validation`.
R was used for this research. 

### Project abstract
The purpose of this project is to `assess the predictive performance of statistical methods and machine learning methods using spatial cross validation`. Two machine learning methods, namely `decision trees` and `random Forrest trees`, are used to model lively woody biomass in an area in the Democratic Republic of the Congo. `Random sampling` and `density based clustering` is used to  sample from the main data set in order to obtain samples of different sizes and spatial distribution. With these samples, the level of `spatial auto-correlation` is measured using the `empirical variogram` method before split this sample into 10-fold for training and validating the model. `10-fold cross validation` is used for conventional cross validation and `10-fold block cross validation` is used for spatial cross validation. `Root mean square error`, `mean absolute error`, `coefficient of determination` and `Pearson correlation coefficient` are used to assess the performance of the machine learning models. This research concludes that `spatial cross validation assess the predictive performance better than the conventional cross validation`. This is because when using spatial data spatial cross validation splits data into `folds that are spatially and statistically independent` unlike the conventional cross validation. Therefore, when spatial data is used to  train machine learning models, spatial cross validation is recommended to assess predictive performance.

### To run 
#### Download the code V1 file 
##### Data.R (optional)

1. Download the data files `feb.zip` and `data.zip` from https://drive.google.com/drive/folders/1lPJkx10lKGL0TtcyfjGTQzL1jIcTKaFV?usp=sharing
   The two folders contain data for the covariates. They are saved differently to enable stacking as they have different extents. 
2. These folders should be unzipped and saved in their respective folders in the `code V1` folder
3. The script can be run normmarly to save the raster layers to  `processed_data` that is in the same extent as the biomass data and cropped to the required Congo region

##### Script.R

1. If `data.R` was run then continue to run `script.R`, else
2. Download `processed_data.zip` from https://drive.google.com/drive/folders/1lPJkx10lKGL0TtcyfjGTQzL1jIcTKaFV?usp=sharing that contains the procesed data. 
3. The script can be run normally collect the resulst in csv files that are later used for plots
