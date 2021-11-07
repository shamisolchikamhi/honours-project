
## honours-project

This repo conatins my final year honours project in Statistics;
`Assessing the Predictive Performance of StatisticalMethods and Machine Learning Algorithms Using Spatial Cross Validation`.
R was used for this research. 

## Project abstract
The purpose of this project is to `assess the predictive performance of statistical methods and machine learning methods using spatial cross validation`. Two machine learning methods, namely `decision trees` and `random Forrest trees`, are used to model lively woody biomass in an area in the Democratic Republic of the Congo. `Random sampling` and `density based clustering` is used to  sample from the main data set in order to obtain samples of different sizes and spatial distribution. With these samples, the level of `spatial auto-correlation` is measured using the `empirical variogram` method before split this sample into 10-fold for training and validating the model. `10-fold cross validation` is used for conventional cross validation and `10-fold block cross validation` is used for spatial cross validation. `Root mean square error`, `mean absolute error`, `coefficient of determination` and `Pearson correlation coefficient` are used to assess the performance of the machine learning models. This research concludes that `spatial cross validation assess the predictive performance better than the conventional cross validation`. This is because when using spatial data spatial cross validation splits data into `folds that are spatially and statistically independent` unlike the conventional cross validation. Therefore, when spatial data is used to  train machine learning models, spatial cross validation is recommended to assess predictive performance.

## To run 

