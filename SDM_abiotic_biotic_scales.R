# This codes was readapted from https://cran.r-project.org/web/packages/biomod2/vignettes/examples_1_mainFunctions.html 

#Read all the data 
setwd("directory")

# Read the shp where you have the presence and absence points
pres_abs <- terra::vect("PresAbsSericea.shp")
names(pres_abs)
pres_abs$Latitude
pa <- data.frame(pres_abs$Northing, pres_abs$Easting, pres_abs$PA)
names(pa) <- c("pres_abs.Northing", "pres_abs.Easting", "pres_abs.PA")
pres <- dplyr::filter(pa, pres_abs.PA=="1")
abs <- dplyr::filter(pa, pres_abs.PA=="0")

library(terra)
# Read the variables to be used as predictors in SDMs
#####read DEM
dtm_1 <- rast("DEM_north.tif")
plot(dtm)

# if aggregate is needed apply the following line
dtm <- aggregate(dtm_1, fact=3, fun=mean) 

#Calculate slope based on dtm or dtm_1 depending on the spatial resolution 
slope <- terrain(dtm, v = 'slope', unit = 'degrees', neighbors = 8)

plot(slope)

# Read annual mean land surface temperature 
lst_mean <- rast("mean_LST_3m_resampled.tif") # all months

# Read Time since fire
TSF_one <- rast('Clip_TSM_ma.tif')

# aggregate time since fire if needed
TSF <- aggregate(TSF_one,3)
plot(TSF)

#  Read the raster of the functional traits

TN <- rast("Clip_mean_50N.tif")
P <- rast("Clip_mean_50P.tif")
K <- rast("Clip_mean_50K.tif")
Height <- rast("Clip_mean_50Height.tif")

######## make them same extent if not the same 
lst_mean <- resample(lst_mean, slope)
TSF <- resample(TSF, slope)
TN <- resample(TN, slope)
P <- resample(P, slope)
K <- resample(K, slope)
Height <- resample(Height, slope)


# check also the crs and make sure they are the same
# crs(lst_mean) <- crs(slope)
# crs(TSF) <- crs(slope)
# crs(TN) <- crs(slope)
# crs(P) <- crs(slope)
# crs(K) <- crs(slope)
# crs(Height) <- crs(slope)

###############################NOW PUT TOGETHER THE PREDICTORS THAT WE WANT TO USE FOR SDMs 
######################################RUN THE FOLLOWING LINES BASED ON NEEDED
## Together abiotic alone 
bioclim_world <- c(slope,lst_mean,  TSF) # abiotic alone 
names(bioclim_world) <- c("slope","lst_mean", "TSF")
plot(bioclim_world)


## Together biotic and abiotic  
bioclim_world <- c(slope,lst_mean, TSF, TN, P, K, Height)
names(bioclim_world) <- c("slope","lst_mean", "TSF", "TN", "P", "K", "Height")
plot(bioclim_world)


## Biotic alone
bioclim_world <- c(TN, P, K, Height)
names(bioclim_world) <- c("TN", "P", "K", "Height")
plot(bioclim_world)

getwd()
writeRaster(bioclim_world, "Env_variable.tif", overwrite=T)

# format data 
# Build and run a range of models using biomod2
# the species occurrences and associated coordinates
# the environmental conditions
# the name of the species of interest
# build pseudo-absence data
## Biomod2 Formatting
ProLau_data <- biomod2::BIOMOD_FormatingData(resp.var = pa$pres_abs.PA #presence only
                                             ,expl.var = bioclim_world       #environmental variables
                                             , resp.xy =  pa[, c('pres_abs.Easting'              #long
                                                                 , 'pres_abs.Northing')]              #latitude
                                             #species name
                                             , PA.nb.rep = 2  #numbers of iterations cross validation 
                                             ,resp.name="sericea"    #to set up the pseudo absence use random
)
terra::plot(ProLau_data)


# tune model parameter # if only running RF then just keep RF in this tuning
ProLau_opt <- biomod2::BIOMOD_ModelingOptions(GLM = list(type = 'quadratic'
                                                         , interaction.level = 1
)
, GBM = list(n.trees = 1000)
, GAM = list( algo = 'GAM_mgcv')
, RF= list(ntree= 500, mtry= 4)
)


# run model
# https://cran.r-project.org/web/packages/biomod2/vignettes/examples_1_mainFunctions.html

# Build models with 50 iterations
ProLau_models <- biomod2::BIOMOD_Modeling(
  bm.format = ProLau_data                  #Input data to be used(from the formatting data)
  , models = 'RF' #c('GLM', 'GBM', 'RF', 'GAM') can also add this all
  , bm.options = ProLau_opt                #model parameters that you just set earlier  
  , nb.rep = 50 # number of repetitions to be done #split the data
  # for calibration/validation splitting 
  , data.split.perc = 80 # 80% of the data for training    
  , var.import = 3 # number of permutations to be done for each
  # variable to estimate variable importance
  , metric.eval = c('TSS','ROC')
  , do.full.models = FALSE
)


# save models
save(ProLau_models, file = 'ProLau_models.Rdata')

# model accuracy
####################################
# load model
load(file = 'ProLau_models.Rdata')

# Decomposing the models?? variability / accuracy 
## get model evaluation scores
## a 5 dimension array containing the scores for the models
ProLau_models_scores <- biomod2::get_evaluations(ProLau_models)
dimnames(ProLau_models_scores)
getwd()

# save the scores
write.csv(ProLau_models_scores, "RF_get_evaluations.csv")

## assess the influence of the different choices made when parameterizing the models
### algorithm
### RF models seem to be the most accurate
# biomod2::models_scores_graph(ProLau_models
#                             , by = 'models'
#                             , metrics = c('ROC', 'TSS')
#                             , xlim = c(0.5,1)
#                             , ylim = c(0.5,1)
# )
biomod2::bm_PlotEvalBoxplot(bm.out = ProLau_models
                            , group.by = c('algo', 'algo')
)
### cross-validation run
# models_scores_graph(ProLau_models
#                     , by = 'cv_run'
#                     , metrics = c('ROC', 'TSS')
#                     , xlim = c(0.5,1)
#                     , ylim = c(0.5,1)
#                    )
biomod2::bm_PlotEvalBoxplot(bm.out = ProLau_models        #ACCROSS DIFFERENT WAY TO RUN THE DATA
                            , group.by = c('run', 'run')
)
### pseudo-absences sampling
# models_scores_graph(ProLau_models
#                     , by = 'data_set'
#                     , metrics = c('ROC', 'TSS')
#                     , xlim = c(0.5, 1)
#                     , ylim = c(0.5, 1)
#                    )
biomod2::bm_PlotEvalBoxplot(bm.out = ProLau_models       # VARIATION IN ACCURACY ACCROSS PA CHOOSEN
                            , group.by = c('PA', 'PA')
)
####################################


# variable importance + response curve
####################################
# load model
load(file = 'ProLau_models.Rdata')

# mean of variable importance by algorithm
## get vi  # Relative importance of one variables, higher value more important. To know which variables driving patterns of variation
ProLau_models_var_import <- biomod2::get_variables_importance(ProLau_models)
## summarize
var_im <- tapply(ProLau_models_var_import$var.imp
                 , list(ProLau_models_var_import$expl.var
                        , ProLau_models_var_import$algo
                 )
                 , mean)

write.csv(ProLau_models_var_import, "variable_importance_FT.csv")
## plot
biomod2::bm_PlotVarImpBoxplot(bm.out = ProLau_models
                              , group.by = c('expl.var', 'algo', 'algo'))

# response curve of each variable
## RF
biomod2::bm_PlotResponseCurves(
  bm.out = ProLau_models
  , models.chosen = biomod2::get_built_models(ProLau_models, algo = 'RF')
  , fixed.var = 'median'
)

#################################### projections, run the models across the whole study area for all the 50 models
getwd()
setwd("outputdirectory")
load(file = 'ProLau_models.Rdata')
bioclim_world <- rast("Env_variable.tif")

for (i in 1:9){
  mods <- biomod2::get_built_models(ProLau_models
                                    , run = paste('RUN', i, sep="")
                                    , algo = 'RF' 
                                    #  "GBM"
                                    , PA = 'PA1')
  
  ## run projection
  ProLau_models_proj_current_1mdl <- biomod2::BIOMOD_Projection(
    bm.mod = ProLau_models
    , proj.name = paste('current',i,sep="")
    , new.env = bioclim_world
    , models.chosen = mods
    , metric.binary = 'all'# evaluation metric
    # selected to transform
    # predictions into
    # binary values
    , metric.filter = 'all'
    , build.clamping.mask = TRUE
    #, output.format = '.img'
    #, do.stack = FALSE
  )
}



##########SOME PLOTS

library(RColorBrewer)
airborne <- raster("D:/Applications/AGU/2023/airborne_but_bandsuntil927/sericea/proj_currentRF1/proj_current_sericea.tif")
mypal <- brewer.pal(n = 20, name = "Accent")
plot(airborne)

getwd()
# set to folder where we have the projection folers
setwd("outputdirectory/sericea") 

library(raster)

pred_proj1 <- stack("./proj_current1/proj_current1_sericea.tif")
pred_proj2 <- stack("./proj_current2/proj_current2_sericea.tif")
pred_proj3 <- stack("./proj_current3/proj_current3_sericea.tif")
pred_proj4 <- stack("./proj_current4/proj_current4_sericea.tif")
pred_proj5 <- stack("./proj_current5/proj_current5_sericea.tif")
pred_proj6 <- stack("./proj_current6/proj_current6_sericea.tif")
pred_proj7 <- stack("./proj_current7/proj_current7_sericea.tif")
pred_proj8 <- stack("./proj_current8/proj_current8_sericea.tif")
pred_proj9 <- stack("./proj_current9/proj_current9_sericea.tif")

# stack the 50 runs
pred_proj <- stack(pred_proj1, pred_proj2, pred_proj3, pred_proj4, pred_proj5, pred_proj6,pred_proj7, pred_proj8, pred_proj9)

# calculate mean of the 50 runs
calc_pred <- calc(pred_proj, mean)

plot(calc_pred_val)
calc_pred_val <- calc_pred*100/1000

getwd()
# write raster final
writeRaster(calc_pred_val, "predicted_sericea_mean_50run_perc.tif")

