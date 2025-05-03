library(terra)       
library(biomod2)      
library(raster)      
library(ggplot2)      
library(rgdal)       

##### Using Alcippe hueti as an example; the same code can be applied to the other two species

crdref <- CRS('+proj=longlat +datum=WGS84')
# Selected bioclimatic variables
pre.var = c("bio_10", "bio_4", "bio_11", "bio_17")
# Load current climate layers
current = stack("./chelsa.v2.1.bio.tif")
names(current) = paste0("bio_", 1:19)
current_sub = raster::subset(current, pre.var)
# 2070 climate layers
f2070_ipsl = stack("./2070/IPSL/chelsa.v2.1.2070.IPSL.ssp585.bio1_19.resample.tif")
names(f2070_ipsl) = paste0("bio_", 1:19)
f2070_ipsl = raster::subset(f2070_ipsl, pre.var)

f2070_mri = stack("./2070/MRI/chelsa.v2.1.2070.MRI.ssp585.bio1_19.resample.tif")
names(f2070_mri) = paste0("bio_", 1:19)
f2070_mri = raster::subset(f2070_mri, pre.var)

f2070_585 = (f2070_ipsl + f2070_mri) / 2

# 2100 climate layers
f2100_ipsl = stack("./2100/IPSL/chelsa.v2.1.2100.IPSL.ssp585.bio1_19.resample.tif")
names(f2100_ipsl) = paste0("bio_", 1:19)
f2100_ipsl = raster::subset(f2100_ipsl, pre.var)

f2100_mri = stack("./2100/MRI/chelsa.v2.1.2100.MRI.ssp585.bio1_19.resample.tif")
names(f2100_mri) = paste0("bio_", 1:19)
f2100_mri = raster::subset(f2100_mri, pre.var)

f2100_585 = (f2100_ipsl + f2100_mri) / 2

### Load filtered occurrence points
acsel = read.csv("./filter.10km.hueti.csv", header = TRUE)[, c("x", "y")]

# Define projection extent
extent = extent(90, 123, 8, 38)
current_crop = crop(current_sub, extent)

coast_i = shapefile("./coast/coast_pa_range/coast_i.shp")
current_crop_coast = raster::mask(current_crop, coast_i)

######## Run biomod2
acsel$presence = 1
pa_df = acsel
myRespName <- "alcippe"
myResp <- as.numeric(pa_df$presence)
myRespXY <- data.frame(pa_df[,1:2])
myExpl = stack(current_crop_coast)
names(myExpl)

myBiomodData.r <- BIOMOD_FormatingData(resp.var = myResp, # presence points
                                       expl.var = stack(myExpl), # environmental layers
                                       resp.xy = myRespXY, # coordinates of presence points
                                       resp.name = myRespName, # species name
                                       PA.nb.rep = 3,
                                       PA.nb.absences = 10000,
                                       PA.strategy = 'random',
                                       seed.val = 345)

plot(myBiomodData.r)

###### Initial modeling to select appropriate algorithms
allModels <- c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM',
               'MAXENT', 'MAXNET', 'RF', 'RFd', 'SRE', 'XGBOOST')

# Default model options
myBiomodOptions <- bm_ModelingOptions(data.type = 'binary',
                                      models = allModels,
                                      strategy = 'default')

myBiomodModelOut_test <- BIOMOD_Modeling(bm.format = myBiomodData.r,
                                         modeling.id = "ALLmodels",
                                         models = allModels,
                                         nb.rep = 3,
                                         data.split.perc = 80,
                                         metric.eval = c('TSS','ROC'), # evaluation metrics
                                         nb.cpu = 10)

bm_PlotEvalMean(bm.out = myBiomodModelOut_test) # Evaluate models

########### Fine-tuning selected models

tuned.glm <- bm_Tuning(model = 'GLM',
                       tuning.fun = 'glm',
                       do.formula = FALSE,
                       bm.options = myBiomodOptions@options$GLM.binary.stats.glm,
                       bm.format = myBiomodData.r)

tuned.gbm <- bm_Tuning(model = 'GBM',
                       tuning.fun = 'gbm',
                       do.formula = FALSE,
                       bm.options = myBiomodOptions@options$GBM.binary.gbm.gbm,
                       bm.format = myBiomodData.r)

tuned.maxent <- bm_Tuning(model = 'MAXENT',
                          tuning.fun = 'ENMevaluate',
                          do.formula = FALSE,
                          metric.eval = "AICc",
                          bm.options = myBiomodOptions@options$MAXENT.binary.MAXENT.MAXENT,
                          bm.format = myBiomodData.r,
                          params.train = list(MAXENT.algorithm = 'maxnet', MAXENT.parallel = FALSE))

tuned.xgboost <- bm_Tuning(model = 'XGBOOST',
                           tuning.fun = 'xgbTree',
                           do.formula = FALSE,
                           bm.options = myBiomodOptions@options$XGBOOST.binary.xgboost.xgboost,
                           bm.format = myBiomodData.r)

new.options <- bm_ModelingOptions(
  GLM = tuned.glm$best.params,
  GBM = tuned.gbm$best.params,
  MAXENT = tuned.maxent$best.params,
  XGBOOST = tuned.xgboost$best.params
)

###
n = 4623
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData.r,
                                    modeling.id = 'hueti',
                                    models = c("MAXENT", "GBM", "XGBOOST", "GLM"),
                                    nb.rep = 5,
                                    models.options = new.options,
                                    data.split.perc = 70,
                                    metric.eval = c('TSS','ROC'),
                                    do.full.models = FALSE,
                                    nb.cpu = 8,
                                    seed.val = n)

bm_PlotEvalMean(bm.out = myBiomodModelOut)


# Ensemble Modeling
myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      metric.select = c('ROC','TSS'),
                                      metric.select.thresh = c(0.8, 0.6),
                                      metric.eval = c('TSS', 'ROC'),
                                      em.algo = c('EMwmean'),
                                      EMwmean.decay = 'proportional',
                                      nb.cpu = 12,
                                      seed.val = n)

##### Predict current distribution
myBiomodEMProj_current <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                                     proj.name = 'central_current',
                                                     new.env = current_crop_coast,
                                                     models.chosen = 'all',
                                                     metric.binary = c('TSS'),
                                                     metric.filter = c('TSS'),
                                                     seed.val = n)

#### Predict future distribution

########### 2070 under SSP585
f2070_crop = crop(f2070_585, extent)
f2070_crop_coast = raster::mask(f2070_crop, coast_i)

Proj_2070 <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                        proj.name = 'hueti_2070_585',
                                        build.clamping.mask = TRUE,
                                        new.env = f2070_crop_coast,
                                        models.chosen = "all",            
                                        metric.binary = c('TSS'),
                                        metric.filter = c('TSS'))

############# 2100 under SSP585
f2100_crop = crop(f2100_585, extent)
f2100_crop_coast = raster::mask(f2100_crop, coast_i)

Proj_2100 <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                        proj.name = 'hueti_2100_585',
                                        new.env = f2100_crop_coast,
                                        models.chosen = 'all',  
                                        metric.binary = 'TSS',
                                        metric.filter = 'TSS')





