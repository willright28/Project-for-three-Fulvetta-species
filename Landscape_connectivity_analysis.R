library(ResistanceGA)
library(raster)
library(parallel)
library(doParallel)
library(JuliaCall)
library(dplyr)
library(geosphere)

#####################################################
#############Load landscape features layers###########
#####################################################


# Land use layers 
land_list=list.files(path="./resistance/land_cover",pattern = "^current.*tif$",full.names = T)
land_raster=stack(land_list)
land_layer=c("current_primf","current_secdf","current_urban")
land_raster_select=subset(land_raster,land_layer)
land_cover=land_raster_select

# Current climate layers
current=raster::stack("./current/chelsa.v2.1.bio.tif")
names(current)=paste0("bio_",1:19)
climate_layer=subset(current,c("bio_4","bio_17","bio_10","bio_11"))

# Topographic layers: elevation and slope
elevation=raster("./resistance/topo/elevation.tif")
slope=raster("./topo/slope.tif")
topo=stack(elevation,slope)

# Resample to unify resolution
topo=resample(topo,climate_layer, method='bilinear')
land_cover=resample(land_cover,climate_layer, method='bilinear') 
range=shapefile("./resistance/range/data_0.shp")

climate_layer=crop(climate_layer,range)
climate_layer=mask(climate_layer,range)

land_cover=crop(land_cover,range)
land_cover=mask(land_cover,range)

topo=crop(topo,range)
topo=mask(topo,range)

# Fst matrix
fst_raw=read.table("./fst.pop.res.txt",header = F)
fst=fst_raw[,3]
fst=fst%>%as.matrix()%>% as.numeric()
matrix=matrix(rep(0,13*13),nrow=13)
matrix[lower.tri(matrix)]=fst
rownames(matrix)=c("2","3","5","8","10","11","14","18","21","24","25","26","27")
colnames(matrix)=rownames(matrix)
fst_dist=matrix



# Sample information
info=read.csv("./info.csv")[,c("long","lat","site_number","pop")]%>%unique()
rownames(info)=info$site_number  
info=info[rownames(fst_dist),]

dist_matrix_m <- distm(info[(rownames(fst_dist)),c("long","lat")], fun = distHaversine)
spdf <- sp::SpatialPoints(coords = info[rownames(fst_dist),c("long","lat")],proj4string = crs(range))
info[rownames(fst_dist),c("long","lat")]

#######################################
####Optimization with CIRCUITSCAPE#####
########################################

JULIA_HOME <- "/home/prunella/HIC_data_flow/phylonet/julia-1.7.0/julia-1.7.0/bin"
JuliaCall::julia_setup(JULIA_HOME)

cl <- makePSOCKcluster(30)
registerDoParallel(cl)

jl.inputs <- jl.prep(n.Pops = length(spdf),
                     response = lower(fst_dist),
                     CS_Point.File = spdf,
                     JULIA_HOME = JULIA_HOME,
                     run_test = F)

GA.inputs <- GA.prep(ASCII.dir = stack(climate_layer,topo,land_cover),
                     Results.dir = "./resistance/result/",
                     min.cat = 1,max.cat = 100,max.cont = 100,
                     method = "AIC",seed = 555,parallel=cl,quiet = TRUE,maxiter=30)

# Export info to cluster
clusterExport(cl=cl,varlist=c("jl.inputs","GA.inputs","climate_layer","topo","land_cover","fst_dist","river_layer")) # list everything you call in ro GA.inputs and gdist
clusterEvalQ(cl=cl, .libPaths("/home/sparrow/.conda/envs/rgdal/lib/R/library")) # set path to where your R library is
clusterCall(cl=cl, library, package = "ResistanceGA", character.only = TRUE)


jl.optim.m <- MS_optim(jl.inputs = jl.inputs,
                    GA.inputs = GA.inputs)

stopCluster(cl)


###########################
######Model Comparison#####
###########################
library(lme4) #v1.1-27.1
library(MuMIn) #1.43.1
library(dplyr)
library(ResistanceGA)
library(raster)
library(parallel)
library(doParallel)


fst_raw=read.table("./fst.pop.res.txt",header = F)
fst=fst_raw[,3]
fst=fst%>%as.matrix()%>% as.numeric()
matrix=matrix(rep(0,13*13),nrow=13)
matrix[lower.tri(matrix)]=fst
rownames(matrix)=c("2","3","5","8","10","11","14","18","21","24","25","26","27")
colnames(matrix)=rownames(matrix)
fst_dist=matrix

#Climate resistance surface
bio4_res <- read.csv("./resistance/result/Results/bio_4_jlResistMat.csv",header = F)
bio4_res <- bio4_res[ lower.tri(bio4_res)]%>% as.data.frame()
bio17_res <- read.csv("./resistance/result/Results/bio_17_jlResistMat.csv",header = F)
bio17_res <- bio17_res[ lower.tri(bio17_res)]%>% as.data.frame()
bio10_res <- read.csv("./resistance/result/Results/bio_10_jlResistMat.csv",header = F)
bio10_res <- bio10_res[ lower.tri(bio10_res)]%>% as.data.frame()
bio11_res <- read.csv("./resistance/result/Results/bio_11_jlResistMat.csv",header = F)
bio11_res <- bio11_res[ lower.tri(bio11_res)]%>% as.data.frame()

#Land use resistance surface
urban_res <- read.csv("./resistance/result/Results/current_urban_jlResistMat.csv",header = F)
urban_res <- urban_res[ lower.tri(urban_res)]%>% as.data.frame()
primf_res <- read.csv("./resistance/result/Results/current_primf_jlResistMat.csv",header = F)
primf_res <- primf_res[ lower.tri(primf_res)]%>% as.data.frame()
secdf_res <- read.csv("./resistance/result/Results/current_secdf_jlResistMat.csv",header = F)
secdf_res <- secdf_res[ lower.tri(secdf_res)]%>% as.data.frame()

#Topology resistance surface
alt_res <- read.csv("./resistance/result/Results/elevation_jlResistMat.csv",header = F)
alt_res <- alt_res[ lower.tri(alt_res)]%>% as.data.frame()
slope_res <- read.csv("./resistance/result/Results/slope_jlResistMat.csv",header = F)
slope_res <- slope_res[ lower.tri(slope_res)]%>% as.data.frame()
slope_res

#Geographic distance 
dist_matrix_m <- distm(info[(rownames(fst_dist)),c("long","lat")], fun = distHaversine)
dis <-  dist_matrix_m 
dis <- dis[ lower.tri(dis)]%>% as.data.frame()
dis

#########

ID <- To.From.ID(13)
Zl <-lapply(c("pop1", "pop2"), function(nm)
      Matrix::fac2sparse(ID[[nm]], "d", drop = FALSE))
ZZ <- Reduce("+", Zl[-1], Zl[[1]])

#########
MLPEnoREML <- function(variables, data) {
  mod2 <- lme4::lFormula(variables, data = data, REML = FALSE)#REML = FALSE
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  mod_2 <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  mod2$reTrms$Zt <- ZZ
  
  # Refit the model
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  modelout <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  return(modelout)
}

df <- cbind(fst,bio4_res,bio17_res,bio10_res,bio11_res,
            urban_res,primf_res,
            secdf_res,alt_res,slope_res,dis)

colnames(df)=c("Genetic_distance","bio4","bio17","bio10","bio11","urban","primf","secdf","alt","slope","dis")
df <- cbind(ID,df)
df <- df %>% mutate(across(c(3:13),scale))#standardization

model_lc_1_aic <- MLPEnoREML(Genetic_distance ~ bio4 + bio17 + bio10 + bio11 + (1|pop1),df) #cliamte
model_lc_2_aic <- MLPEnoREML(Genetic_distance ~ urban + primf + secdf + (1|pop1),df) #land_cover
model_lc_3_aic <- MLPEnoREML(Genetic_distance ~ slope +alt+ (1|pop1),df) #topo
model_lc_4_aic <- MLPEnoREML(Genetic_distance ~ dis + (1|pop2),df) #geographic distance
model_lc_5_aic <- MLPEnoREML(Genetic_distance ~ bio4 + bio17 + bio10 + bio11 + urban + primf + secdf+ alt + slope + dis + (1|pop1),df) #all

aic_1=extractAIC(model_lc_1_aic);aic_1
aic_2=extractAIC(model_lc_2_aic);aic_2
aic_3=extractAIC(model_lc_3_aic);aic_3
aic_4=extractAIC(model_lc_4_aic);aic_4
aic_5=extractAIC(model_lc_5_aic);aic_5

        
##############################################
###Prediction under current conditions###
##############################################

GA.inputs <- GA.prep(ASCII.dir = stack(climate_layer,topo),
                     Results.dir = "./resistance/result/Results/",
                     min.cat = 1,max.cat = 100,max.cont = 100,
                     method = "AIC",seed = 555,parallel=cl,quiet = TRUE,maxiter=30)
jl.optim.m <- MS_optim(jl.inputs = jl.inputs,
                     GA.inputs = GA.inputs)         
jl.out <- Run_CS.jl(jl.inputs = jl.inputs,
                    r = "./resistance/result/Results/elevation.slope.asc",
                    CurrentMap = TRUE,
                    output = 'raster',
                    EXPORT.dir = "./resistance/result/Results/")

###############################################
###Prediction under future climate scenarios ##
###############################################

PARM=c(7,4.086941,93.58588,1,2.031101,56.33411,1,7.116264,46.27791,3,6.664063,78.34327,3,7.209217,99.44295,1,8.109725,23.13104)
IPSL_2070_585=raster::stack("./2070/IPSL/chelsa.v2.1.2070.IPSL.ssp585.bio1_19.resample.tif")
IPSL_2070_585
var=c("bio_4","bio_10","bio_11","bio_17")
names(IPSL_2070_585)=paste0("bio_",1:19)
IPSL_2070_585_pre=subset(IPSL_2070_585,var)
IPSL_2070_585_pre_c=crop(IPSL_2070_585_pre,range)
IPSL_2070_585_pre_m=mask(IPSL_2070_585_pre_c,range)


GA.inputs.IPSL_2070_585 <- GA.prep(ASCII.dir = stack(IPSL_2070_585_pre_m,topo),
                                   Results.dir = "./resistance/result/Results/ipsl_2070_585/",
                                   maxiter = 20,
                                   method = "AIC",
                                   parallel = 4)  


trans.IPSL_2070_585=Combine_Surfaces(PARM=PARM, 
                                     jl.inputs = jl.inputs,
                                     GA.inputs=GA.inputs.IPSL_2070_585, 
                                     out='./resistance/result/Results/ipsl_2070_585/',
                                     File.name='ipsl_2070_585', 
                                     rescale = TRUE, 
                                     p.contribution = T)

#RunCS
jl.out.IPSL_2070_585 <- Run_CS.jl(jl.inputs = jl.inputs,
                                  r = "./resistance/result/Results/ipsl_2070_585/ipsl_2070_585.asc",
                                  CurrentMap = TRUE,
                                  output = 'raster',
                                  EXPORT.dir = "./resistance/result/Results/ipsl_2070_585/")

#MRI
MRI_2070_585=raster::stack("./2070/MRI/chelsa.v2.1.2070.MRI.ssp585.bio1_19.resample.tif")
names(MRI_2070_585)=paste0("bio_",1:19)
MRI_2070_585_pre=subset(MRI_2070_585,var)
MRI_2070_585_pre_c=crop(MRI_2070_585_pre,range)
MRI_2070_585_pre_m=mask(MRI_2070_585_pre_c,range)


GA.inputs.MRI_2070_585 <- GA.prep(ASCII.dir = stack(MRI_2070_585_pre_m,topo),
                                  Results.dir = "./resistance/result/Results/mri_2070_585/",
                                  maxiter = 5,
                                  method = "AIC",
                                  parallel = 4)  


trans.MRI_2070_585=Combine_Surfaces(PARM=PARM,
                                    jl.inputs = jl.inputs,
                                    GA.inputs=GA.inputs.MRI_2070_585, 
                                    out='./resistance/result/Results/mri_2070_585/',
                                    File.name='mri_2070_585', 
                                    rescale = TRUE, 
                                    p.contribution = T)

#RunCS
jl.out.MRI_2070_585 <- Run_CS.jl(jl.inputs = jl.inputs,
                                 r = "./resistance/result/Results/mri_2070_585/mri_2070_585.asc",
                                 CurrentMap = TRUE,
                                 output = 'raster',
                                 EXPORT.dir = "./resistance/result/Results/mri_2070_585/")



######################
###2100
###################
IPSL_2100_585=raster::stack("./2100/IPSL/chelsa.v2.1.2100.IPSL.ssp585.bio1_19.resample.tif")
IPSL_2100_585
var=c("bio_4","bio_10","bio_11","bio_17")
names(IPSL_2100_585)=paste0("bio_",1:19)
IPSL_2100_585_pre=subset(IPSL_2100_585,var)
IPSL_2100_585_pre_c=crop(IPSL_2100_585_pre,range)
IPSL_2100_585_pre_m=mask(IPSL_2100_585_pre_c,range)


GA.inputs.IPSL_2100_585 <- GA.prep(ASCII.dir = stack(IPSL_2100_585_pre_m,topo),
                                   Results.dir = "./resistance/result/Results/ipsl_2100_585/",
                                   maxiter = 20,
                                   method = "AIC",
                                   parallel = 4)  


trans.IPSL_2100_585=Combine_Surfaces(PARM=PARM, 
                                     jl.inputs = jl.inputs,
                                     GA.inputs=GA.inputs.IPSL_2100_585, 
                                     out='./resistance/result/Results/ipsl_2100_585/',
                                     File.name='ipsl_2100_585', 
                                     rescale = TRUE, 
                                     p.contribution = T)

#RunCS
jl.out.IPSL_2100_585 <- Run_CS.jl(jl.inputs = jl.inputs,
                                  r = "./resistance/result/Results/ipsl_2100_585/ipsl_2100_585.asc",
                                  CurrentMap = TRUE,
                                  output = 'raster',
                                  EXPORT.dir = "./resistance/result/Results/ipsl_2100_585/")

#MRI
MRI_2100_585=raster::stack("./2100/MRI/chelsa.v2.1.2100.MRI.ssp585.bio1_19.resample.tif")
names(MRI_2100_585)=paste0("bio_",1:19)
MRI_2100_585_pre=subset(MRI_2100_585,var)
MRI_2100_585_pre_c=crop(MRI_2100_585_pre,range)
MRI_2100_585_pre_m=mask(MRI_2100_585_pre_c,range)


GA.inputs.MRI_2100_585 <- GA.prep(ASCII.dir = stack(MRI_2100_585_pre_m,topo),
                                  Results.dir = "./resistance/result/Results/mri_2100_585/",
                                  maxiter = 5,
                                  method = "AIC",
                                  parallel = 4)  


trans.MRI_2100_585=Combine_Surfaces(PARM=PARM, 
                                    jl.inputs = jl.inputs,
                                    GA.inputs=GA.inputs.MRI_2100_585, 
                                    out='./resistance/result/Results/mri_2100_585/',
                                    File.name='mri_2100_585', 
                                    rescale = TRUE, 
                                    p.contribution = T)

#RunCS
jl.out.MRI_2100_585 <- Run_CS.jl(jl.inputs = jl.inputs,
                                 r = "./resistance/result/Results/mri_2100_585/mri_2100_585.asc",
                                 CurrentMap = TRUE,
                                 output = 'raster',
                                 EXPORT.dir = "./resistance/result/Results/mri_2100_585/")

