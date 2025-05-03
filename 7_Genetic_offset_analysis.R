library(dplyr)
library(raster)
library(LEA)
require(gradientForest)
library(ggplot2)
library(ggpubr)
library(RColorBrewer) 

# This script calculates genetic offset using two methods: Gradient ForestGO and Geometric GO approach

#########################################
####### Gradient Forest Method #########
#########################################

# Convert raster to data frame
raster2df=function(raster,pre_var,aggregate_factor=FALSE){ 
    raster_s = raster::subset(raster, pre_var)
    raster_c = crop(raster_s, range)
    raster_c = mask(raster_c, range)
    
    if (aggregate_factor){
        raster_c = aggregate(raster_c, aggregate_factor) # Adjust raster resolution
    }
    
    print(paste0("res:", res(raster_c)))
    
    raster_data = as.data.frame(rasterToPoints(raster_c))
    raster_env = na.omit(raster_data)
    print(dim(raster_env))
    return(raster_env)
}

# Build gradient forest model
run_gf <- function(Y, X, confounding_var=c()){
  nb_confound <- dim(confounding_var)[2]
  if (is.null(nb_confound)){
    confound_name <- c()
    df_gf <- data.frame(Y, X)
    colnames(df_gf) <- gsub("\\.", ":", colnames(df_gf))
    pred_name <- c(colnames(X))
  } else {
    confound_name <- colnames(confounding_var)
    df_gf <- data.frame(Y, X, confounding_var)
    colnames(df_gf) <- gsub("\\.", ":", colnames(df_gf))
    pred_name <- c(colnames(X), confound_name)
  }
  gf <- gradientForest(data = df_gf, predictor.vars = pred_name, response.vars = colnames(Y), ntree = 500, trace = T, check.names = FALSE)
}

# Predict genetic offset using GF model
gf_pred_surface <- function(gf, X_surface, X_surface_pred, pred_name){
  df_cur_var <- data.frame(X_surface)[,pred_name]
  df_fut_var <- data.frame(X_surface_pred)[,pred_name]
  currentcumimp <- predict(gf, df_cur_var)
  futurcumimp <- predict(gf, df_fut_var)
  genetic_offset <- sqrt(rowSums((currentcumimp - futurcumimp)^2))
  return(genetic_offset)
}


range=shapefile("/disk3/hkqm/shp/pop_sample_site_range/range_modify.shp")

# Load climate layers
current = raster::stack("./current/chelsa.v2.1.bio.tif")
names(current) = paste0("bio_", 1:19)

# 2070 ssp585 projections
IPSL_2070_585 = raster::stack("./2070/IPSL/chelsa.v2.1.2070.IPSL.ssp585.bio1_19.resample.tif")
MRI_2070_585 = raster::stack("./2070/MRI/chelsa.v2.1.2070.MRI.ssp585.bio1_19.resample.tif")
names(IPSL_2070_585) = names(MRI_2070_585) = paste0("bio_", 1:19)
mean_2070_585 = (IPSL_2070_585 + MRI_2070_585) / 2

# 2100 ssp585 projections
IPSL_2100_585 = raster::stack("./2100/IPSL/chelsa.v2.1.2100.IPSL.ssp585.bio1_19.resample.tif")
MRI_2100_585 = raster::stack("./2100/MRI/chelsa.v2.1.2100.MRI.ssp585.bio1_19.resample.tif")
names(IPSL_2100_585) = names(MRI_2100_585) = paste0("bio_", 1:19)
mean_2100_585 = (IPSL_2100_585 + MRI_2100_585) / 2

# Load genotype matrix 
load("./lfmm.012.RData") 
snps = read.table("./ai.txt", header = F) %>% as.matrix() %>% as.vector()
snp_geno = snp_geno[, snps]
snp_geno_df = data.frame(sample = rownames(snp_geno), snp_geno)
colnames(snp_geno_df) <- gsub("\\.", ":", colnames(snp_geno_df))

# Load sample information
all_info = read.csv("./info.csv", header = T)
rownames(all_info) = all_info$sample
sample_info = all_info[snp_geno_df$sample, c("sample", "long", "lat", "site_number", "sample_number")]

# Load population structure covariate (PC1)
eigenvec = read.table("./neu.pca.eigenvec", header = F, row.names = 1) %>% data.frame()
colnames(eigenvec) = paste0("PC_str_", 1:(ncol(eigenvec)))
eigenvec$sample = rownames(eigenvec)
pc_1 = eigenvec[snp_geno_df$sample, c("sample", "PC_str_1")]

sample_geno_info = merge(sample_info, snp_geno_df, "sample")
sample_geno_info = merge(sample_geno_info, pc_1, "sample")
sample_geno_info_filter = filter(sample_geno_info, sample_number >= 3)

pc1_mean = sample_geno_info_filter[, c("site_number", "sample_number", "PC_str_1")] %>%
  group_by(site_number) %>% summarise(sample_number = mean(sample_number), pc_1 = mean(PC_str_1))
condound_var = as.matrix(pc1_mean$pc_1)
colnames(condound_var) = "pc_1"


# Calculate allele frequency per site

info_data = sample_geno_info_filter[, 1:5] %>% distinct(long, .keep_all = TRUE)
allele_count = sample_geno_info_filter %>% group_by(site_number) %>% summarise_if(is.numeric, ~sum(.))
allele_count = allele_count[, c(-2, -3, -4)]
allele_count = inner_join(info_data, allele_count, "site_number")
allele_count = allele_count[, -1]
allele_fre = allele_count[, 5:ncol(allele_count)] / (allele_count$sample_number * 2)

# Extract environmental variables at sample sites
ind_current = raster::extract(current, allele_count[, c("long", "lat")])
ind_2070_585 = (raster::extract(IPSL_2070_585, allele_count[, c("long", "lat")]) + 
                raster::extract(MRI_2070_585, allele_count[, c("long", "lat")])) / 2
ind_2100_585 = (raster::extract(IPSL_2100_585, allele_count[, c("long", "lat")]) + 
                raster::extract(MRI_2100_585, allele_count[, c("long", "lat")])) / 2

# Prepare input matrices
Y = allele_fre
X = ind_current


# Select environmental predictors
keep.bio = c("bio_10", "bio_4", "bio_11", "bio_17")

# Train GF model
gf = run_gf(Y, X[, keep.bio], condound_var); gf

# Project future GF-based genetic offset
X.sur.current = raster2df(current, paste0("bio_", 1:19))
scenarios = list("2070.ssp585" = mean_2070_585, "2100.ssp585" = mean_2100_585)
gf_sur_raster_list = list()

for (name in names(scenarios)) {
  scenario_raster = scenarios[[name]]
  X.sur = raster2df(scenario_raster, paste0("bio_", 1:19))  
  gf_sur = data.frame(X.sur[, c("x", "y")], go = gf_pred_surface(gf, X.sur.current, X.sur, keep.bio))
  gf_sur_raster = rasterFromXYZ(gf_sur)
  gf_sur_raster_list[[name]] = gf_sur_raster
}


#########################################
####### Geometric Genetic Offset ########
#########################################

allele_count_geom = sample_info[, c("sample", "long", "lat")]
ind_current = raster::extract(current, allele_count_geom[, c("long", "lat")])
ind_2070_585 = (raster::extract(IPSL_2070_585, allele_count_geom[, c("long", "lat")]) + 
                raster::extract(MRI_2070_585, allele_count_geom[, c("long", "lat")])) / 2
ind_2100_585 = (raster::extract(IPSL_2100_585, allele_count_geom[, c("long", "lat")]) + 
                raster::extract(MRI_2100_585, allele_count_geom[, c("long", "lat")])) / 2

# Run geometric GO
for (a in 1:length(names(scenarios))){
  label = names(scenarios)[a]
  pred_env = list(ind_2070_585, ind_2100_585)[[a]]
  
  g.gap <- genetic.gap(input = Y, 
                       env = X[, keep.bio], 
                       new.env = ind_current[, keep.bio],
                       pred.env = pred_env[, keep.bio],
                       scale = TRUE,
                       K = 3)
  
  g.gap.data = cbind(sample_info[, c("sample", "long", "lat")], go = g.gap$offset)
  write.csv(g.gap.data, file = paste0("./result/g.go.", label, ".csv"), row.names = F, quote = F)
}
