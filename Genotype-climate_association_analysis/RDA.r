library(ggplot2)
library(raster)
library(RColorBrewer)
library(qvalue)
library(rnaturalearth)
library(rnaturalearthdata)
library(vegan)
library(sf)
library(dplyr)

# Load bioclimatic layers
current <- stack("/disk3/hkqm/climate.tif/chelsa2.1/resample/current/chelsa.v2.1.1981_2010.bio1_19.resample.tif")
names(current) <- paste0("bio_", 1:19)

# Load genotype data
load("./lfmm.012.RData")  # SNP genotype matrix from LFMM2 output
dim(snp_geno)

# Load sampling site information
info <- read.csv("./info.csv", header = TRUE)  # Sample coordinates (latitude, longitude)
samplename <- read.table("./sample.txt")       # Sample IDs in the same order as snp_geno columns
colnames(samplename) <- "sample"

# Match sample coordinates with sample IDs
sample_info <- inner_join(samplename, info)
# Extract environmental variables at sample locations
ind_clim <- raster::extract(current, sample_info[, 2:3])
ind_clim_data <- cbind(sample_info[, 2:3], ind_clim)        

# Function to identify outlier SNPs using RDA loadings
rdadapt <- function(rda, K){
  library(robust)
  library(qvalue)
  zscores <- as.matrix(rda$CCA$v[, 1:as.numeric(K)])        # Extract loadings of canonical axes
  resscale <- apply(zscores, 2, scale)                      # Standardize
  resmaha <- covRob(resscale, distance = TRUE, 
                    na.action = na.omit, estim = "pairwiseGK")$dist  # Mahalanobis distance
  lambda <- median(resmaha) / qchisq(0.5, df = K)          
  reschi2test <- pchisq(resmaha / lambda, K, lower.tail = FALSE)  # Chi-square test
  padj_BH <- p.adjust(reschi2test, method = "BH")           # FDR correction (Benjamini-Hochberg)
  qval <- qvalue(reschi2test)$qvalues                       # q-values
  padj <- p.adjust(reschi2test, method = "bonferroni")      # Bonferroni correction
  return(data.frame(p.values = reschi2test, 
                    padj_BH = padj_BH,
                    padj_q = qval,
                    padj = padj))
}

# Run RDA using selected environmental predictors
RDA <- rda(snp_geno ~ bio_11 + bio_10 + bio_4 + bio_17, data = ind_clim_data)

# Detect candidate loci using Mahalanobis distance on RDA loadings
rdadapt_res <- rdadapt(RDA, 2)

# Load SNP IDs
library(data.table)
snpid <- fread("./lfmm.vcfsnp", fill = TRUE, sep = " ", header = FALSE)[, 3]

# Identify outliers with FDR-adjusted p-value threshold
thres_env <- 0.01
outliers <- data.frame(
  Loci = snpid[which(rdadapt_res$padj_BH < thres_env)],
  p.value = rdadapt_res$p.values[which(rdadapt_res$padj_BH < thres_env)]
)