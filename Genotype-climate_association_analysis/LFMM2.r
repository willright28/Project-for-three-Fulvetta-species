library(LEA)
library(doParallel)
library(foreach)
library(data.table)
library(xlsx)
library(raster)
library(dplyr)

##### Convert VCF to geno format
output = vcf2geno("./maf.0.05.vcf", output.file = "./lfmm.geno")
snp_geno = read.geno("./lfmm.geno")

# Set SNP IDs
snpid = fread("./lfmm.vcfsnp", fill = TRUE, sep = " ", header = FALSE)[, 3]
colnames(snp_geno) = as.vector(as.matrix(snpid))

# Set sample IDs
system("bcftools query -l ./maf.0.05.vcf > ./sample.txt")
samplename = read.table("./sample.txt")
rownames(snp_geno) = as.vector(as.matrix(samplename))

# Load environmental variables
current = stack(".CHELSA_V.2.1.bio.tif")
names(current) = paste0("bio_", 1:19)
current = subset(current, c("bio_10", "bio_4", "bio_11", "bio_17"))

# Sample latitude and longitude
all_info = read.csv("./info.csv", header = TRUE)
samplename = read.table("./sample.txt")
colnames(samplename) = "sample"
sample_info = inner_join(samplename, all_info)
ind_clim = raster::extract(current, sample_info[, 2:3])

#### Define LFMM running function
lfmmfunction <- function(i, input_lfmm, K) {
  environ <- as.matrix(environments[, i])
  mod <- lfmm2(input = input_lfmm, env = environ, K = K)
  pv <- lfmm2.test(object = mod, input = input_lfmm, env = environ, linear = TRUE)
  p = pv$pvalues
  padj_BH <- p.adjust(pv$pvalues, method = "BH")
  padj <- p.adjust(pv$pvalues, method = "bonferroni")
  pv_df = data.frame(pos = c(1:dim(snp_geno)[2]), p = p, padj_BH = padj_BH, padj = padj, row.names = NULL)
  write.table(pv_df, paste("./lfmm/vcf71_keep/", K, ".biov_colname_", i, ".txt", sep = ""), row.names = FALSE, quote = FALSE)
}

# Parameters
environments = ind_clim
K = 3
p = 0.01

### Run in parallel
cl <- makeCluster(20)
registerDoParallel(cl)
setwd("./lfmm/vcf71_keep")
lfmm2res <- foreach(i = c(1:4), .combine = "c", .packages = "LEA") %dopar% lfmmfunction(i, snp_geno, K)
stopCluster(cl)
closeAllConnections()
