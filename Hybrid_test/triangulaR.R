library(triangulaR)
# library(vcfR)
# devtools::install_github("omys-omics/triangulaR")
# Read in data
######
sample=read.csv('/mnt/disk1/zhangshuai/hkqm/genetic_offset/pop.0.9.csv',header=T)
sample=sample[sample$sample %in% ind.71,]
head(sample)
sample_mix=sample[,c("sample","mix")]
colnames(sample_mix)=c("id","pop")
sample_mix
#####
sample_mix_ec=filter(sample_mix,pop!="Western")
write.table(sample_mix_ec$id,"/mnt/disk1/zhangshuai/hkqm/af_shifts/e_c.id",col.names=F,row.names = F,quote=F)
sample_mix_wc=filter(sample_mix,pop!="Eastern")
write.table(sample_mix_wc$id,"/mnt/disk1/zhangshuai/hkqm/af_shifts/w_c.id",col.names=F,row.names = F,quote=F)


system("
plink2 --vcf /mnt/disk1/zhangshuai/hkqm/gghybrid/vcf71.8.60.maf.0.05.ld.nomissing.vcf   --aec --keep /mnt/disk1/zhangshuai/hkqm/af_shifts/e_c.id --export vcf --out /mnt/disk1/zhangshuai/hkqm/gghybrid/ec
 ")

system("
plink2 --vcf /mnt/disk1/zhangshuai/hkqm/gghybrid/vcf71.8.60.maf.0.05.ld.nomissing.vcf   --aec --keep /mnt/disk1/zhangshuai/hkqm/af_shifts/w_c.id --export vcf --out /mnt/disk1/zhangshuai/hkqm/gghybrid/wc")

options(vsc.dev.args = list(width = 600, height = 500))


data_ec <- read.vcfR("/mnt/disk1/zhangshuai/hkqm/gghybrid/ec.vcf", verbose = F)
data_ec
# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
example.vcfR.diff_ec <- alleleFreqDiff(vcfR = data_ec , pm = sample_mix_ec, p1 = "Central", p2 = "Eastern", difference = 0.95)

hi.het.ec <- hybridIndex(vcfR = example.vcfR.diff_ec, pm = sample_mix_ec, p1 = "Central", p2 = "Eastern")
hi.het.ec
pdf("/mnt/disk1/zhangshuai/hkqm/gghybrid/output.pdf", width = 7, height = 5) 
triangle.plot(hi.het.ec,cex=3, max.overlaps = 1)
dev.off() 


data_wc <- read.vcfR("/mnt/disk1/zhangshuai/hkqm/gghybrid/wc.vcf", verbose = F)
data_wc
sample_mix_wc
# Create a new vcfR object composed only of sites above the given allele frequency difference threshold
example.vcfR.diff_wc <- alleleFreqDiff(vcfR = data_wc, pm = sample_mix_wc, p1 = "Central", p2 = "Western", difference = 0.95)

hi.het.wc <- hybridIndex(vcfR = example.vcfR.diff_wc, pm = sample_mix_wc,  p1 = "Central", p2 = "Western")
hi.het.wc
pdf("/mnt/disk1/zhangshuai/hkqm/gghybrid/output2.pdf", width = 7, height = 5) 
triangle.plot(hi.het.wc,cex=3, max.overlaps = 3)
dev.off() 
