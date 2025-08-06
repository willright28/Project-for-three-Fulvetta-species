options(vsc.dev.args = list(width = 600, height = 500))
library(triangulaR)
library(vcfR)
library(dplyr)
######
sample=read.csv('./hybrid_test_pop_info.csv',header=T)
sample_mix=sample[,c("sample","label")]
colnames(sample_mix)=c("id","pop")

#############test 1##########################
sample_mix_ec=filter(sample_mix,pop!="fratercula")
write.table(sample_mix_ec$id,"./e_c.id",col.names=F,row.names = F,quote=F)
system("
plink2 --vcf /mnt/disk1/zhangshuai/hkqm/gghybrid/vcf71.8.60.maf.0.05.ld.nomissing.vcf   --aec --keep ./e_c.id --export vcf --out ./hueti_davidi")

data_hueti_davidi <- read.vcfR("/mnt/disk1/zhangshuai/hkqm/gghybrid/hueti_davidi.vcf", verbose = F)
example.vcfR.diff.1 <- alleleFreqDiff(vcfR =data_hueti_davidi , pm = sample_mix_ec, p1 = "davidi", p2 = "hueti", difference = 0.95)
hi.het.1 <- hybridIndex(vcfR = example.vcfR.diff.1 , pm = sample_mix_ec, p1 = "davidi", p2 = "hueti")

#pdf("/mnt/disk1/zhangshuai/hkqm/gghybrid/output.pdf", width = 7, height = 5) 
triangle.plot(hi.het.1,cex=3, max.overlaps = 1)
#dev.off() 


#############test 2##########################
sample_mix_wc=filter(sample_mix,pop!="hueti")
write.table(sample_mix_wc$id,"./w_c.id",col.names=F,row.names = F,quote=F)

system("
plink2 --vcf /mnt/disk1/zhangshuai/hkqm/gghybrid/vcf71.8.60.maf.0.05.ld.nomissing.vcf   --aec --keep ./w_c.id --export vcf --out ./fratercula_davidi")

data_fratercula_davidi <- read.vcfR("/mnt/disk1/zhangshuai/hkqm/gghybrid/fratercula_davidi.vcf", verbose = F)
example.vcfR.diff.2 <- alleleFreqDiff(vcfR = data_fratercula_davidi, pm = sample_mix_wc, p1 = "davidi", p2 = "fratercula", difference = 0.95)
hi.het.2 <- hybridIndex(vcfR = example.vcfR.diff.2, pm = sample_mix_wc,  p1 = "davidi", p2 = "fratercula")
hi.het.2
#pdf("/mnt/disk1/zhangshuai/hkqm/gghybrid/output2.pdf", width = 7, height = 5) 
triangle.plot(hi.het.2,cex=3, max.overlaps = 1)
#dev.off() 
