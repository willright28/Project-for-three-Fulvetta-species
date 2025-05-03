library(raster)
library(vcfR)
library(dplyr)
library(tidyr)
library(adegenet)
library(AlleleShift)
library(plyr)


##################################################
## Load adaptive introgressed loci
intro=read.vcfR(file="./alleleshift/ai.vcf") #VCF containing adaptive introgressed loci
intro.genind <- vcfR2genind(intro)
info=read.csv("./hkqm/alleleshift/info.csv") #sample information
rownames(info)=info$sample
info=info[rownames(intro.genind@tab),] 
pop(intro.genind)=info$pop
intro.genpop.pop=adegenet::genind2genpop(intro.genind)


## Calculate current allele frequency per locus and assign to different introgression events
fre=makefreq(intro.genpop.pop)%>%data.frame()
fre_t=t(fre)%>%data.frame()
fre_t$id=rownames(fre_t)
fre_t=tidyr::separate(data=fre_t,col=id,into=c("chr","pos","allele"),sep="\\.")
fre_t=tidyr::unite(fre_t,"loci",6:7,sep=":")

fin_res=read.csv('./filet_res.csv',header = T)[,c("loci","status","set")]  ### Info on adaptive introgressed loci from FILET analysis
fre_t_join=left_join(fre_t,fin_res,)
head(fre_t_join)

# davidi to hueti
hueti_davidi_1=dplyr::filter(fre_t_join, status == 1 & set == "hueti_davidi"&  davidi > hueti & davidi > fratercula )
hueti_davidi_1_id=hueti_davidi_1[,c("loci","allele")]
# huetiern to davidiral
hueti_davidi_2=dplyr::filter(fre_t_join, status == 2 & set == "hueti_davidi"  & hueti > davidi & hueti > fratercula )
hueti_davidi_2_id=hueti_davidi_2[,c("loci","allele")]
# davidiral to fraterculaern
fratercula_davidi_1=dplyr::filter(fre_t_join, status == 1 & set == "fratercula_davidi"  & davidi > hueti & davidi > fratercula)
fratercula_davidi_1_id=fratercula_davidi_1[,c("loci","allele")]
# fraterculaern to davidiral
fratercula_davidi_2=dplyr::filter(fre_t_join, status == 2 & set == "fratercula_davidi" & fratercula > davidi &  fratercula > hueti )
fratercula_davidi_2_id=fratercula_davidi_2[,c("loci","allele")]


########################################
pop(intro.genind)=info$pop_site_number
intro.genpop=adegenet::genind2genpop(intro.genind)

#Environmental variables
pop_site=unique(info[,c("long","lat","pop_site_number")])
rownames(pop_site)=pop_site$pop_site_number
pop_site=pop_site[rownames(intro.genpop@tab),]

baseline.stack=stack("./chelsa.v2.1.tif")
names(baseline.stack)=paste0("bio_",1:19)
baseline.env=raster::extract(baseline.stack,pop_site[,c("long","lat")],df=TRUE)
rownames(baseline.env)=rownames(pop_site)

#ssp585#
# 2070
# MRI
mri_2070_585=stack("./2070/MRI/chelsa.v2.1.2070.MRI.ssp585.bio1_19.resample.tif")
names(mri_2070_585)=paste0(rep("bio_",19),c(1:19))
mri_2070_585_env=raster::extract(mri_2070_585,pop_site[,c("long","lat")],df=TRUE)
rownames(mri_2070_585_env)=rownames(pop_site)

# IPSL
ipsl_2070_585=stack("./2070/IPSL/chelsa.v2.1.2070.IPSL.ssp585.bio1_19.resample.tif")
names(ipsl_2070_585)=paste0(rep("bio_",19),c(1:19))
ipsl_2070_585_env=raster::extract(ipsl_2070_585,pop_site[,c("long","lat")],df=TRUE)
rownames(ipsl_2070_585_env)=rownames(pop_site)

f2070.env=(mri_2070_585_env+ipsl_2070_585_env)/2

# 2100
# MRI
mri_2100_585=stack("./2100/MRI/chelsa.v2.1.2100.MRI.ssp585.bio1_19.resample.tif")
names(mri_2100_585)=paste0(rep("bio_",19),c(1:19))
mri_2100_585_env=raster::extract(mri_2100_585,pop_site[,c("long","lat")],df=TRUE)
rownames(mri_2100_585_env)=rownames(pop_site)

# IPSL
ipsl_2100_585=stack("./2100/IPSL/chelsa.v2.1.2100.IPSL.ssp585.bio1_19.resample.tif")
names(ipsl_2100_585)=paste0(rep("bio_",19),c(1:19))
ipsl_2100_585_env=raster::extract(ipsl_2100_585,pop_site[,c("long","lat")],df=TRUE)
rownames(ipsl_2100_585_env)=rownames(pop_site)

f2100.env=(mri_2100_585_env+ipsl_2100_585_env)/2

mean(f2100.env$bio_10)
mean(baseline.env$bio_10)

bio=c("bio_10","bio_11","bio_17","bio_4")


alcippe=read.csv("./info.csv",header = T)
alcippe=unique(alcippe[,c("long","lat")])

occ.current.env=raster::extract(baseline.stack,alcippe)
occ.2070.585.env=(raster::extract(mri_2070_585,alcippe)+raster::extract(ipsl_2070_585,alcippe))/2
occ.2100.585.env=(raster::extract(mri_2100_585,alcippe)+raster::extract(ipsl_2100_585,alcippe))/2

baseline.env.data=occ.current.env[,bio]
future.env.data=occ.2070.585.env[,bio]
future.env.2.data=occ.2100.585.env[,bio]

###############################################################
# Model and predict allele frequencies

lfmm.count.model <- count.model(intro.genpop, 
                                env.data=baseline.env[,bio], 
                                ordistep=F)
lfmm.pred.baseline <- count.pred(lfmm.count.model, env.data=baseline.env[,bio])

lfmm.freq.model <- freq.model(lfmm.pred.baseline)
lfmm.freq.baseline <- freq.pred(lfmm.freq.model,
                                count.predicted=lfmm.pred.baseline)

# Difference between predicted and observed values
plotA1 <- freq.ggplot(lfmm.freq.baseline, plot.best=TRUE,ylim=c(0.0, 0.8))

# Predict for future
# 2070
lfmm.pred.2070 <- count.pred(lfmm.count.model, env.data=f2070.env[,bio])
lfmm.freq.2070 <- freq.pred(lfmm.freq.model, count.predicted=lfmm.pred.2070)
times=nrow(lfmm.freq.2070)/length(unique(lfmm.freq.2070$Pop))

# 2100
lfmm.pred.2100 <- count.pred(lfmm.count.model, env.data=f2100.env[,bio])

lfmm.freq.2100 <- freq.pred(lfmm.freq.model,
                            count.predicted=lfmm.pred.2100)
