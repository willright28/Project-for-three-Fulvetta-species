library("devtools")
library("gghybrid")
library("coda")

info_pop=read.csv("./sampling_site_pop_info.csv",header=F)

#Convert the VCF file into the gghybrid input format
data2 <- read.data(file="/mnt/disk1/zhangshuai/hkqm/gghybrid/hkqm.east_central.recode.53.strct_in",nprecol=2,NUMINDS=53,MISSINGVAL=NA,NUMLOCI =151255,ONEROW=1,precol.headers=0,MAPDISTANCE=1) 

prepdata_2= data.prep(
 data=data2$data,
 loci=data2$loci,
 sourceAbsent = FALSE,        #default is FALSE, for the existence of parental reference samples in the dataset#
 #marker.info=dat$marker.info,#this time we're not uploading any extra marker info#
 alleles=data2$alleles,
 S0=c(12,16,6,7,14,17),                                    #first parental reference set; must be specified when sourceAbsent = FALSE#
 S1=c(8,1,3,2,18),                          #second parental reference set; must be specified when sourceAbsent = FALSE#
 precols=data2$precols,
 return.genotype.table=TRUE,  #This returns an unmelted table of diploid (or other ploidy) genotypes, coded as number of copies of the allele with higher frequency in S1#
 return.locus.table=TRUE      #Returns a table with one row per locus and all the individual marker data, including those already uploaded plus values calculated within this function#
)


hindlabel_2=esth(
 data.prep.object=prepdata_2$data.prep,
 read.data.precols=data2$precols,
 include.Source=TRUE,	                 #Leave at default TRUE if you want hybrid indices for the parental reference individuals, which is often useful#
plot.ind = c(4,5,9,10,11,13,15),  #Optionally plot some individuals in real time. Merely shows how well the adaptive burnin is working#
plot.col = c("blue", "green", "cyan", "purple", "magenta", "red", "orange"),
 nitt=3000,                              #Testing suggests nitt=3000,burnin=1000 are sufficient for accurate posterior estimation#
 burnin=1000
)





