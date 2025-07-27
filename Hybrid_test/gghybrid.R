library("devtools")
library("gghybrid")
library("coda")

info=read.csv("/mnt/disk1/zhangshuai/hkqm/gghybrid/hkqm.info.2.csv",header=T)
id=read.table("/mnt/disk1/zhangshuai/hkqm/gghybrid/id.txt")
rownames(info)=info$sample
info=info[id$V1,]
info_pop=info[,c("site_pop","sample")]
info_pop
write.table(info_pop,"/mnt/disk1/zhangshuai/hkqm/gghybrid/info_pop.csv",col.names=F,row.names=F,quote=F,sep=" ")

data1 <- read.data(file="/mnt/disk1/zhangshuai/hkqm/gghybrid/hkqm.west_central.recode.58.strct_in",nprecol=2,NUMINDS=58,MISSINGVAL=NA,NUMLOCI =151255,ONEROW=1,precol.headers=0,MAPDISTANCE=1) 
save(data1,file="/mnt/disk1/zhangshuai/hkqm/gghybrid/data1.RData")
str(data1$data$POPID)

str(data1)
prepdata= data.prep(
 data=data1$data,
 loci=data1$loci,
 sourceAbsent = FALSE,        #default is FALSE, for the existence of parental reference samples in the dataset#
 #marker.info=dat$marker.info,#this time we're not uploading any extra marker info#
 alleles=data1$alleles,
 S0=c("11","14"),                                    #first parental reference set; must be specified when sourceAbsent = FALSE#
S1=c("8","13","3","4","10","15"),                         #second parental reference set; must be specified when sourceAbsent = FALSE#
 precols=data1$precols,
#  return.genotype.table=TRUE,  #This returns an unmelted table of diploid (or other ploidy) genotypes, coded as number of copies of the allele with higher frequency in S1#
#  return.locus.table=TRUE      #Returns a table with one row per locus and all the individual marker data, including those already uploaded plus values calculated within this function#
)

prepdata

prepdata_1$sourceInfo

hindlabel_1=esth(
 data.prep.object=prepdata_1$data.prep,
 read.data.precols=data1$precols,
 include.Source=TRUE,	                 #Leave at default TRUE if you want hybrid indices for the parental reference individuals, which is often useful#
plot.ind = c(1,2,5,6,7,9,12),  #Optionally plot some individuals in real time. Merely shows how well the adaptive burnin is working#
plot.col = c("blue", "green", "cyan", "purple", "magenta", "red", "orange"),
 nitt=3000,                              #Testing suggests nitt=3000,burnin=1000 are sufficient for accurate posterior estimation#
 burnin=1000
)
setkey(hindlabel_1$hi,beta_mean)

hindlabel_1
getwd()
hist(hindall_1$h_gelman)

chr=c("1","2","3")

hindlabel_1$hi$Source=c(rep("T",58))


abc = plot_h(data=hindlabel_1$hi,
 test.subject=hindlabel_1$test.subject,
 mean.h.by="POPID",			             #Calculate the mean hybrid index for each value of the "POPID" column#
 sort.by=c("mean_h","POPID","h_posterior_mode"),  #Order test subjects along the x axis by the mean hybrid index calculated above and also by individual hybrid index ("POPID" is included as some population pairs may have identical mean hi).
 col.group="POPID",
 group.sep="POPID",
 fill.source=FALSE,
 basic.lines=FALSE,
 source.col=c("blue","red"),
 source.limits=c("blue","red"),
 #custom.abline=abline(h=0.76,col="grey",lty=2),  #The lowest hybrid index estimate for (Sardinian) Spanish sparrows#
 cex=1,pch=16,
 cex.lab=1.5,cex.main=1.5,ylim=c(0,1))

#########################################


data2 <- read.data(file="/mnt/disk1/zhangshuai/hkqm/gghybrid/hkqm.east_central.recode.53.strct_in",nprecol=2,NUMINDS=53,MISSINGVAL=NA,NUMLOCI =151255,ONEROW=1,precol.headers=0,MAPDISTANCE=1) 
save(data2,file="/mnt/disk1/zhangshuai/hkqm/gghybrid/data2.RData")

load("/mnt/disk1/zhangshuai/hkqm/gghybrid/data2.RData")

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



####filter




