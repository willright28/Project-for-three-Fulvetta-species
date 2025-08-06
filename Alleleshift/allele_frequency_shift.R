library(vcfR)
library(poppr)
library(dplyr)
library(raster)
library(AlleleShift)
library(mgcv)
library(ggplot2)
library(reshape2)

################################################################
#########################allele frequency shift#########################
################################################################

#Load vcf file
intro=read.vcfR(file="./filet_intro.vcf")
intro.genind <- vcfR2genind(intro)
#Load sampling information
info=read.csv("./pop_info.csv") 
rownames(info)=info$sample
info=info[rownames(intro.genind@tab),]


#####Filter and identify alleles specific to different directions of gene flow#####
pop(intro.genind)=info$mix2 #population grouping information
intro.genpop.pop=adegenet::genind2genpop(intro.genind)
fre=makefreq(intro.genpop.pop)%>%data.frame() #calculate Allele frequency
fre_t=t(fre)%>%data.frame()
fre_t$id=rownames(fre_t)
fre_t=tidyr::separate(data=fre_t,col=id,into=c("chr","pos","allele"),sep="\\.")
fre_t=tidyr::unite(fre_t,"loci",6:7,sep=":")

fin_res=read.csv('./info_filet_all.csv',header = T) #derived from FILET results and contains information on introgression loci
fin_res=fin_res[,c("loci","direction","pair")]
fre_t_join=left_join(fre_t,fin_res,relationship = "many-to-many")
dim(fre_t_join)/2


#Loci introgressed from davidi to hueti
hueti_davidi_1=dplyr::filter(fre_t_join, direction == 1 & pair == "hueti_davidi" & 
                               davidi_pure > 0.9 & hueti_pure < 0.1 & hueti_hybrid > hueti_pure & davidi_pure > fratercula_pure )

hueti_davidi_1=hueti_davidi_1[,c("loci","allele")]
dim(hueti_davidi_1)

#Loci introgressed from hueti to davidi
hueti_davidi_2=dplyr::filter(fre_t_join, direction == 2 & pair == "hueti_davidi" & 
                               hueti_pure > 0.9 & davidi_pure < 0.1 & davidi_hybrid > davidi_pure & hueti_pure > fratercula_pure )

hueti_davidi_2=hueti_davidi_2[,c("loci","allele")]
dim(hueti_davidi_2)

#Loci introgressed from davidi to fratercula
fratercula_davidi_1=dplyr::filter(fre_t_join, direction == 1 & pair == "fratercula_davidi" & 
                                    davidi_pure > 0.9 & fratercula_pure < 0.1 & davidi_hybrid > fratercula_pure & davidi_pure > hueti_pure )
fratercula_davidi_1=fratercula_davidi_1[,c("loci","allele")]
dim(fratercula_davidi_1)

#Loci introgressed from fratercula to davidi
fratercula_davidi_2=dplyr::filter(fre_t_join, direction == 2 & pair == "fratercula_davidi" & 
                                    fratercula_pure > 0.9 & davidi_pure < 0.1 & davidi_hybrid > davidi_pure & fratercula_pure >  hueti_pure )
fratercula_davidi_2=fratercula_davidi_2[,c("loci","allele")]
dim(fratercula_davidi_2)

intro_allele=rbind(hueti_davidi_1,hueti_davidi_2,fratercula_davidi_1,fratercula_davidi_2)
View(intro_allele)

#####Predict allele frequency#####
pop(intro.genind)=info$site_number #use local sampling population label
intro.genpop=adegenet::genind2genpop(intro.genind)
pop_site=unique(info[,c("long","lat","site_number")])
rownames(pop_site)=pop_site$site_number
pop_site=pop_site[rownames(intro.genpop$tab),]
pop_site
intro.genind

#Load current biocliam layers
baseline.stack=stack("D:/climate.resample/resample/current/chelsa.v2.1.1981_2010.bio1_19.resample.tif")
names(baseline.stack)=paste0("bio_",1:19)
baseline.env=raster::extract(baseline.stack,pop_site[,c("long","lat")],df=TRUE)
rownames(baseline.env)=rownames(pop_site)

#Load future biocliam layers(ssp585)
#use two GCM here
#2070
#MRI
mri_2070_585=stack("D:/climate.resample/resample/2070/MRI/chelsa.v2.1.2070.MRI.ssp585.bio1_19.resample.tif")
names(mri_2070_585)=paste0(rep("bio_",19),c(1:19))
mri_2070_585_env=raster::extract(mri_2070_585,pop_site[,c("long","lat")],df=TRUE)
rownames(mri_2070_585_env)=rownames(pop_site)
# IPSL
ipsl_2070_585=stack("D:/climate.resample/resample/2070/IPSL/chelsa.v2.1.2070.IPSL.ssp585.bio1_19.resample.tif")
names(ipsl_2070_585)=paste0(rep("bio_",19),c(1:19))
ipsl_2070_585_env=raster::extract(ipsl_2070_585,pop_site[,c("long","lat")],df=TRUE)
rownames(ipsl_2070_585_env)=rownames(pop_site)

f2070.env=(mri_2070_585_env+ipsl_2070_585_env)/2

#2100
#MRI
mri_2100_585=stack("D:/climate.resample/resample/2100/MRI/chelsa.v2.1.2100.MRI.ssp585.bio1_19.resample.tif")
names(mri_2100_585)=paste0(rep("bio_",19),c(1:19))
mri_2100_585_env=raster::extract(mri_2100_585,pop_site[,c("long","lat")],df=TRUE)
rownames(mri_2100_585_env)=rownames(pop_site)

#ISPL
ipsl_2100_585=stack("D:/climate.resample/resample/2100/IPSL/chelsa.v2.1.2100.IPSL.ssp585.bio1_19.resample.tif")
names(ipsl_2100_585)=paste0(rep("bio_",19),c(1:19))
ipsl_2100_585_env=raster::extract(ipsl_2100_585,pop_site[,c("long","lat")],df=TRUE)
rownames(ipsl_2100_585_env)=rownames(pop_site)

f2100.env=(mri_2100_585_env+ipsl_2100_585_env)/2

mean(f2100.env$bio_10)
mean(baseline.env$bio_10)

bio=c("bio_10","bio_11","bio_17","bio_4") #specify which predictor variables are used

#alleleshift model
###
count.model <- count.model(intro.genpop, 
                                env.data=baseline.env[,bio], 
                                ordistep=F,cca.model =F)
pred.baseline <- count.pred(count.model, env.data=baseline.env[,bio])
###
freq.model <- freq.model(pred.baseline)
freq.baseline <- freq.pred(freq.model,
                                count.predicted=pred.baseline)

plotA1 <- freq.ggplot(freq.baseline, 
                      plot.best=TRUE,
                      ylim=c(0.0, 0.8)) #select population with R2 > 0.5


#alleleshift predict
####2070
pred.2070 <- count.pred(count.model, env.data=f2070.env[,bio])

freq.2070 <- freq.pred(freq.model,
                            count.predicted=pred.2070)


freq.2070=tidyr::separate(freq.2070,Allele,c("chr","pos","status"),sep="\\.")
freq.2070=tidyr::unite(freq.2070,"id",c("chr","pos"), sep=":", remove = T)
freq.2070=arrange(freq.2070,id)

n=nrow(freq.2070)/length(unique(freq.2070$Pop))
n
res.2070.list=list()

for(i in c(1:n)){
  res.2070.list[i]=list(freq.2070[(c(1:length(unique(freq.2070$Pop)))+13*(i-1)),][,c("id","Pop","status","Allele.freq","Bp","Ap","Freq.e2","increasing")])
  res.2070.list[[i]]=arrange(res.2070.list[[i]],Pop)
  names(res.2070.list)[i]=unique(res.2070.list[[i]]$id)
}
length(res.2070.list)

res.2070.list=res.2070.list[intro_allele[,1]]
res.2070.list[1]

for(i in 1:length(res.2070.list)){
  print(i)
  name=names(res.2070.list)[i]
  row.names(res.2070.list[[i]])=res.2070.list[[i]]$Pop
  status=res.2070.list[[i]]$status%>%unique()%>%as.integer();status
  status_1=intro_allele[i,2]%>%as.integer();status_1
  if( status_1 !=  status ){  
    data=data.frame(id=1:13)
    data$id=res.2070.list[[i]][,"id"]
    data$Pop=res.2070.list[[i]][,"Pop"]
    data$status=abs(1-status)
    data$Allele.freq=1-res.2070.list[[i]][,"Allele.freq"]
    data$Ap=res.2070.list[[i]][,"Bp"]
    data$Bp=res.2070.list[[i]][,"Ap"]
    data$Freq.e2=1-res.2070.list[[i]][,"Freq.e2"]
    data$increasing=!res.2070.list[[i]][,"increasing"]
    res.2070.list[[i]]=data
    names( res.2070.list)[i]=name
  }
  res.2070.list[[i]]$time="2070"
}


####2100
pred.2100 <- count.pred(count.model, env.data=f2100.env[,bio])

freq.2100 <- freq.pred(freq.model,
                       count.predicted=pred.2100)


freq.2100=tidyr::separate(freq.2100,Allele,c("chr","pos","status"),sep="\\.")
freq.2100=tidyr::unite(freq.2100,"id",c("chr","pos"), sep=":", remove = T)
freq.2100=arrange(freq.2100,id)

n=nrow(freq.2100)/length(unique(freq.2100$Pop))
n
res.2100.list=list()

for(i in c(1:n)){
  res.2100.list[i]=list(freq.2100[(c(1:length(unique(freq.2100$Pop)))+13*(i-1)),][,c("id","Pop","status","Allele.freq","Bp","Ap","Freq.e2","increasing")])
  res.2100.list[[i]]=arrange(res.2100.list[[i]],Pop)
  names(res.2100.list)[i]=unique(res.2100.list[[i]]$id)
}
length(res.2100.list)

res.2100.list=res.2100.list[intro_allele[,1]]
res.2100.list[1]

for(i in 1:length(res.2100.list)){
  print(i)
  name=names(res.2100.list)[i]
  row.names(res.2100.list[[i]])=res.2100.list[[i]]$Pop
  status=res.2100.list[[i]]$status%>%unique()%>%as.integer();status
  status_1=intro_allele[i,2]%>%as.integer();status_1
  if( status_1 !=  status ){  
    data=data.frame(id=1:13)
    data$id=res.2100.list[[i]][,"id"]
    data$Pop=res.2100.list[[i]][,"Pop"]
    data$status=abs(1-status)
    data$Allele.freq=1-res.2100.list[[i]][,"Allele.freq"]
    data$Ap=res.2100.list[[i]][,"Bp"]
    data$Bp=res.2100.list[[i]][,"Ap"]
    data$Freq.e2=1-res.2100.list[[i]][,"Freq.e2"]
    data$increasing=!res.2100.list[[i]][,"increasing"]
    res.2100.list[[i]]=data
    names( res.2100.list)[i]=name
  }
  res.2100.list[[i]]$time="2100"
}


#Summarize allele frequencies across the three time periods for each locus 

af_summ=data.frame()

for(i in (1:length(res.2100.list)) ){
  res.lgm.list.mod.12=list(res.2070.list[[i]],res.2100.list[[i]])
  res.lgm.list.mod.12.df=do.call(rbind,res.lgm.list.mod.12)
  res.lgm.list.mod.12.current=res.lgm.list.mod.12.df
  res.lgm.list.mod.12.current=res.lgm.list.mod.12.current[,c("Pop","Allele.freq","id","time","increasing","status")]
  res.lgm.list.mod.12.current$time="Current"
  res.lgm.list.mod.12.current$increasing=NA
  colnames(res.lgm.list.mod.12.current)=c("Pop","Freq.e2","id","time","increasing","status")
  res.lgm.list.mod.12.df=res.lgm.list.mod.12.df[,c("Pop","Freq.e2","id","time","increasing","status")]
  res.lgm.list.mod.12.df$status=as.numeric(res.lgm.list.mod.12.df$status)
  res.lgm.list.mod.12.df=rbind(res.lgm.list.mod.12.current,res.lgm.list.mod.12.df)
  res.lgm.list.mod.12.df$time=factor(res.lgm.list.mod.12.df$time,level=c("Current","2070","2100"),labels=c("2000","2070","2100"))
  label_text=unique(res.lgm.list.mod.12.df$Allele)
  
  #mean allele frequency
  res.lgm.list.mod.12.df.mod=res.lgm.list.mod.12.df%>%group_by(time)%>%dplyr::summarise(m_freq=mean(Freq.e2),id=id,status=status,time=time)%>%unique()
  res.lgm.list.mod.12.df.mod$status=as.numeric(res.lgm.list.mod.12.df.mod$status)
  res.lgm.list.mod.12.df.mod.m=dcast( res.lgm.list.mod.12.df.mod%>%as.data.frame(),id + status ~time,value.var="m_freq")
  af_summ=rbind( af_summ,res.lgm.list.mod.12.df.mod.m)
}

summary(af_summ)

af_summ=reshape2::melt(data=af_summ,id.vars=c("id","status"),measure.vars=c("2000","2070","2100"),variable.name="time",value.name="m_freq")
af_summ$time=as.character(af_summ$time)
af_summ$time=as.numeric(af_summ$time)

###GAM fit

ggplot()+
geom_smooth(data = af_summ,
              aes(x = time, y = m_freq),
              method = "gam",
              formula = y ~ s(x,k = 3),     
              size = 2,
              color = "#7A96A7",
              fill = "#7A96A7",
              alpha = 0.2,level=0.9)+
  theme_bw(base_size = 20) +
  theme(aspect.ratio=3/4)+
  scale_x_continuous(limits = c(2000,2100),breaks = seq(2000,2100,50))+
  coord_cartesian(ylim = c(0.38,0.48))+
  
  scale_y_continuous(breaks = seq(0.38, 0.48, 0.02)) +  
  coord_cartesian(ylim = c(0.38, 0.48)) +  
  ylab("Allele frequency")+xlab("")

################################################################
#########################alleles absent#########################
################################################################
#load vcf file
intro=read.vcfR(file="./filet_intro.vcf")
intro.genind <- vcfR2genind(intro)
info=read.csv("./pop_info.csv") 
rownames(info)=info$sample
info=info[rownames(intro.genind@tab),] 
pop(intro.genind)=info$mix
intro.genpop=adegenet::genind2genpop(intro.genind)

#calculate allele frequency
fre=makefreq(intro.genpop)%>%data.frame()

count.model <- count.model(intro.genpop, 
                                env.data=baseline.env[,bio]) 
pred.baseline <- count.pred(count.model, env.data=baseline.env[,bio])
pred.baseline

#glm
freq.model <- freq.model(pred.baseline)
freq.baseline <- freq.pred(freq.model,
                                count.predicted=pred.baseline)

pred.2100 <- count.pred(count.model, env.data=f2100.env[,bio])

freq.2100 <- freq.pred(freq.model,
                            count.predicted=pred.2100)

data_2100=filter(freq.2100, (A==0 | B==0) & (Ap>=0 & Bp>=0) )

summ=table(data_2100$Pop)
#calculate for each species and hybrid base on result 'summ'
hybrid_2100=(1378+70+245+220)/4/(2458*2)*100;hybrid_2100
davidi_2100=(1702+1249+2291+1809)/4/(2458*2)*100;davidi_2100
hueti_2100=(909+1054+492)/3/(2458*2)*100;hueti_2100
fratercula_2100=(1050+2009)/2/(2458*2)*100;fratercula_2100

#######################
data <- data.frame(
  sp = c("fratercula", "davidi", "hueti","Hybrids"),
  Value = c(31.12,35.86,16.65,9.73))

data$sp=factor(data$sp,level=c("fratercula", "davidi", "hueti","Hybrids"))

library(scales)
ggplot(data, aes(x = sp, y = Value,fill=sp)) +
  scale_fill_manual(values=c("#fdc695","#b9ccef", "#9ae5cf","#f7a4a4"))+
  geom_bar(position = position_dodge(width = 0.6), stat = "identity", width = 0.6) +
  labs(
    x = "Categories",
    y = "Values",
    fill = NULL
  ) +
  theme_classic()+xlab(NULL)+ylab("Percentage of alleles absent")+theme(text = element_text(size = 20))




