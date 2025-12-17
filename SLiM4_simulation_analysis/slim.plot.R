library(dplyr)
library(ggplot2)
library(ggnewscale)
library(patchwork)
# library(extrafont)
# loadfonts(device = "pdf")
burnin=20000

####################################
#########fitness###################
####################################

#####Gene flow model fitness data####
file_list=list.files("D:/slim4/slim_code/slim/1_1.5_1.5/geneflow",full.names = T,pattern ="^pop" )
geneflow_data=data.frame(times=c(),fitness=c(),run=c())
gf_fitness=data.frame()
for ( i in 1:length(file_list)){
  gf_fitness_df=read.table(file_list[i])
  colnames(gf_fitness_df)=c("times","pop","fitness","size","phenotypes","n","GV")
  gf_fitness_df$pop=factor(gf_fitness_df$pop)
  gf_fitness_df$run=i
  gf_fitness=rbind(gf_fitness,gf_fitness_df)
  gf_fitness_df_mean=gf_fitness_df%>%group_by(times)%>%dplyr::summarise(fitness=mean(fitness))
  gf_fitness_df_mean= gf_fitness_df_mean %>% dplyr::summarise(times=times,fitness=fitness/max(fitness))
  gf_fitness_df_mean$run=i
  geneflow_data=rbind(geneflow_data,gf_fitness_df_mean)
}

gf_fitness$type="Gene flow model"
gf_fitness_mean=gf_fitness %>% group_by(pop,times) %>% dplyr::summarise(fitness=mean(fitness))

geneflow_data$group="Gene flow model"
geneflow_data_mean=geneflow_data%>%group_by(times)%>%dplyr::summarise(times=mean(times),fitness=mean(fitness))
geneflow_data_mean$group="Gene flow model"


#####Non-gene flow model fitness data#####
file_list=list.files("D:/slim4/slim_code/slim/1_1.5_1.5/no_geneflow/",full.names = T,pattern ="^pop")
no_geneflow_data=data.frame(times=c(),fitness=c(),run=c())
no_gf_fitness=data.frame()

for ( i in 1:length(file_list)){
  gf_fitness_df=read.table(file_list[i])
  colnames(gf_fitness_df)=c("times","pop","fitness","size","phenotypes","n","GV")
  gf_fitness_df$pop=factor(gf_fitness_df$pop)
  gf_fitness_df$run=i
  no_gf_fitness=rbind(no_gf_fitness,gf_fitness_df)
  gf_fitness_df_mean=gf_fitness_df%>%group_by(times)%>%dplyr::summarise(fitness=mean(fitness))
  gf_fitness_df_mean= gf_fitness_df_mean %>% dplyr::summarise(times=times,fitness=fitness/max(fitness))
  gf_fitness_df_mean$run=i
  no_geneflow_data=rbind(no_geneflow_data,gf_fitness_df_mean)
}

no_gf_fitness$type="Non-gene flow model"
no_gf_fitness_mean=no_gf_fitness %>% group_by(pop,times) %>% dplyr::summarise(fitness=mean(fitness))

no_geneflow_data$group="Non-gene flow model"
no_geneflow_data_mean=no_geneflow_data%>%group_by(times)%>%dplyr::summarise(times=mean(times),fitness=mean(fitness))
no_geneflow_data_mean$group="Non-gene flow model"
scale=max(no_geneflow_data$fitness)

data_all=rbind(geneflow_data,no_geneflow_data)
data_all$group=factor(data_all$group,levels = c("Non-gene flow model","Gene flow model"),labels=c("Non-gene flow model","Gene flow model"))
data_mean=rbind(no_geneflow_data_mean,geneflow_data_mean)
data_mean$group=factor(data_mean$group,levels = c("Non-gene flow model","Gene flow model"),labels=c("Non-gene flow model","Gene flow model"))

##########plot############
fitness=ggplot()+
  geom_line(data=data_all,aes(x=times-burnin,y=fitness/scale,group=interaction(group,run),col=group,alpha=group),size=1,show.legend = F)+
  scale_color_manual(values=c("#e8b89a","#5a9bd7"))+
  scale_alpha_manual(values=c(0.2,0.2))+
  new_scale_color()+
  geom_line(data=data_mean,aes(x=times-burnin,y=fitness/scale,col=group),size=1.5)+
  scale_color_manual(values=c("#e8b89a","#5a9bd7"))+
  theme_classic(base_size = 15) +
  theme(legend.position = c(0.9, 0),       
        legend.justification = c(0.8, 0),        
        legend.text = element_text(size = 18), 
        legend.key.size = unit(2, "lines"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key.height = unit(1.5, "lines"),
        aspect.ratio = 0.5)+ 
  ylab("Relative Fitness")+xlab("Generations")+
  coord_cartesian(xlim=c(-2,200 ),ylim=c(0,1))+theme(legend.position = "none")+ 
  ggtitle(expression("σ"["e"]==1~","~italic("Esd"==1.5)~","~"σ"["t"]==1.5))+theme(plot.title = element_text(size = 15))
fitness


###########################################################
#################Introgression allele rate#################
#The fraction of climate‐adaptive loci shared among species
###########################################################

#Gene flow model
file_list=list.files("D:/slim4/slim_code/slim/1_1.5_1.5/geneflow/",full.names = T,pattern = "^mutation")
geneflow_df=data.frame(times=c(),rate=c(),run=c())
for(i in 1:length(file_list)){
  gf_sites_df=read.table(file_list[i])
  colnames(gf_sites_df)=c("times","standing","introgression","rate","mutation","fitness")
  gf_sites_df$rate=gf_sites_df$introgression/gf_sites_df$standing
  gf_sites_df$run=i
  geneflow_df=rbind(geneflow_df,gf_sites_df)
}

geneflow_df$group="Gene flow model"
geneflow_df_mean=geneflow_df%>%group_by(times)%>%dplyr::summarise(times=mean(times),rate=mean(rate))
geneflow_df_mean$group="Gene flow model"

#Non-gene flow model
file_list=list.files("D:/slim4/slim_code/slim/1_1.5_1.5/no_geneflow/",full.names = T,pattern = "^mutation")
no_geneflow_df=data.frame(times=c(),rate=c(),run=c())

for(i in 1:length(file_list)){
  gf_sites_df=read.table(file_list[i])
  colnames(gf_sites_df)=c("times","standing","introgression","rate","mutation","fitness")
  gf_sites_df$rate=gf_sites_df$introgression/gf_sites_df$standing
  gf_sites_df$run=i
  no_geneflow_df=rbind(no_geneflow_df,gf_sites_df)
}
no_geneflow_df$group="Non-gene flow model"
no_geneflow_df_mean=no_geneflow_df%>%group_by(times)%>%dplyr::summarise(times=mean(times),rate=mean(rate))
no_geneflow_df_mean$group="Non-gene flow model"

data_all=rbind(geneflow_df,no_geneflow_df)
data_all$group=factor(data_all$group,levels = c("Non-gene flow model","Gene flow model"),labels=c("Non-gene flow model","Gene flow model"))
data_mean=rbind(geneflow_df_mean,no_geneflow_df_mean)
data_mean$group=factor(data_mean$group,levels = c("Non-gene flow model","Gene flow model"),labels=c("Non-gene flow model","Gene flow model"))

rate=ggplot()+
  geom_line(data=data_all,aes(x=times-burnin,y=rate,group=interaction(run,group),col=group,alpha=group),alpha=0.2,size=1)+
  scale_color_manual(values=c("#E8B89A","#5A9BD7"))+
  scale_alpha_manual(values=c(0.2,0.2))+
  new_scale_color()+
  geom_line(data=data_mean,aes(x=times-burnin,y=rate,group=group,col=group),size=1.5)+
  scale_color_manual(values=c("#E8B89A","#5A9BD7"))+
  theme_classic(base_size = 15)+
  theme(legend.position = c(0.9, 0),      
        legend.justification = c(0.8, 0),       
        legend.text = element_text(size = 18),  
        legend.key.size = unit(2, "lines"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key.height = unit(1.5, "lines"))+  
  ylab("Introgression Allele Rate")+xlab("Generations")+
  coord_cartesian(xlim=c(-12,200 ),ylim=c(-0.1,0.75))+
  theme(aspect.ratio = 0.5)+
  theme(legend.position = "none")
rate

fitness+rate

