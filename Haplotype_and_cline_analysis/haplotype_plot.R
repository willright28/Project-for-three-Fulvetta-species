library(dplyr)
library(LEA)
library(pheatmap)

#######
output=vcf2geno("/mnt/disk1/zhangshuai/hkqm/af_shifts/hd_intro.vcf",output.file ="/mnt/disk1/zhangshuai/hkqm/af_shifts/hd_intro.geno" ) #only include adaptive introgressed SNPs in the VCF file  
lines <- readLines("/mnt/disk1/zhangshuai/hkqm/af_shifts/hd_intro.geno")

char_list <- strsplit(lines, "")
char_matrix <- do.call(rbind, char_list)        
num_matrix <- matrix(as.integer(unlist(char_list)), nrow = length(char_list),  byrow = TRUE)
id=read.csv("/mnt/disk1/zhangshuai/hkqm/wza/id.txt",head=F)

colnames(num_matrix)=id$V1
num_matrix=t(num_matrix)
# View(num_matrix)
num_matrix=data.frame(num_matrix)
num_matrix$id=id$V1

#载入分组信息
info=read.csv("/mnt/disk1/zhangshuai/hkqm/gghybrid/pop.0.9.csv",header=T)
info=info[,c("sample","pop_mix","site_number","sample_number","long","lat")]
num_matrix=left_join(num_matrix,info,by = c("id" = "sample"))
num_matrix$pop_mix=factor(num_matrix$pop_mix,levels=c("Central","C_mix","E_mix","Eastern","Western"))
num_matrix=arrange(num_matrix,pop_mix)
info=num_matrix[,(ncol(num_matrix)-5):ncol(num_matrix)]
num_matrix=num_matrix[,1:(ncol(num_matrix)-6)]
dim(num_matrix)


#确定纯合状态
central=colSums(num_matrix[1:22,]%>%as.matrix())
ancest=central>22

num_matrix[,ancest][num_matrix[,ancest] == 0] = 3
num_matrix[,ancest][num_matrix[,ancest] == 2] = 0
num_matrix[,ancest][num_matrix[,ancest] == 3] = 2
east=colSums(num_matrix[41:53,]%>%as.matrix())
central=colSums(num_matrix[1:22,]%>%as.matrix())
west=colSums(num_matrix[54:71,]%>%as.matrix())
mix=colSums(num_matrix[23:34,]%>%as.matrix())
a=(central==0 & east==26   )%>%as.vector()
table(a)
num_matrix[c(41:53),a]
colnames(num_matrix[c(41:53),a])%>%unique()


###绘制gene plot
snp_info=read.table("/mnt/disk1/zhangshuai/hkqm/af_shifts/hd_intro.vcfsnp")[,3]
snp_info=snp_info[a]
snp_info

pdf("/mnt/disk1/zhangshuai/hkqm/hap/filet.hueti.2.pdf",width = 6, height = 3)
pheatmap(num_matrix[1:53,a],fontsize=2, cluster_rows = F, cluster_cols = F,   gaps_row = c(22,40),gaps_col =c(3,65,145,160,372,476,500,505,508,600,525,530,540,545,585,632,698,699,702),
          color = c("#bab7d7",  "#ff9292","#8accc0")   )  
dev.off()
