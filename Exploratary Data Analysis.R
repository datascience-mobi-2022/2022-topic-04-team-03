#please copy paste to our main project

#Exploratary Data Analysis
##Boxplot for distribution of all TRAs that are expressed in our chip

fusion.fusion.ensembl.data.frame=as.data.frame(fusion.fusion.ensembl)

###generate a list containing all tissue name with number of tissue TRAs, which is detected in our chips
tissue=c("Brain","Esophagus","Heart","Liver","Cervix","Muscle","Ovary","Colon","Breast","Kidney","Pituitary","Testis","Whole Blood","Cells","Spleen","Lung","Stomach","Pancreas","Small Intestine","Skin","Artery","Adrenal","Gland","Pituitary","Nerve","Minor","Salivary","Bladder","Adipose","Thyroid","Prostate","Vagina")
tissue.number=c()

for (i in 1:32) {
  tissue.number[[i]]=nrow(fusion.fusion.ensembl.data.frame %>% filter(grepl(tissue[i],V22)))
}
print(tissue.number)
tissue.distribution=cbind(tissue,tissue.number)
tissue.distribution=as.data.frame(tissue.distribution)

###arranging the list according to the descending TRAs number
tissue.distribution.arranged=arrange(tissue.distribution,desc(as.numeric(tissue.number)))
library(RColorBrewer)
coul <- brewer.pal(5, "Set2") 

barplot(height=as.numeric(tissue.distribution.arranged$tissue.number[1:10]),names.arg=tissue.distribution.arranged$tissue[1:10],cex.names=0.8,col=coul)

setwd("/Users/yaxin/Desktop/Bioinformatik/Projekt")

dev.copy2pdf(file="tissue.distribution.arranged.barplot.pdf", width = 10, height = 7)

###now maybe a boxplot and a heatmap could help to show our data better. 
###Before doing that, I thought it might be helpful to calculate the mean for every 3 replicates
###(I didn't pay attention to the deleted chip, we can adapt the code later
###I defined "human.vsnrma.new" which also does not contain the 150 genes that the ensembl table does not contain.
human.vsnrma.new.mean=c()
for (i in 1:6)
  human.vsnrma.new.mean[[i]]=((human.vsnrma.new[,i]+human.vsnrma.new[,(i+1)]+human.vsnrma.new[,(i+2)])/3)
boxplot(human.vsnrma.new.mean)

human.vsnrma.new.mean=matrix(unlist(human.vsnrma.new.mean), ncol =6, nrow =24783)
rownames(human.vsnrma.new.mean)=rownames(human.vsnrma.new)
pheatmap(human.vsnrma.new.mean,cluster_cols=FALSE)