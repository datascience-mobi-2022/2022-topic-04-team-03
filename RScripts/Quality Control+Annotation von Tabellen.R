# 1) Load libraries

library(affy)
library(vsn)
library(AnnotationDbi)
library(hgu133plus2hsenstcdf)
library(hgu133plus2hsenstprobe)
library(limma)
library(pheatmap)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(gplots) #for heatmap
library(tidyverse) # data manipulation
library(cluster) # clustering algorithms
library(factoextra)
library(gridExtra)
library(Rcpp)
library(tidyverse)
library(affy)
library(vsn)
library(AnnotationDbi)
library(ggplot2)
library(readr)
library(hgu133plus2hsenstcdf)
library(hgu133plus2hsenstprobe)
library(hexbin)


# 2) Read in .CEL files
setwd("/Users/yaxin/Desktop/Bioinformatik/Projekt/Rawdata")

data.human=ReadAffy()
data.human@cdfName <- "HGU133Plus2_Hs_ENST"

setwd("/Users/yaxin/Desktop/Bioinformatik/Projekt")
save.image(file="rawdata_human_18290.rda")


# 3) Quality control
# 3.1) Single chip control

# GSM456643
image(data.human[,1], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456644
image(data.human[,2], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456645
image(data.human[,3], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456646
image(data.human[,4], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456647
image(data.human[,5], col=rainbow(100, start=0, end=0.75)[100:1]) 

# GSM456648
image(data.human[,6], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456649
image(data.human[,7], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456650
image(data.human[,8], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456651
image(data.human[,9], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456652
image(data.human[,10], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456653
image(data.human[,11], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456654
image(data.human[,12], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456655
image(data.human[,13], col=rainbow(100, start=0, end=0.75)[100:1]) 

# GSM456656
image(data.human[,14], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456657
image(data.human[,15], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456658
image(data.human[,16], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456659
image(data.human[,17], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456660
image(data.human[,18], col=rainbow(100, start=0, end=0.75)[100:1])

# 3.2) Normalization

human.vsnrma <- vsnrma(data.human)

setwd("/Users/yaxin/Desktop/Bioinformatik/Projekt/Rawdata")
save.image(file="normalized_human_data_18290.rda")

# 3.3) meanSdPlot

meanSdPlot(human.vsnrma)

setwd("/Users/yaxin/Desktop/Bioinformatik/Projekt/Rawdata/Plots") # was ist das?#
dev.copy2eps(file="meanSdPlot_human_vsnrma_normalized.eps")

# 3.4) Boxplot

# Before normalization:
par(las=2)
mmi=c(1,0.7,1.0477939,0.5366749)
par(mai=mmi)
boxplot(data.human, col= rainbow(35), cex.axis=0.5, main="Gene expression in human embryogenesis data before normalization")

dev.copy2eps(file="boxplot_rawdata_18290.eps")

# After Normalization:
boxplot(exprs(human.vsnrma), col= rainbow(35), cex.axis=0.5, main="Gene expression in human embryogenesis data after normalization")
dev.copy2eps(file="boxplot_normalizeddata_18290.eps")
dev.off() # what does this function do?

# 3.5) Density plot

# Before Normalization:
hist(data.human, col=rainbow(35), main="Density function of log Intensity of human embryogenesis raw data")
dev.copy2pdf(file="Histogram_rawdata_39897.pdf", width = 10, height = 7)

# After Normalization:

eset=exprs(human.vsnrma)

plot(density(eset[,1]), type="n", xlab="log Intensity", ylim=c(0,1), main="Density function of log Intensity of human embryogenesis normalized")
for (i in 1:ncol(eset)) {
  lines(density(eset[,i]), col=rainbow(35)[i])
}

dev.copy2pdf(file="Histogram_NormalizedData_39897.pdf", width = 10, height = 7)

# 3.6) RNA Degradation Plot

par(mfrow = c(2,1))
rnadeg.raw = AffyRNAdeg(data.human)

plotAffyRNAdeg(rnadeg.raw, col=rainbow(35))
title(sub="human embroyogenesis rawdata")

plotAffyRNAdeg(rnadeg.raw, col=rainbow(35), transform= "shift.only")
title(sub="human embryogenesis rawdata - shifted")

dev.copy2pdf(file="RNAdegrad_plot.pdf", width = 12.5, height = 20)
dev.off()

# 3.7) Scatter plot
eset <- exprs(human.vsnrma)

for(i in 1:9){
  plot(eset[,c(i,i+1)], pch=".", cex=2)
  abline(0, 1, col="red")               # 45 degree dividing line
  
  title(main = paste("Scatterplot of probe", 
                     substr(colnames(human.vsnrma)[i], 1, nchar(colnames(human.vsnrma)[i])), "and", 
                     substr(colnames(human.vsnrma)[i+1], 1, nchar(colnames(human.vsnrma)[i+1])), 
                     sep=" ", collapse = NULL))
  
  file.name <- paste("/Users/yaxin/Desktop/Bioinformatik/Projekt/Rawdata", 
                     as.character(substr(colnames(human.vsnrma)[i], 1, nchar(colnames(human.vsnrma)[i]))), "_",
                     as.character(substr(colnames(human.vsnrma)[i+1], 1, nchar(colnames(human.vsnrma)[i+1]))),
                     ".pdf", sep="")
  
  dev.copy2pdf(file = file.name)
  dev.off()
}


# 4) Annotation of genes with Ensemble Biomart
#----------------------------------------------

#annotation table from ensemble biomart with the given features is downloaded:
# ensemble genes 105 -> Mouse genes (GRCm39) -> "Gene.stable.ID", "Gene.stable.ID.version", "Transcript.stable.ID" 
# "Transcript.stable.ID.version", "Chromosome.scaffold.name", "Gen description"       
#"Gene.name".   

setwd("/Users/yaxin/Documents/GitHub/2022-topic-04-team-03/Tables")
a <- read.csv("mart_export.txt")


# exclude Affymetrix control genes which begin with "AFFX"
eset=exprs(human.vsnrma)
eset1 =eset[62:95721,] #the first 62 rows of the data are control genes, there are all together 95721 genes.

ensemble.genes= a[,1]
ensemble.transcripts= a[,2]
ensemble.chromosome= a[,3]
ensemble.description =a[,4]
ensemble.gene.name=a[,5]
#filter each variable

#organize the transcript names in eset1
tr.ID=rownames(eset1)
tr= substr(tr.ID,0,15) #The transkript ID is way to long, this function extracts the first 15 numbers/alphabets of the ID
rownames(eset1)= tr

#read in tra data for human
b <- read.table("tra.2017.human.gtex.5x.table.tsv",header=TRUE)

TRAs= as.character(b[,3])
#60131

TRAs.unique=unique(TRAs)
#19473

# find TRAs  symbols in ensemble data
i=which(ensemble.gene.name %in% TRAs.unique)
TRAs.symbols= ensemble.gene.name[i]
#203610

#extract corresponding TRA transcript id from ensemble data
TRAs.transcripts=ensemble.transcripts[i] #203610
names(TRAs.symbols) = TRAs.transcripts 

# find the index of rownames that contain transcript Ids from TRA-ensemble comparision
i1=which(rownames(eset1) %in% TRAs.transcripts) 
TRA.transcripts2= rownames(eset1)[i1] #77996

TRA.symbols2=TRAs.symbols[TRA.transcripts2] #77996 

#subset of Affymetrix expression data, that contains only transcripts from TRA-Ensemble comparision
eset1.TRA= eset1[i1,]  #77996

# bind the tables
data.TRA.info= cbind(eset1.TRA, TRA.symbols2) #

length(unique(TRA.symbols2)) 
# ---> 13406 out of 19473 TRA genes are present on Affymetrix chips

setwd("/Users/yaxin/Documents/GitHub/2022-topic-04-team-03/Tables")
write.csv(data.TRA.info, file="TRA_Exp_GenName_39897.csv")

########################################
# 5.1) comparision between 1 Cell stage and 2 Cell stage
########################################
M = cbind(eset1[,1]-eset1[,4], 
          eset1[,2]-eset1[,5], 
          eset1[,3]-eset1[,6])

colnames(M) = paste(rep("1 cell stage vs 2 cell stage",3), as.character(1:3), "M", sep=".")

# Create and empty matrix, 
design=as.matrix(rep(1,3))
colnames(design)= "1 cell-2 cell"

#calculate the fit and thus p-values
fit1= lmFit(M,design)
fit1= eBayes(fit1)

pvalue01= sum(p.adjust(fit1$p.value,"BH")< 0.01)
pvalue01 #0

pvalue05= sum(p.adjust(fit1$p.value,"BH")< 0.05)
pvalue05 #20

pvalue1= sum(p.adjust(fit1$p.value,"BH")< 0.1)
pvalue1 #28

pvalue2= sum(p.adjust(fit1$p.value,"BH")< 0.2)
pvalue2 #225

pvalue3= sum(p.adjust(fit1$p.value,"BH")< 0.3)
pvalue3 #811

pvalue4= sum(p.adjust(fit1$p.value,"BH")< 0.4)
pvalue4 #1719

pvalue5= sum(p.adjust(fit1$p.value,"BH")< 0.5)
pvalue5 #5387

pvalue01= sum(fit1$p.value< 0.01)
pvalue01 #0

pvalue05= sum(fit1$p.value< 0.05)
pvalue05 #20

pvalue1= sum(fit1$p.value< 0.1)
pvalue1 #28

pvalue2= sum(fit1$p.value< 0.2)
pvalue2 #225

pvalue3= sum(fit1$p.value< 0.3)
pvalue3 #811

pvalue4= sum(fit1$p.value< 0.4)
pvalue4 #1719

pvalue5= sum(fit1$p.value< 0.5)
pvalue5 #5387
#extracting the differentially expressed transcripts with pvalue=0.5 and lfc=1
top_table1= topTable(fit1, number = pvalue5, lfc=1, p.value = 0.5,sort.by ="logFC")
# -> 18 transcripts differentially expressed

setwd("/Users/nazliaybikeboldemir/Desktop/MoBi Data/ders MASTER/Praktika/Bioinfo/Tables")
write.csv(top_table1, file="topTable1_EggTS01_39897.csv")


#annotation with ensemble data
#-------------------------------
i2= which(ensemble.transcripts %in% rownames(top_table1))
tr.symbols= ensemble.symbols[i2]
trEggTS01 = ensemble.transcripts[i2]
ensem.EggTS01= cbind(trEggTS01,tr.symbols)

#remove duplicates of transcript Ids in ensem.EggTS01 
ensem.EggTS01=as.data.frame(ensem.EggTS01)
ensem.EggTS01=ensem.EggTS01[!duplicated(as.vector(ensem.EggTS01$trEggTS01)), ]
rownames(ensem.EggTS01)= ensem.EggTS01$trEggTS01

# control a random transcript, if the info is correct
#a[a$Transcript.stable.ID =="ENSMUST00000149936",]

#merge two tables
top_table1_anno= merge(top_table1, ensem.EggTS01, by = 'row.names', all = TRUE)

#sort the table regarding the lfc values
top_table1_anno= top_table1_anno[order(abs(top_table1_anno$logFC), decreasing= TRUE),]
rownames(top_table1_anno)= top_table1_anno$Row.names
top_table1_anno$Row.names <- NULL
top_table1_anno$transcripts <- NULL

setwd("/Users/nazliaybikeboldemir/Desktop/MoBi Data/ders MASTER/Praktika/Bioinfo/Tables")
write.csv(top_table1_anno, file="topTable1_anno_EggTS01_39897.csv")

