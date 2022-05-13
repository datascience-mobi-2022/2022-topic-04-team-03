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
# ensemble genes 105 -> Human genes (GRCm39) -> "Gene.stable.ID", "Gene.stable.ID.version", "Transcript.stable.ID" 
# "Transcript.stable.ID.version", "Chromosome.scaffold.name", "Gen description"       
#"Gene.name".   

setwd("/Users/yaxin/Documents/GitHub/2022-topic-04-team-03/Tables")
ensembl.data <- read.csv("ensembl.human.txt")


# exclude Affymetrix control genes which begin with "AFFX"
human.vsnrma.df=exprs(human.vsnrma)
human.vsnrma.df2 =human.vsnrma.df[63:95721,] #the first 62 rows of the data are control genes, there are all together 95721 genes.

ensembl.genes= ensembl.data[,1]
ensembl.transcripts= ensembl.data[,3]
ensembl.chromosome= ensembl.data[,7]
ensembl.description =ensembl.data[,8]
ensembl.symbols=ensembl.data[,9]
#filter each variable

#organize the transcript names in eset1
rownames(human.vsnrma.df2)=substr(rownames(human.vsnrma.df2),1,15)

#read in tra data for human
tra.data <- read.table("tra.2017.human.gtex.5x.table.tsv",header=TRUE)

#TRAs= as.character(b[,3])
#60131

#TRAs.unique=unique(TRAs)
#19473

# find TRAs  symbols in ensemble data
#i=which(ensembl.symbols %in% TRAs.unique)
#TRAs.symbols= ensemble.gene.name[i]
#203610

#extract corresponding TRA transcript id from ensemble data
#TRAs.transcripts=ensemble.transcripts[i] #203610
#names(TRAs.symbols) = TRAs.transcripts 

# find the index of rownames that contain transcript Ids from TRA-ensemble comparision
#We extracted the TRA genes out of the three table(human.vsnrma, tra, ensembl)

j=which(rownames(human.vsnrma.df2) %in% tra.data$ensembl.transcript) 
tra.extracted = rownames(human.vsnrma.df2)[j] #24,783
engength
k=which(tra.data$ensembl.transcript %in% tra.extracted)
tra.new=tra.data[k,] #77996 

#150 genes weniger als bei anderen zwei Tabellen
c=which(ensembl.transcripts %in% tra.extracted)
ensembl.new=ensembl.data[c,]


human.vsnrma.new=human.vsnrma.df2[j,]


#reorder the rows of ensembl.new and tra.new
tra.new=arrange(tra.new,ensembl.transcript)
ensembl.new=arrange(ensembl.new,Transcript.stable.ID)

#fusion of two tables!
fusion.tra.human=cbind(human.vsnrma.new,tra.new$ensembl.gene,tra.new$tiss.number,tra.new$tissues,tra.new$max.tissue)
#colnames(fusion.tra.human[,19:22])=c("gene.name","tissue.number","tissue","max.tissue")
#rename(fusion.tra.human, V19 = gene.name, V20 = tissue.number, V21=tissue, V22=max.tissue)#funktioniert nicht


p=which(rownames(fusion.tra.human) %in% ensembl.transcripts)
fusion.tra.human.extracted=fusion.tra.human[p,]
fusion.fusion.ensembl=cbind(fusion.tra.human.extracted,ensembl.new$Gene.description,ensembl.new$Chromosome.scaffold.name)

setwd("/Users/yaxin/Documents/GitHub/2022-topic-04-team-03/Tables")
write.csv(data.TRA.info, file="fusion.fusion.ensembl.csv")
write.csv(data.TRA.info, file="fusion.tra.human.csv")


#ensembl.new$Gene.description,ensembl.new$Chromosome.scaffold.name,of Affymetrix expression data, that conmrtains only transcripts from TRA-Ensemble comparision
#eset1.TRA= eset1[i1,]  #77996

# bind the tables
#data.TRA.info= cbind(eset1.TRA, TRA.symbols2) #

#length(unique(TRA.symbols2)) 
# ---> 13406 out of 19473 TRA genes are present on Affymetrix chips

setwd("/Users/yaxin/Documents/GitHub/2022-topic-04-team-03/Tables")
write.csv(data.TRA.info, file="TRA_Exp_GenName_39897.csv")

########################################
# 5.1) comparision between 1 Cell stage and 2 Cell stage
########################################
M = cbind(human.vsnrma.new[,1]-human.vsnrma.new[,4], 
          human.vsnrma.new[,2]-human.vsnrma.new[,5], 
          human.vsnrma.new[,3]-human.vsnrma.new[,6])

colnames(M) = paste(rep("1 cell stage vs 2 cell stage",3), as.character(1:3), "M", sep=".")

# Create and empty matrix, 
design=as.matrix(rep(1,3))
colnames(design)= "1 cell-2 cell"

#calculate the fit and thus p-values
fit1= lmFit(M,design)
fit2= eBayes(fit1)


M = cbind(human.vsnrma.new[,4]-human.vsnrma.new[,7], 
          human.vsnrma.new[,5]-human.vsnrma.new[,8], 
          human.vsnrma.new[,6]-human.vsnrma.new[,9])

colnames(M) = paste(rep("1 cell stage vs 2 cell stage",3), as.character(1:3), "M", sep=".")

# Create and empty matrix, 
design=as.matrix(rep(1,3))
colnames(design)= "1 cell-2 cell"

#calculate the fit and thus p-values
fit1= lmFit(M,design)
fit2= eBayes(fit1)


M = cbind(human.vsnrma.new[,7]-human.vsnrma.new[,10], 
          human.vsnrma.new[,8]-human.vsnrma.new[,11], 
          human.vsnrma.new[,9]-human.vsnrma.new[,12])

colnames(M) = paste(rep("1 cell stage vs 2 cell stage",3), as.character(1:3), "M", sep=".")

# Create and empty matrix, 
design=as.matrix(rep(1,3))
colnames(design)= "1 cell-2 cell"

#calculate the fit and thus p-values
fit1= lmFit(M,design)
fit2= eBayes(fit1)

M = cbind(human.vsnrma.new[,10]-human.vsnrma.new[,13], 
          human.vsnrma.new[,11]-human.vsnrma.new[,14], 
          human.vsnrma.new[,12]-human.vsnrma.new[,15])

colnames(M) = paste(rep("1 cell stage vs 2 cell stage",3), as.character(1:3), "M", sep=".")

# Create and empty matrix, 
design=as.matrix(rep(1,3))
colnames(design)= "1 cell-2 cell"

#calculate the fit and thus p-values
fit1= lmFit(M,design)
fit2= eBayes(fit1)


M = cbind(human.vsnrma.new[,13]-human.vsnrma.new[,16], 
          human.vsnrma.new[,14]-human.vsnrma.new[,17], 
          human.vsnrma.new[,15]-human.vsnrma.new[,18])

colnames(M) = paste(rep("1 cell stage vs 2 cell stage",3), as.character(1:3), "M", sep=".")

# Create and empty matrix, 
design=as.matrix(rep(1,3))
colnames(design)= "1 cell-2 cell"

#calculate the fit and thus p-values
fit1= lmFit(M,design)
fit2= eBayes(fit1)


i=which(fit2$p.value<0.05)
fit3=fit2[i,]
#5230


