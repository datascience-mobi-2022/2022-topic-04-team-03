

# 1) Load libraries

library(affy)
library(vsn)
library(AnnotationDbi)
library(hgu133plus2hsenstcdf)
library(hgu133plus2hsenstprobe)
library(limma)
library(pheatmap)
#library(dplyr)
library(tidyverse)
library(hexbin)


# 2) Read in .CEL files
setwd("/Users/Clara/Documents/Studium/Bioinformatik/Projekt/rawdata/GSE18290_RAW")

data.human=ReadAffy()
data.human@cdfName <- "HGU133Plus2_Hs_ENST"

setwd("/Users/Clara/Documents/Studium/Bioinformatik/Projekt/sessions/rda")
save.image(file="rawdata_human_18290.rda")


# 3) Quality control
# 3.1) Single chip control

# GSM456643
image(data.human[,1], col=rainbow(100, start=0, end=0.75)[100:1]) 

# GSM456644
image(data.human[,2], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456645
image(data.human[,3], col=rainbow(100, start=0, end=0.75)[100:1]) #-> a little smudge on the bottom edge

# GSM456646
image(data.human[,4], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456647
image(data.human[,5], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456648
image(data.human[,6], col=rainbow(100, start=0, end=0.75)[100:1])

# GSM456649
image(data.human[,7], col=rainbow(100, start=0, end=0.75)[100:1]) #->very bright?

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
image(data.human[,17], col=rainbow(100, start=0, end=0.75)[100:1]) #-> little dots and little line at the bottom edge

# GSM456660
image(data.human[,18], col=rainbow(100, start=0, end=0.75)[100:1])

#exclude the broken chip (GSM456645) from data
data.human <- data.human[,-c(3)]

# 3.2) Normalization

human.vsnrma <- vsnrma(data.human)

setwd("/Users/Clara/Documents/Studium/Bioinformatik/Projekt/sessions/rda")
save.image(file="normalized_human_data_18290.rda")


# 3.3) meanSdPlot

meanSdPlot(human.vsnrma)

setwd("/Users/Clara/Documents/Studium/Bioinformatik/Projekt/plots")
dev.copy2eps(file="meanSdPlot_human_vsnrma_normalized.eps")


# 3.4) Boxplot

# Before normalization:
par(las=2)
mmi=c(1,0.7,1.0477939,0.5366749)
par(mai=mmi)
boxplot(data.human, col= rainbow(35), cex.axis=0.5, main="Gene expression in human embroyogenesis data before normalization")

dev.copy2eps(file="boxplot_rawdata_18290.eps")

# After Normalization:
boxplot(exprs(human.vsnrma), col= rainbow(35), cex.axis=0.5, main="Gene expression in human embroyogenesis data after normalization")
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


# 3.7) Scatter Plot

expression.data <- exprs(human.vsnrma)

for(i in 1:9){
  plot(expression.data[,c(i,i+1)], pch=".", cex=2)
  abline(0, 1, col="red")               # 45 degree dividing line
  
  title(main = paste("Scatterplot of probe", 
                     substr(colnames(human.vsnrma)[i], 1, nchar(colnames(human.vsnrma)[i])), "and", 
                     substr(colnames(human.vsnrma)[i+1], 1, nchar(colnames(human.vsnrma)[i+1])), 
                     sep=" ", collapse = NULL))
  
  file.name <- paste("/Users/Clara/Documents/Studium/Bioinformatik/Projekt/rawdata", 
                     as.character(substr(colnames(human.vsnrma)[i], 1, nchar(colnames(human.vsnrma)[i]))), "_",
                     as.character(substr(colnames(human.vsnrma)[i+1], 1, nchar(colnames(human.vsnrma)[i+1]))),
                     ".pdf", sep="")
  
  dev.copy2pdf(file = file.name)
  dev.off()
}



# 4) Data Analysis

# read in ensembl table with following attributes: "Gene stable ID", 
#"Gene stable ID version", "Transcript stable ID", "Transscript stable ID version", 
#"Gene.name", "Transcript name", "Chromosome.scaffold.name", "Gene.description", "HGNC.symbol"

setwd("/Users/Clara/Documents/Studium/Bioinformatik/Projekt/rawdata")

ensembl.data <- read.csv("ensembl.human.txt")

# create variables that contain gene ID, transcript ID, 
#the chromosome name and the gene symbol
ensembl.genes = ensembl.data[,1]
ensembl.transcripts = ensembl.data[,3]
ensembl.chromosome = ensembl.data[,7]
ensembl.symbols = ensembl.data[,9]

# Create a data frame out of the expression data from the normalized data
human.vsnrma.df = data.frame(exprs(human.vsnrma))

# Check dimensions
dim(human.vsnrma.df)
#[1] 95,721    17

# exclude Affymetrix control genes which begin with "AFFX"
human.vsnrma.df2 = human.vsnrma.df[63:95721,]

dim(human.vsnrma.df2)
# [1] 95,659    17

#remove ".xx_at" from the rownames
rownames(human.vsnrma.df2) = substr(rownames(human.vsnrma.df2), 1,15)

# read in TRA data for human
tra.data <- read.table("tra.2017.human.gtex.5x.table.tsv",header=TRUE)


# find the index of rownames of our dataset that contain transcript Ids from the tra dataset
#We extracted the TRA genes out of the three table(human.vsnrma, tra, ensembl)

j = which(rownames(human.vsnrma.df2) %in% tra.data$ensembl.transcript) 
tra.extracted = rownames(human.vsnrma.df2)[j] #24,783
human.vsnrma.new = human.vsnrma.df2[j,]

k = which(tra.data$ensembl.transcript %in% tra.extracted)
tra.new = tra.data[k,] #24,183 

#150 genes weniger als bei anderen zwei Tabellen
c = which(ensembl.transcripts %in% tra.extracted)
ensembl.new = ensembl.data[c,]


#reorder the rows of ensembl.new and tra.new
tra.new = arrange(tra.new,ensembl.transcript)
ensembl.new = arrange(ensembl.new,Transcript.stable.ID)

#fusion of two tables! tra.new and human.vsnrma.new
fusion.tra.human = cbind(human.vsnrma.new,tra.new$ensembl.gene,tra.new$tiss.number,tra.new$tissues,tra.new$max.tissue)
#colnames(fusion.tra.human[,18:21])=c("gene.name","tissue.number","tissue","max.tissue")
#rename(fusion.tra.human, V19 = gene.name, V20 = tissue.number, V21=tissue, V22=max.tissue)#funktioniert nicht

# fusion of fusion.tra.human and ensembl.new because ensembl.new contains less rows than the other two tables
p = which(rownames(fusion.tra.human) %in% ensembl.transcripts)
fusion.tra.human.extracted = fusion.tra.human[p,]
fusion.fusion.ensembl = cbind(fusion.tra.human.extracted,ensembl.new$Gene.description,ensembl.new$Chromosome.scaffold.name)

setwd("/Users/Clara/Documents/Studium/Bioinformatik/Projekt/Tables")
write.csv(fusion.fusion.ensembl, file="fusion.fusion.ensembl.csv")


# 4.1) Exploratory data analysis

#GSM456643 human embryo at 1 cell stage, biological rep 1
#GSM456644 human embryo at 1 cell stage, biological rep 2

#GSM456646	human embryo at 2 cell stage, biological rep 1
#GSM456647	human embryo at 2 cell stage, biological rep 2
#GSM456648	human embryo at 2 cell stage, biological rep 3
#GSM456649	human embryo at 4 cell stage, biological rep 1
#GSM456650	human embryo at 4 cell stage, biological rep 2
#GSM456651	human embryo at 4 cell stage, biological rep 3
#GSM456652	human embryo at 8 cell stage, biological rep 1
#GSM456653	human embryo at 8 cell stage, biological rep 2
#GSM456654	human embryo at 8 cell stage, biological rep 3
#GSM456655	human embryo at morula stage, biological rep 1
#GSM456656	human embryo at morula stage, biological rep 2
#GSM456657	human embryo at morula stage, biological rep 3
#GSM456658	human blastocyst, biological rep 1
#GSM456659	human blastocyst, biological rep 2
#GSM456660	human blastocyst, biological rep 3





###########################################
# 4.2) comparison between 1 Cell stage and 2 Cell stage
########################################
M = cbind(human.vsnrma.new[,1]-human.vsnrma.new[,3], 
          human.vsnrma.new[,2]-human.vsnrma.new[,4])

colnames(M) = paste(rep("1 cell stage vs 2 cell stage",3), as.character(1:3), "M", sep=".")

# Create and empty matrix, 
design = as.matrix(rep(1,3))
colnames(design) = "1 cell-2 cell"

#calculate the fit and thus p-values
fit1= lmFit(M,design)
fit2= eBayes(fit1)

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







###################################









