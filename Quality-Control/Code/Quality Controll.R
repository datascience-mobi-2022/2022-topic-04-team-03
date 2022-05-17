# 1) Load libraries

library(AnnotationDbi)
library(limma)
library(pheatmap)
library(dplyr)
library("Rcpp")
library("tidyverse")
library("affy")
library("vsn")
library("AnnotationDbi")
library("ggplot2")
library("readr")
library("hgu133plus2hsenstcdf")
library("hgu133plus2hsenstprobe")
library("hexbin")


# 2) Read in .CEL files
setwd("C:\\Users\\LKaup\\OneDrive\\Dokumente\\Data Analysis\\2022-topic-04-team-03\\Data\\rawData")

data.human=ReadAffy()
data.human@cdfName <- "HGU133Plus2_Hs_ENST"

setwd("C:\\Linda\\Uni\\Kurse\\SoSe22\\Data Analysis\\Quality Control")
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
image(data.human[,13], col=rainbow(100, start=0, end=0.75)[100:1]) #->very bright?

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

setwd("C:\\Users\\LKaup\\OneDrive\\Dokumente\\Data Analysis\\2022-topic-04-team-03\\Data\\rawData")
save.image(file="normalized_human_data_18290.rda")

# 3.3) meanSdPlot

meanSdPlot(human.vsnrma)

setwd("C:\\Users\\LKaup\\OneDrive\\Dokumente\\Data Analysis\\2022-topic-04-team-03\\Data\\rawData") # was ist das?#
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
  
  file.name <- paste("C:\\Users\\LKaup\\OneDrive\\Dokumente\\Data Analysis\\2022-topic-04-team-03\\Data\\rawData", 
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

setwd("C:/Users/LKaup/OneDrive/Dokumente/Data Analysis/2022-topic-04-team-03/Tables")
table <- read.csv("fusion.fusion.ensembl.csv")
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
#[1] 95721    18

# exclude Affymetrix control genes which begin with "AFFX"
human.vsnrma.df2 = human.vsnrma.df[63:9571,]

dim(human.vsnrma.df2)
# [1] 9509   18

#remove ".xx_at" from the rownames
rownames(human.vsnrma.df2) = substr(rownames(human.vsnrma.df2), 1,15)

# read in TRA data for human
tra.data <- read.table("tra.2017.human.gtex.5x.table.tsv",header=TRUE)

# Exclude the following:
# haben wir weggelassen
# TRAs= as.character(tra.data[,3])

# Extract unique TRA-Symbols
#tra.symbols = tra.data[,3]

#tra.unique = unique(tra.symbols)
# -> 18473 unique tras

# find TRAs  symbols in ensembl data -> evtl nach transcript filtern nicht nach symbol
#i = which(ensembl.symbols %in% tra.unique)
#tra.symbols = ensembl.symbols[i]
# -> 204,514
# there are still symbols without name in tra.symbols! filter before?



# extract corresponding TRA transcript id from ensemble data
#tra.transcripts.TE = ensembl.transcripts[i] 
# -> 204,514
#names(tra.symbols) = tra.transcripts.TE
u = which(ensembl.transcripts %in% tra.data)

# find the index of rownames of our dataset that contain transcript Ids from TRA-ensembl comparison
## hier statt tra.transcripts.TE nach tra.data[,1] (=Transcript namen) filtern
j = which(rownames(human.vsnrma.df2) %in% tra.transcripts.TE) 
tra.transcripts.TEH = rownames(human.vsnrma.df2)[j] 
# -> 6427

# assign symbols to extracted tra transcript IDs
tra.transcripts.symbols.TEH = tra.symbols[tra.transcripts.TEH]

#Take dataset and extract tras
human.vsnrma.tra = human.vsnrma.df2[j,]
# -> 7,427

## hier weitermachen
# bind the tables
data.TRA.info = cbind(human.vsnrma.tra, tra.transcripts.symbols.TEH)


## wahrscheinlich auch unnötig?:
# Extract tras from our dataset in tra-dataset
k = which(tra.data[,1] %in% rownames(human.vsnrma.tra))
tra.extracted = tra.data[k,]

tra.extracted = arrange(tra.extracted, ensembl.transcript)


## hier noch unsere Daten einfügen:
length(unique(TRA.symbols2)) 
# ---> 3785 out of 4154 TRA genes are present on Affymetrix chips

setwd("/Users/nazliaybikeboldemir/Desktop/MoBi Data/ders MASTER/Praktika/Bioinfo/Tables")
write.csv(data.TRA.info, file="TRA_Exp_GenName_39897.csv")



colnames(human.vsnrma.df2)

#Info from GEO datenbank
#GSM456643	human embryo at 1 cell stage, biological rep 1
#GSM456644	human embryo at 1 cell stage, biological rep 2
#GSM456645	human embryo at 1 cell stage, biological rep 3
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

########################################
# 5.1) comparision between embryo 1 cell stage and embryo 2 cell stage
########################################
M = cbind(human.vsnrma.df[,1]-human.vsnrma.df[,4], 
          human.vsnrma.df[,2]-human.vsnrma.df[,5], 
          human.vsnrma.df[,3]-human.vsnrma.df[,6])

colnames(M) = paste(rep("1S vs 2S",3), as.character(1:3), "M", sep=".")

# Create and empty matrix, 
design=as.matrix(rep(1,3))
colnames(design)= "1S-2S"

#calculate the fit and thus p-values
fit1= lmFit(M,design)
fit1= eBayes(fit1)

pvalue01= sum(p.adjust(fit1$p.value,"BH")< 0.01)
pvalue01 #0

pvalue05= sum(p.adjust(fit1$p.value,"BH")< 0.05)
pvalue05 #20


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





