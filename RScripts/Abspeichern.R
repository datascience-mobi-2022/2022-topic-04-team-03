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
library(hexbin)
library(cluster)
library(factoextra)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(RColorBrewer)
library(VennDiagram)
library(GO.db)
library(tinytex)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)


# 2) Read in .CEL files
setwd("~//documents//GitHub//2022-topic-04-team-03//Data//rawData")

data.human=ReadAffy()
data.human@cdfName <- "HGU133Plus2_Hs_ENST"

# 3.2) Normalization

human.vsnrma <- vsnrma(data.human)
human.vsnrma.df = data.frame(exprs(human.vsnrma))
human.vsnrma.df2 = human.vsnrma.df[63:95721,]
rownames(human.vsnrma.df2) = substr(rownames(human.vsnrma.df2), 1,15)

setwd("~//documents//GitHub//2022-topic-04-team-03//Tables")

ensembl.data <- read.csv("ensembl.human.txt")
dim(ensembl.data)
##268341      9

# create variables that contain gene ID, transcript ID, 
#the chromosome name and the gene symbol
ensembl.genes = ensembl.data[,1]
ensembl.transcripts = ensembl.data[,3]
ensembl.chromosome = ensembl.data[,7]
ensembl.symbols = ensembl.data[,9]




# read in TRA data for human
tra.data <- read.table("tra.2017.human.gtex.5x.table.tsv",header=TRUE)
dim(tra.data)
#60131    10


# find the index of rownames of our dataset that contain transcript Ids from the tra dataset
# We extracted the TRA genes out of the three tables(human.vsnrma, tra, ensembl)

j = which(rownames(human.vsnrma.df2) %in% tra.data$ensembl.transcript) 
tra.extracted = rownames(human.vsnrma.df2)[j] #24,783
human.vsnrma.only.tra = human.vsnrma.df2[j,]
#dim(24783,18)

k = which(tra.data$ensembl.transcript %in% tra.extracted)
tra.expressed.in.chips = tra.data[k,] #24,783 
#dim(24783,10)

#150 genes weniger als bei anderen zwei Tabellen
c = which(ensembl.transcripts %in% tra.extracted)
ensembl.only.tra = ensembl.data[c,]
#dim(24623,9)

#reorder the rows of ensembl.only.tra and tra.expressed.in.chips
tra.expressed.in.chips = arrange(tra.expressed.in.chips,ensembl.transcript)
ensembl.only.tra = arrange(ensembl.only.tra,Transcript.stable.ID)

#fusion of two tables! tra.expressed.in.chips and human.vsnrma.only.tra
fusion.tra.expression.tra.table = cbind(human.vsnrma.only.tra,tra.expressed.in.chips$ensembl.gene,tra.expressed.in.chips$tiss.number,tra.expressed.in.chips$tissues,tra.expressed.in.chips$max.tissue)
colnames(fusion.tra.expression.tra.table)[19:22]<-c("ensembl.gene","tissue.number","tissues","max.tissue")


# fusion of fusion.tra.expression.tra.table and ensembl.only.tra because ensembl.new contains less rows than the other two tables
p = which(rownames(fusion.tra.expression.tra.table) %in% ensembl.transcripts)
fusion.tra.expression.tra.table.extracted = fusion.tra.expression.tra.table[p,]
fusion.tra.expression.tra.table.ensembl.table = cbind(fusion.tra.expression.tra.table.extracted,ensembl.only.tra$Gene.name,ensembl.only.tra$Gene.description,ensembl.only.tra$Chromosome.scaffold.name)

#rename colnames
colnames(fusion.tra.expression.tra.table.ensembl.table)[23:25]<-c("Gene.name","Gene.description","Chromosome.scaffold.name")


setwd("~//documents//GitHub//2022-topic-04-team-03//Tables")
write.csv(fusion.tra.expression.tra.table.ensembl.table, file="fusion.tra.expression.tra.table.ensembl.table.csv")

