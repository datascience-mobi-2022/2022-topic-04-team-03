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

# 2) Read in .CEL files
setwd("/Users/Clara/Documents/Studium/Bioinformatik/Projekt/rawdata/GSE18290_RAW")
data.human=ReadAffy()
data.human@cdfName <- "HGU133Plus2_Hs_ENST"

# 3) Normalization
human.vsnrma <- vsnrma(data.human)

# 4) Annotation
setwd("/Users/Clara/Documents/Studium/Bioinformatik/Projekt/rawdata")
ensembl.data <- read.csv("ensembl.human.txt")

ensembl.genes = ensembl.data[,1]
ensembl.transcripts = ensembl.data[,3]
ensembl.chromosome = ensembl.data[,7]
ensembl.symbols = ensembl.data[,9]

human.vsnrma.df = data.frame(exprs(human.vsnrma))

human.vsnrma.df2 = human.vsnrma.df[63:95721,]

rownames(human.vsnrma.df2) = substr(rownames(human.vsnrma.df2), 1,15)

tra.data <- read.table("tra.2017.human.gtex.5x.table.tsv",header=TRUE)

j = which(rownames(human.vsnrma.df2) %in% tra.data$ensembl.transcript) 
tra.extracted = rownames(human.vsnrma.df2)[j] #24,783
human.vsnrma.new = human.vsnrma.df2[j,]

k = which(tra.data$ensembl.transcript %in% tra.extracted)
tra.new = tra.data[k,] #24,783 

c = which(ensembl.transcripts %in% tra.extracted)
ensembl.new = ensembl.data[c,]


tra.new = arrange(tra.new,ensembl.transcript)
ensembl.new = arrange(ensembl.new,Transcript.stable.ID)

fusion.tra.human = cbind(human.vsnrma.new,tra.new$ensembl.gene,tra.new$tiss.number,tra.new$tissues,tra.new$max.tissue)

p = which(rownames(fusion.tra.human) %in% ensembl.transcripts)
fusion.tra.human.extracted = fusion.tra.human[p,]
fusion.fusion.ensembl = cbind(fusion.tra.human.extracted,ensembl.new$Gene.description,ensembl.new$Chromosome.scaffold.name)


# 5.3) Limma analysis
# Define the different stages of the chips in a vector
stage = c(1,1,1,2,2,2,4,4,4,8,8,8,"morula","morula","morula", "blastocyst", "blastocyst", "blastocyst")

# Create a design matrix with the stages
design = model.matrix(~0 + stage)

# Create a contrast matrix which compares all the stages with one another 
cm = makeContrasts(stage1.2 = stage1-stage2, stage1.4 = stage1-stage4, stage1.8 = stage1-stage8, stage1.morula = stage1-stagemorula, stage1.blastocyst = stage1-stageblastocyst,
                   stage2.4 = stage2-stage4, stage2.8 = stage2-stage8, stage2.morula = stage2-stagemorula,  stage2.blastocyst = stage2-stageblastocyst,
                   stage4.8 = stage4-stage8, stage4.morula = stage4-stagemorula,  stage4.blastocyst = stage4-stageblastocyst,
                   stage8.morula = stage8-stagemorula, stage8.blastocyst = stage8-stageblastocyst,
                   stagemorula.blastocyst = stagemorula-stageblastocyst, 
                   levels = design)

# Fit the coefficients of the model
fit = lmFit(human.vsnrma.df2, design)
fit2 = contrasts.fit(fit, contrasts = cm)

# Calculate the t-statistics
fit2 = eBayes(fit2)

# Number of differentially expressed genes (p-value = 0.05)
results = decideTests(fit2, p.value = 0.05)
summary(results)

# Number of differentially expressed genes (p-value = 0.01)
results2 = decideTests(fit2, p.value = 0.01)
summary(results2)

# Create fits for every contrast and then create tables with the results of Limma analysis
# Stage 1-2
fit.1.2 = contrasts.fit(fit, contrasts = cm[,1])
fit.1.2 = eBayes(fit.1.2)

pvalue01 = sum(p.adjust(fit.1.2$p.value, "BH") < 0.01)
# 0
pvalue05 = sum(p.adjust(fit.1.2$p.value, "BH") < 0.05)
# 0

# Create fits for every contrast and then create tables with the results of limma analysis

setwd("~//documents//GitHub//2022-topic-04-team-03//Tables")


for(i in 2:15) {
  fit.1 = contrasts.fit(fit, contrasts = cm[,i])
  fit.1 = eBayes(fit.1)
  
  pvalue01 = sum(p.adjust(fit.1$p.value, "BH") < 0.01)
  pvalue05 = sum(p.adjust(fit.1$p.value, "BH") < 0.05)
  
  if (i==2){
    limma.table.1.4 = topTable(fit.1, number = pvalue01)
  }
  
  if (i==3){
    limma.table.1.8 = topTable(fit.1, number = pvalue01)
  }
  
  if (i==4){
    limma.table.1.m = topTable(fit.1, number = pvalue01)
  }
  
  if (i==5){
    limma.table.1.b = topTable(fit.1, number = pvalue01)
  }
  
  if (i==6){
    limma.table.2.4 = topTable(fit.1, number = pvalue01)
  }
  
  if (i==7){
    limma.table.2.8 = topTable(fit.1, number = pvalue01)
  }
  
  if (i==8){
    limma.table.2.m = topTable(fit.1, number = pvalue01)
  }
  
  if (i==9){
    limma.table.2.b = topTable(fit.1, number = pvalue01)
  }
  
  if (i==10){
    limma.table.4.8 = topTable(fit.1, number = pvalue01)
  }
  
  if (i==11){
    limma.table.4.m = topTable(fit.1, number = pvalue01)
  }
  
  if (i==12){
    limma.table.4.b = topTable(fit.1, number = pvalue01)
  }
  
  if (i==13){
    limma.table.8.m = topTable(fit.1, number = pvalue01)
  }
  
  if (i==14){
    limma.table.8.b = topTable(fit.1, number = pvalue01)
  }
  
  if (i==15){
    limma.table.m.b = topTable(fit.1, number = pvalue01)
  }
  
  #file.name = paste("limma.table",colnames(cm)[i], sep=" ")
  #toString(filename) = limma.table.1.i                 
  #(limma.table.1.i,file = file.name)
}

# Annotate topTables
# Define a function for annotation 
# Parameter x is the limma table

limma.annotation = function(x){
  a = which(ensembl.transcripts %in% rownames(x))
  ensembl.new = ensembl.data[a,]
  
  b = which(rownames(x) %in% ensembl.new$Transcript.stable.ID)
  x.new = x[b,]
  
  ensembl.new = arrange(ensembl.new,Transcript.stable.ID)
  x.new = arrange(x.new, rownames(x.new))
  
  annotated.x = cbind(x.new, ensembl.new$Gene.name, ensembl.new$Gene.description,ensembl.new$Chromosome.scaffold.name)
  annotated.x = arrange(annotated.x, adj.P.Val)
  return(annotated.x)
}

annotated.limma.1.4 = limma.annotation(limma.table.1.4)
annotated.limma.1.8 = limma.annotation(limma.table.1.8)
annotated.limma.1.m = limma.annotation(limma.table.1.m)
annotated.limma.1.b = limma.annotation(limma.table.1.b)
# annotated.limma.2.4 = limma.annotation(limma.table.2.4)
# Error
annotated.limma.2.8 = limma.annotation(limma.table.2.8)
# Error
annotated.limma.2.m = limma.annotation(limma.table.2.m)
annotated.limma.2.b = limma.annotation(limma.table.2.b)
annotated.limma.4.8 = limma.annotation(limma.table.4.8)
annotated.limma.4.m = limma.annotation(limma.table.4.m)
annotated.limma.4.b = limma.annotation(limma.table.4.b)
annotated.limma.8.m = limma.annotation(limma.table.8.m)
annotated.limma.8.b = limma.annotation(limma.table.8.b)
annotated.limma.m.b = limma.annotation(limma.table.m.b)