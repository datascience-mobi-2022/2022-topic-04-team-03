

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


# 2) Read in .CEL files
setwd("~//documents//GitHub//2022-topic-04-team-03//Data//rawData")

data.human=ReadAffy()
data.human@cdfName <- "HGU133Plus2_Hs_ENST"

#setwd("/Users/Clara/Documents/Studium/Bioinformatik/Projekt/sessions/rda")#
#save.image(file="rawdata_human_18290.rda")#


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


# 3.2) Normalization

human.vsnrma <- vsnrma(data.human)

#setwd("/Users/Clara/Documents/Studium/Bioinformatik/Projekt/sessions/rda")#
#save.image(file="normalized_human_data_18290.rda")#


# 3.3) meanSdPlot

meanSdPlot(human.vsnrma)

setwd("~//documents//GitHub//2022-topic-04-team-03//Plots")
dev.copy2eps(file="meanSdPlot_human_vsnrma_normalized.eps")


# 3.4) Boxplot

names = c("1 cell stage, rep 1", "1 cell stage, rep 2", "1 cell stage, rep 3",
          "2 cell stage, rep 1", "2 cell stage, rep 2", "2 cell stage, rep 3",
          "4 cell stage, rep 1", "4 cell stage, rep 2", "4 cell stage, rep 3",
          "8 cell stage, rep 1", "8 cell stage, rep 2", "8 cell stage, rep 3",
          "morula stage, rep 1", "morula stage, rep 2", "morula stage, rep 3",
          "blastocyst stage, rep 1", "blastocyst stage, rep 2", "blastocyst stage, rep 3")


# Before normalization:
par(las=2)
mmi=c(1,0.7,1.0477939,0.5366749)
par(mai=mmi)
boxplot(data.human, col= rainbow(35), cex.axis=0.5, main="Gene expression in human embroyogenesis data before normalization", xlim = names)

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
#dev.off()


# 3.7) Scatter Plot

setwd("~//documents//GitHub//2022-topic-04-team-03//Plots")
expression.data <- exprs(human.vsnrma)

for(i in 1:17){
  plot(expression.data[,c(i,i+1)], pch=".", cex=2)
  abline(0, 1, col="red")               # 45 degree dividing line
  
  title(main = paste("Scatterplot of probe", 
                     substr(colnames(human.vsnrma)[i], 1, nchar(colnames(human.vsnrma)[i])), "and", 
                     substr(colnames(human.vsnrma)[i+1], 1, nchar(colnames(human.vsnrma)[i+1])), 
                     sep=" ", collapse = NULL))
  
  file.name <- paste("~//documents//GitHub//2022-topic-04-team-03//Plots", 
                     as.character(substr(colnames(human.vsnrma)[i], 1, nchar(colnames(human.vsnrma)[i]))), "_",
                     as.character(substr(colnames(human.vsnrma)[i+1], 1, nchar(colnames(human.vsnrma)[i+1]))),
                     ".pdf", sep="")
  dev.copy2pdf(file = file.name)
  dev.off()
}




# 4) Data clean-up

# Are there any NAs in our dataset?
sum(apply(human.vsnrma.df,1,function(x){sum(is.na(x))}))
#> No




# 5) Data Analysis

# 5.1) Annotation

# read in ensembl table with following attributes: "Gene stable ID", 
#"Gene stable ID version", "Transcript stable ID", "Transscript stable ID version", 
#"Gene.name", "Transcript name", "Chromosome.scaffold.name", "Gene.description", "HGNC.symbol"

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

# Create a data frame out of the expression data from the normalized data
human.vsnrma.df = data.frame(exprs(human.vsnrma))

#rename the colnames to the stages
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456643.CEL'] <- '1-cell stage.rep.1'
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456644.CEL'] <- '1-cell stage.rep.2' 
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456645.CEL'] <- '1-cell stage.rep.3' 
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456646.CEL'] <- '2-cell stage.rep.1' 
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456647.CEL'] <- '2-cell stage.rep.2' 
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456648.CEL'] <- '2-cell stage.rep.3' 
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456649.CEL'] <- '4-cell stage.rep.1' 
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456650.CEL'] <- '4-cell stage.rep.2' 
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456651.CEL'] <- '4-cell stage.rep.3' 
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456652.CEL'] <- '8-cell stage.rep.1' 
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456653.CEL'] <- '8-cell stage.rep.2' 
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456654.CEL'] <- '8-cell stage.rep.3' 
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456655.CEL'] <- 'morula stage.rep.1' 
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456656.CEL'] <- 'morula stage.rep.2' 
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456657.CEL'] <- 'morula stage.rep.3'
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456658.CEL'] <- 'blastocyst stage.rep.1' 
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456659.CEL'] <- 'blastocyst stage.rep.2' 
names(human.vsnrma.df)[names(human.vsnrma.df) == 'GSM456660.CEL'] <- 'blastocyst stage.rep.3' 

# Check dimensions
dim(human.vsnrma.df)
#95,721    18

# exclude Affymetrix control genes which begin with "AFFX"
human.vsnrma.df2 = human.vsnrma.df[63:95721,]

dim(human.vsnrma.df2)
# 95,659    18

# remove ".xx_at" from the rownames
rownames(human.vsnrma.df2) = substr(rownames(human.vsnrma.df2), 1,15)

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





# 5.2) Exploratory data analysis


# Generate a list containing all tissue names with number of tissue TRAs, which is detected in our chips
tissue = c("Brain","Esophagus","Heart","Liver","Cervix","Muscle","Ovary","Colon","Breast","Kidney","Pituitary","Testis","Whole Blood","Cells - Transformed fibroblasts", "Cells - EBV-transformed lymphocytes", "Spleen","Lung","Stomach","Pancreas","Small Intestine","Skin","Artery","Adrenal","Gland","Pituitary","Nerve","Minor","Salivary","Bladder","Adipose","Thyroid","Prostate","Vagina")
tissue.number = c()

for (i in 1:33) {
  tissue.number[[i]] = nrow(fusion.tra.expression.tra.table.ensembl.table %>% filter(grepl(tissue[i],fusion.tra.expression.tra.table.ensembl.table[,21])))
}
print(tissue.number)
tissue.distribution = cbind(tissue,tissue.number)
tissue.distribution = as.data.frame(tissue.distribution)

# Arranging the list according to the descending TRAs number
tissue.distribution.arranged = arrange(tissue.distribution, desc(as.numeric(tissue.number)))

# Create a barplot
library(RColorBrewer)
coul <- brewer.pal(5, "Set2")
barplot(height = as.numeric(tissue.distribution.arranged$tissue.number[1:10]),names.arg = tissue.distribution.arranged$tissue[1:10],cex.names=0.8,col=coul, main="Frequency of TRAs in our data", ylab="Frequency")

setwd("~//documents//GitHub//2022-topic-04-team-03//Plots")

dev.copy2pdf(file="tissue.distribution.arranged.barplot.pdf" )


# Create a heatmap of the TRA genes in our dataset
## What does unlist mean? -> Transforms list into matrix
human.vsnrma.only.tra.matrix = matrix(unlist(human.vsnrma.only.tra), ncol = 6, nrow = 24783)

colnames(human.vsnrma.only.tra.matrix) = c("1 cell stage", "2 cell stage", "4 cell stage", "8 cell stage", "morula stage", "blastocyst stage")

heatmap = pheatmap(human.vsnrma.only.tra.matrix, 
                   cluster_cols=FALSE, 
                   show_rownames = TRUE, 
                   legend=TRUE, 
                   fontsize_row=0.5,
                   main = "Expression of 24,783 TRA genes in 6 different stages of embryogenesis")

setwd("~//documents//GitHub//2022-topic-04-team-03//Plots")
dev.copy2pdf(file="heatmap.pdf")


# 5.3) Limma analysis

#GSM456643 human embryo at 1 cell stage, biological rep 1
#GSM456645 human embryo at 1 cell stage, biological rep 2
#GSM456644 human embryo at 1 cell stage, biological rep 3
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


# Rename colnames of our dataset: "x stage, rep y" for a better overview
colnames(human.vsnrma.df2) = names

################# Our solution ##############################################

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
    limma.table.1.4 = topTable(fit.1,number=pvalue01,p.value = 0.01)
  }
  
  if (i==3){
    limma.table.1.8 = topTable(fit.1, number = pvalue01,p.value = 0.01)
  }
  
  if (i==4){
    limma.table.1.m = topTable(fit.1, number = pvalue01,p.value = 0.01)
  }
  
  if (i==5){
    limma.table.1.b = topTable(fit.1, number = pvalue01,p.value = 0.01)
  }
  
  if (i==6){
    limma.table.2.4 = topTable(fit.1, number = pvalue01,p.value = 0.01)
  }
  
  if (i==7){
    limma.table.2.8 = topTable(fit.1, number = pvalue01,p.value = 0.01)
  }
  
  if (i==8){
    limma.table.2.m = topTable(fit.1, number = pvalue01,p.value = 0.01)
  }
  
  if (i==9){
    limma.table.2.b = topTable(fit.1, number = pvalue01,p.value = 0.01)
  }
  
  if (i==10){
    limma.table.4.8 = topTable(fit.1, number = pvalue01,p.value = 0.01)
  }
  
  if (i==11){
    limma.table.4.m = topTable(fit.1, number = pvalue01,p.value = 0.01)
  }
  
  if (i==12){
    limma.table.4.b = topTable(fit.1, number = pvalue01,p.value = 0.01)
  }
  
  if (i==13){
    limma.table.8.m = topTable(fit.1, number = pvalue01,p.value = 0.01)
  }
  
  if (i==14){
    limma.table.8.b = topTable(fit.1, number = pvalue01,p.value = 0.01)
  }
  
  if (i==15){
    limma.table.m.b = topTable(fit.1, number = pvalue01,p.value = 0.01)
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
  colnames(annotated.x)[7:9]<-c("gene.name","gene.description","Chromosome.scaffold.name")
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










#setwd("~\\GitHub\\2022-topic-04-team-03")
#save.image(file="human_18290.bis.limma.RData")


# 5.4) Dimension Reduction using PCA

###code von Carl

#topVar = apply(human.vsnrma.df2, 1, var)
#q75 = quantile(topVar, probs = 0.75)
#i.topvar = which(topVar >= q75)
#human.vsnrma.df2.topVar = human.vsnrma.df2[i.topvar,]
#dim(human.vsnrma.df2.topVar)
#pca = prcomp(t(human.vsnrma.df2.topVar), center = F, scale. = F)
#print(pca)

#zeigt an welche PCs wievel standardabweichung erklären
#plot(pca$sdev)
#PCs anteile an gesamt Varianz
#variance = (pca$sdev)^2
#prop.variance = variance/sum(variance)
#names(prop.variance) = 1:length(prop.variance)
#barplot(prop.variance[1:20],ylab='Proportion of variance') # we only plot the first 20 PCs

#color = c(rep("red",3),rep("orange",3),rep("yellow",3),rep("green",3),rep("blue",3),rep("purple",3))

#plot(pca$x[,1], pca$x[,2],col=color, pch=19,xlab='PC1',ylab='PC2')


### code von Data camp
#pca = prcomp(t(human.vsnrma.df2), center = T, scale. = F)

# Variability of each principal component: pr.var
#pca.var <- pca$sdev^2

# Variance explained by each principal component: pve
#pve <- pca.var / sum(pca.var)

# Plot variance explained for each principal component
#plot(pve, xlab = "Principal Component",
     #ylab = "Proportion of Variance Explained",
     #ylim = c(0, 1), type = "b")

# Plot cumulative proportion of variance explained
#plot(cumsum(pve), xlab = "Principal Component",
     #ylab = "Cumulative Proportion of Variance Explained",
     #ylim = c(0, 1), type = "b")
## 4 PCs are enough to show more than 90%

#color = c(rep("red",3),rep("orange",3),rep("yellow",3),rep("green",3),rep("blue",3),rep("purple",3))

#plot(pca$x[,1], pca$x[,2],col=color, pch=19,xlab='PC1',ylab='PC2')
#plot(pca$x[,2], pca$x[,3],col=color, pch=19,xlab='PC2',ylab='PC3')
#plot(pca$x[,3], pca$x[,4],col=color, pch=19,xlab='PC3',ylab='PC4')


#code von Carl

topVar = apply(human.vsnrma.df2, 1, var)
q75 = quantile(topVar, probs = 0.75)
i.topvar = which(topVar >= q75)
human.vsnrma.df2.topVar = human.vsnrma.df2[i.topvar,]
dim(human.vsnrma.df2.topVar)

pca = prcomp(t(human.vsnrma.df2.topVar), center = T, scale. = T)
print(pca)

#zeigt an welche PCs wievel standardabweichung erklären
plot(pca$sdev)

#PCs anteile an gesamt Varianz
variance = (pca$sdev)^2
prop.variance = variance/sum(variance)
names(prop.variance) = 1:length(prop.variance)
barplot(prop.variance[1:20],ylab='Proportion of variance') # we only plot the first 20 PCs

# Plot cumulative proportion of variance explained
plot(cumsum(prop.variance), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     ylim = c(0, 1), type = "b")

abline(h=0.9)

#Zeigen von PC1 and PC2
color = c(rep("red",3),rep("orange",3),rep("yellow",3),rep("green",3),rep("blue",3),rep("purple",3))

plot(pca$x[,1], pca$x[,2],col=color, pch=19,xlab='PC1',ylab='PC2')

#Zeigen von PC2 and PC3
plot(pca$x[,2], pca$x[,3],col=color, pch=19,xlab='PC2',ylab='PC3')

#Zeigen von PC3 and PC4
plot(pca$x[,3], pca$x[,4],col=color, pch=19,xlab='PC3',ylab='PC4')

# 4 PCs are enough to show more than 90%

#Kmeans to cluster the stages
# Use map_dbl to run many models with varying value of k (centers)
tot_withinss <- map_dbl(1:10,  function(k){
  model <- kmeans(x = pca$x[,1:4], centers = k,nstart=30)
  model$tot.withinss
})

# Generate a data frame containing both k and tot_withinss
elbow_df <- data.frame(
  k = 1:10 ,
  tot_withinss = tot_withinss
)

# Plot the elbow plot
ggplot(elbow_df, aes(x = k, y =tot_withinss)) +
  geom_line() +
  scale_x_continuous(breaks = 1:10)
# seems like k=2 is the best

#show cluster
model <- kmeans(x = pca$x[,1:4], centers = 2,nstart=30)
model$cluster

#Plot the Siloutte plot
library(cluster)
library(factoextra)

# Use map_dbl to run many models with varying value of k
sil_width <- map_dbl(2:10,  function(k){
  model <- pam(pca$x[,1:4], k = k)
  model$silinfo$avg.width
})

# Generate a data frame containing both k and sil_width
sil_df <- data.frame(
  k = 2:10,
  sil_width = sil_width
)

# Plot the relationship between k and sil_width
ggplot(sil_df, aes(x = k, y = sil_width)) +
  geom_line() +
  scale_x_continuous(breaks = 2:10)

#it seems like that k=4 is a good choice

# Generate a k-means model using the pam() function with a k = 4
pam_k4 <- pam(pca$x[,1:13],k=4)

# Show the clustering
pam_k4$clustering

#clustering without PCA
# Use map_dbl to run many models with varying value of k (centers)
tot_withinss <- map_dbl(1:10,  function(k){
  model <- kmeans(x = t(human.vsnrma.df, centers = k,nstart=30)
                  model$tot.withinss
})
  
  # Generate a data frame containing both k and tot_withinss
  elbow_df <- data.frame(
    k = 1:10 ,
    tot_withinss = tot_withinss
  )
  
  # Plot the elbow plot
  ggplot(elbow_df, aes(x = k, y =tot_withinss)) +
    geom_line() +
    scale_x_continuous(breaks = 1:10)
  # seems like k=2 is the best
  
  #show cluster
  model <- kmeans(x = t(human.vsnrma.df2), centers = 2,nstart=30)
  model$cluster
  
  #Plot the Siloutte plot
  library(cluster)
  library(factoextra)
  
  # Use map_dbl to run many models with varying value of k
  sil_width <- map_dbl(2:10,  function(k){
    model <- pam(t(human.vsnrma.df2), k = k)
    model$silinfo$avg.width
  })
  
  # Generate a data frame containing both k and sil_width
  sil_df <- data.frame(
    k = 2:10,
    sil_width = sil_width
  )
  
  # Plot the relationship between k and sil_width
  ggplot(sil_df, aes(x = k, y = sil_width)) +
    geom_line() +
    scale_x_continuous(breaks = 2:10)
  
  #it seems like that k=2 is a good choice
  
  # Generate a k-means model using the pam() function with a k = 4
  pam_k2 <- pam(t(human.vsnrma.df2),k=2)
  
  # Show the clustering
  pam_k2$clustering
  
  ###Finding Variance Genes
  pc1=pca$rotation[,1]
  top_genes_pca1 = names(sort(abs(pc1), decreasing = T)[1:10])
  k = which(rownames(fusion.tra.expression.tra.table.ensembl.table)%in% top_genes_pca1)
  fusion.tra.expression.tra.table.ensembl.table[k,]
  
  j = which(ensembl.data$Transcript.stable.ID%in% top_genes_pca1)
  ensembl.data[j,]
  
  pc2=pca$rotation[,2]
  top_genes_pca2 = names(sort(abs(pc2), decreasing = T)[1:10])
  k = which(rownames(fusion.tra.expression.tra.table.ensembl.table)%in% top_genes_pca2)
  fusion.tra.expression.tra.table.ensembl.table[k,]
  
  j = which(ensembl.data$Transcript.stable.ID%in% top_genes_pca2)
  ensembl.data[j,]
  
  pc3=pca$rotation[,3]
  top_genes_pca3 = names(sort(abs(pc3), decreasing = T)[1:10])
  k = which(rownames(fusion.tra.expression.tra.table.ensembl.table)%in% top_genes_pca3)
  fusion.tra.expression.tra.table.ensembl.table[k,]
  
  j = which(ensembl.data$Transcript.stable.ID%in% top_genes_pca3)
  ensembl.data[j,]
  
  pc4=pca$rotation[,4]
  top_genes_pca4 = names(sort(abs(pc4), decreasing = T)[1:10])
  k = which(rownames(fusion.tra.expression.tra.table.ensembl.table)%in% top_genes_pca4)
  fusion.tra.expression.tra.table.ensembl.table[k,]
  
  j = which(ensembl.data$Transcript.stable.ID%in% top_genes_pca4)
  ensembl.data[j,]
  
  #1-2-4 cell vergleich mit 8 cell
  pca.1.8=prcomp(t(human.vsnrma.df2[,1:12]))
  #zeigt an welche PCs wievel standardabweichung erklären
  plot(pca.1.8$sdev)
  
  #PCs anteile an gesamt Varianz
  variance = (pca.1.8$sdev)^2
  prop.variance = variance/sum(variance)
  names(prop.variance) = 1:length(prop.variance)
  barplot(prop.variance[1:20],ylab='Proportion of variance') # we only plot the first 20 PCs
  
  # Plot cumulative proportion of variance explained
  plot(cumsum(prop.variance), xlab = "Principal Component",
       ylab = "Cumulative Proportion of Variance Explained",
       ylim = c(0, 1), type = "b")
  
  abline(h=0.9)
  
  color = c(rep("red",3),rep("orange",3),rep("yellow",3),rep("green",3),rep("blue",3),rep("purple",3))
  
  plot(pca.1.8$x[,1], pca.1.8$x[,2],col=color, pch=19,xlab='PC1',ylab='PC2')
  
  ###Finding Variance Genes
  pc1=pca.1.8$rotation[,1]
  top_genes_pca1 = names(sort(abs(pc1), decreasing = T)[1:10])
  k = which(rownames(fusion.tra.expression.tra.table.ensembl.table)%in% top_genes_pca1)
  fusion.tra.expression.tra.table.ensembl.table[k,]
  
  j = which(ensembl.data$Transcript.stable.ID%in% top_genes_pca1)
  View(ensembl.data[j,])
  
  #8-cell vergleicht mit m 
  pca.8.m=prcomp(t(human.vsnrma.df2[,10:15]))
  #zeigt an welche PCs wievel standardabweichung erklären
  plot(pca.8.m$sdev)
  
  #PCs anteile an gesamt Varianz
  variance = (pca.8.m$sdev)^2
  prop.variance = variance/sum(variance)
  names(prop.variance) = 1:length(prop.variance)
  barplot(prop.variance[1:20],ylab='Proportion of variance') # we only plot the first 20 PCs
  
  # Plot cumulative proportion of variance explained
  plot(cumsum(prop.variance), xlab = "Principal Component",
       ylab = "Cumulative Proportion of Variance Explained",
       ylim = c(0, 1), type = "b")
  
  abline(h=0.9)
  
  color = c(rep("red",3),rep("orange",3),rep("yellow",3),rep("green",3),rep("blue",3),rep("purple",3))
  
  plot(pca.8.m$x[,1], pca.8.m$x[,2],col=color, pch=19,xlab='PC1',ylab='PC2')
  
  ###Finding Variance Genes
  pc1=pca.8.m$rotation[,1]
  top_genes_pca1 = names(sort(abs(pc1), decreasing = T)[1:10])
  k = which(rownames(fusion.tra.expression.tra.table.ensembl.table)%in% top_genes_pca1)
  fusion.tra.expression.tra.table.ensembl.table[k,]
  
  j = which(ensembl.data$Transcript.stable.ID%in% top_genes_pca1)
  View(ensembl.data[j,])
  
  #m vergleicht mit b
  pca.m.b=prcomp(t(human.vsnrma.df2[,13:18]))
  #zeigt an welche PCs wievel standardabweichung erklären
  plot(pca.m.b$sdev)
  
  #PCs anteile an gesamt Varianz
  variance = (pca.m.b$sdev)^2
  prop.variance = variance/sum(variance)
  names(prop.variance) = 1:length(prop.variance)
  barplot(prop.variance[1:20],ylab='Proportion of variance') # we only plot the first 20 PCs
  
  # Plot cumulative proportion of variance explained
  plot(cumsum(prop.variance), xlab = "Principal Component",
       ylab = "Cumulative Proportion of Variance Explained",
       ylim = c(0, 1), type = "b")
  
  abline(h=0.9)
  
  color = c(rep("red",3),rep("orange",3),rep("yellow",3),rep("green",3),rep("blue",3),rep("purple",3))
  
  plot(pca.m.b$x[,1], pca.m.b$x[,2],col=color, pch=19,xlab='PC1',ylab='PC2')
  
  ###Finding Variance Genes
  pc1=pca.8.m.b$rotation[,1]
  top_genes_pca1 = names(sort(abs(pc1), decreasing = T)[1:10])
  k = which(rownames(fusion.tra.expression.tra.table.ensembl.table)%in% top_genes_pca1)
  fusion.tra.expression.tra.table.ensembl.table[k,]
  
  j = which(ensembl.data$Transcript.stable.ID%in% top_genes_pca1)
  ensembl.data[j,]
  
  
  # 6) filter DE genes
  ###between 1 cell - 8 cell stadium
  
  hist(annotated.limma.1.8$logFC)
  #abline(v=mean(abs(limma.table.1.8$logFC)))#
  #c=which(limma.table.1.8$logFC< 0)
  #limma.1.8.down=limma.table.1.8[c,]
  #limma.1.8.down.sig=rownames(annotated.limma.1.8.down[which(annotated.limma.1.8.down$P.Value < 0.05),])
  #limma.table.1.8.arranged=arrange(limma.table.1.8,desc(abs(logFC)))
  #limma.1.8.down.sig=arrange(limma.1.8.down.sig,)
 
  #filter the top 20 DE upregulated genes
c=which(annotated.limma.1.8$logFC > 2)
limma.table.1.8.up=annotated.limma.1.8[c,]
limma.1.8.up.sig=limma.table.1.8.up[which(limma.table.1.8.up$P.Value < 0.05),]
limma.1.8.up.sig=arrange(limma.1.8.up.sig,desc(limma.1.8.up.sig$logFC))  
interesting_genes.1.8=rownames(limma.1.8.up.sig[1:20,])

k = which(rownames(fusion.tra.expression.tra.table.ensembl.table)%in% interesting_genes.1.8)
  View(fusion.tra.expression.tra.table.ensembl.table[k,])
  
j = which(ensembl.data$Transcript.stable.ID%in%interesting_genes.1.8 )
  View(ensembl.data[j,])
  
  #between 8 cell -  m stadium
  hist(limma.table.8.m$logFC)
  ###filter the top20 DE upregulated genes
  
  #c=which(limma.table.8.m$logFC< -3)
  #limma.table.8.m.down=limma.table.8.m[c,]
  #limma.8.m.down.sig=rownames(limma.table.8.m.down[which(limma.table.1.8.down$P.Value < 0.05),])
  #c=which(limma.table.8.m$logFC > 4)
  
  c=which(annotated.limma.8.m$logFC > 2)
  limma.table.8.m.up=annotated.limma.8.m[c,]
  limma.8.m.up.sig=limma.table.8.m.up[which(limma.table.8.m.up$P.Value < 0.05),]
  limma.8.m.up.sig=arrange(limma.8.m.up.sig,desc(limma.8.m.up.sig$logFC))  
  
  #Transcript of Top 20 genes
  interesting_genes.8.m=rownames(limma.8.m.up.sig[1:20,])
  #Name of Top 20 genes
  limma.8.m.up.sig[1:20,]$`ensembl.new$Gene.name`
  
  k = which(rownames(fusion.tra.expression.tra.table.ensembl.table)%in% interesting_genes.8.m)
  View(fusion.tra.expression.tra.table.ensembl.table[k,])
  
  j = which(ensembl.data$Transcript.stable.ID%in%interesting_genes.8.m )
  View(ensembl.data[j,])
  
  #between m - b stadium
  
  ###filter the top 20 regulated genes
  
  c=which(annotated.limma.m.b$logFC > 2)
  limma.table.m.b.up=annotated.limma.m.b[c,]
  limma.m.b.up.sig=limma.table.m.b.up[which(limma.table.m.b.up$P.Value < 0.05),]
  limma.m.b.up.sig=arrange(limma.m.b.up.sig,desc(limma.m.b.up.sig$logFC))  
  interesting_genes.m.b=rownames(limma.m.b.up.sig[1:20,])
  
  k = which(rownames(fusion.tra.expression.tra.table.ensembl.table)%in% interesting_genes.m.b)
  View(fusion.tra.expression.tra.table.ensembl.table[k,])
  
  j = which(ensembl.data$Transcript.stable.ID%in%interesting_genes.m.b )
  View(ensembl.data[j,])
  
  ###Volcanoplot
  library(EnhancedVolcano)
  EnhancedVolcano(annotated.limma.4.8[,1:6],lab=annotated.limma.4.8[,7],x = 'logFC',y = 'P.Value', pCutoff = 0.05,
                  FCcutoff = 3,)
  EnhancedVolcano(annotated.limma.8.m[,1:6],lab=annotated.limma.8.m[,7],x = 'logFC',y = 'P.Value', pCutoff = 10e-8,
                  FCcutoff = 2,)
  EnhancedVolcano(annotated.limma.m.b[,1:6],lab=annotated.limma.m.b[,7],x = 'logFC',y = 'P.Value', pCutoff = 10e-8,
                  FCcutoff = 2,)
  
  
  
  
  
  fit.8.m = contrasts.fit(fit, contrasts = cm[,13])
  fit.8.m = eBayes(fit.8.m)
  volcanoplot(fit.8.m,highlight = 10)
  
  fit.m.b = contrasts.fit(fit, contrasts = cm[,15])
  fit.m.b = eBayes(fit.m.b)
  volcanoplot(fit.8.m,highlight = 10)
  