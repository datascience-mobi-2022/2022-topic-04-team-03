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

### 8. Data analysis

## Create dataframe

df.norm <- as.data.frame(exprs(human.vsnrma))

#Remove .cel from rows
rownames(df.norm) <- gsub("\\..*$" ,"", rownames(df.norm)) # Strip file endings

#Import IL Gene table from ex 2
table.il <- DataFrame(read.csv("../../Tables/exercise_2_table_ILgenes.csv", row.names = 3))
#rownames(table.il) <- table.il$transID


## Filter microarray data for IL genes

df.norm.filt <- df.norm[rownames(df.norm)%in%rownames(table.il),]

# Translocate
df.norm.filt.t <- t(df.norm.filt)

# Regive lost column name
colnames(df.norm.filt.t) <- table.il[rownames(df.norm.filt), "trans_name"]


## Order column names
df.norm.filt.t <- df.norm.filt.t[, order(colnames(df.norm.filt.t))]


## Create column for facetting

map.align <- tibble(Transcript = colnames(df.norm.filt.t), align =(1:dim(df.norm.filt.t)[2] %/% 100)) # A column with separators: 0, 1, 2
df.norm.filt.longer <- gather(as_tibble(df.norm.filt.t), key="Transcript", value = "Expression") # Merging all transcripts into a single column with a separate key column
df.norm.merge <- merge(df.norm.filt.longer, map.align, on="Transcript") # Merging our data frame with the map.align column. We lose chip donor but its ok,

## Plot

fig1 <- ggplot(filter(df.norm.merge, align==0 | align==1)) +
  geom_boxplot(aes(Transcript, Expression, fill=Transcript), show.legend = F, outlier.size = 0.6, outlier.shape= 16, outlier.stroke = 0, lwd=0.1) +
  facet_wrap(.~align, scales="free_x", ncol = 1) + #nur eine Spalte
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 7),
    strip.text = element_blank() # Delete facetting text
  ) + 
  expand_limits(x=100) +
  labs(x = "Transcript",
       title = "Distribution of interleukin expression in breast cancer",
       subtitle = "GEO dataset: GSE27830")

print(fig1)


fig2 <- ggplot(filter(df.norm.merge, align==2 | align==3)) +
  geom_boxplot(aes(Transcript, Expression, fill=Transcript), show.legend = F, outlier.size = 0.6, outlier.shape= 16, outlier.stroke = 0, lwd=0.1) +
  facet_wrap(.~align, scales="free_x", ncol = 1) + #nur eine Spalte
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = 7),
    strip.text = element_blank() # Delete facetting text
  ) + 
  expand_limits(x=100) +
  labs(x = "Transcript",
       title = "Distribution of interleukin expression in breast cancer",
       subtitle = "GEO dataset: GSE27830")

print(fig2)

## Save plot
pdf(file="../../Plots/exercise_7_Gene_expression_IL_breast_cancer.pdf", height = 7, width = 10)
print(fig1)
print(fig2)
dev.off()

### Save data

# Expression table with annotations
hg104 <- read.csv("../../Raw-Data/ensembl_human_104.csv") #human genome

## csv: All microarray genes
gene.exp <- data.frame(exprs(breast.vsnrma))
gene.exp["Transcript.stable.ID"] <- gsub("[\\._].*$" ,"", rownames(gene.exp)) # Prep for merging

write.csv(
  right_join(hg104, gene.exp, by="Transcript.stable.ID"), 
  "../../Tables/exercise_7_Microarray_Genes.csv",
  quote=F,
  row.names = F
)

## csv: Only IL microarray genes
il.exp <- data.frame(t(df.norm.filt.t)[order(rownames(t(df.norm.filt.t))),])

# Prepare extra column for merging
il.exp["Transcript.name"] <- rownames(il.exp) 

write.csv(right_join(hg104, il.exp, by="Transcript.name"),
          "../../Tables/exercise_7_Micorarray_Ilgenes.csv",
          quote=F,
          row.names=F
)

save.image("../RDA-Files/exercise_7.rda")
