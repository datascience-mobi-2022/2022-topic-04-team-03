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

# 3.2) Normalization

human.vsnrma <- vsnrma(data.human)
human.vsnrma.df = data.frame(exprs(human.vsnrma))
human.vsnrma.df2 = human.vsnrma.df[63:95721,]
rownames(human.vsnrma.df2) = substr(rownames(human.vsnrma.df2), 1,15)
