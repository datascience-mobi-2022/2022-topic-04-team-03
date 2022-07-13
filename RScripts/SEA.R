BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

setwd("~//documents//GitHub//2022-topic-04-team-03")

# reading in data from deseq2
df = read.csv("drosphila_example_de.csv", header=TRUE)

# we want the log2 fold change 
GO.1.8.genes<- annotated.limma.1.8$logFC

# name the vector
names(GO.1.8.genes) <- annotated.limma.1.8$Gene.stable.ID


# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(GO.1.8.genes, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)




