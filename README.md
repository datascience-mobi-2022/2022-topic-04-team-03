# 2022-topic-04-team-03

This is our project

# 1.Abstract
Research, especially in the field of early embryogenesis is important in many ways. It is not only required for a better understanding of genetic defects of the embryos in an early stage but also important for the improvement of assisted human conception. In this report we investigated gene expression data of the early stages of embryogenesis up to the formation of the blastocyst. As the formation of tissues only occurs after the stages displayed by our dataset, we expected to not observe any tissue specific development and therefore no differentially expressed tissue-restricted antigens (TRAs) during this period. However, we found numerous differentially expressed TRAs that display the developmental processes in early embryogenesis accurately. For instance, the inquired genes reflect the activation of transcription of the embryonic genome as well as the switch of energy metabolism in the early embryo.
Furthermore, we found three differentially expressed chemokines that play a role in implantation and embryonic development. 
These insights may contribute to a deeper understanding of embryo implantation and therefore successful conception.

# 2. Used data & packages 

* Geo accession to the analysed data set: GSE18290 ([link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18290))
* TRA data set: GTEX RNAseq human dataset provided by Dr. Dinkelacker in 2019 (60,131 transcripts)


List of all R and Bioconducter packages used in our project
  Packages         |                      |                 |                            |
-------------      | -------------        |-------------    | -------------              |
ggplot2 3.3.6      | AnnotationDbi 1.58.0 | dplyr 1.0.9     | hgu133plus2hsenstcdf 25.0.0|
-------------      | -------------        |-------------    | -------------              |
affy 1.74.0        | vsn 3.64.2.          | tidyverse 1.3.1 | hgu133plus2hsenstprobe     |
-------------      | -------------        |-------------    | -------------              |
VennDiagram 1.7.3  | limma 3.52.2         | hexbin 1.28.2   | org.Hs.eg.d 3.15.0         |
-------------      | -------------        |-------------    | -------------              |
cluster 2.1.3      | factoextra 1.0.7     | go.db 3.15.0    | pheatmap 1.0.12            |
-------------      | -------------        |-------------    | -------------              |
imager 0.42.13     | RColorBrewer 1.1-3   | tinytex 0.40    | clusterProfiler 4.4.4      |
-------------      | -------------        |-------------    | -------------              |
enrichplot 1.16.1  |
