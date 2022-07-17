# 2022-topic-04-team-03
## The role of tissue-restricted antigens in the early embryonic development
-------------  
**Team:** Clara Certa, Yaxin Chen, Linda Kaupp, Alewtina Towara
**Supervisor:** Dr. Maria Dinkelacker, Dr. Carl Herrmann 
**Tutor:** Ian Dirk Fichtner 

# 1. Abstract
Research, especially in the field of early embryogenesis is important in many ways. It is not only required for a better understanding of genetic defects of the embryos in an early stage but also important for the improvement of assisted human conception. In this report we investigated gene expression data of the early stages of embryogenesis up to the formation of the blastocyst. As the formation of tissues only occurs after the stages displayed by our dataset, we expected to not observe any tissue specific development and therefore no differentially expressed tissue-restricted antigens (TRAs) during this period. However, we found numerous differentially expressed TRAs that display the developmental processes in early embryogenesis accurately. For instance, the inquired genes reflect the activation of transcription of the embryonic genome as well as the switch of energy metabolism in the early embryo.
Furthermore, we found three differentially expressed chemokines that play a role in implantation and embryonic development. 
These insights may contribute to a deeper understanding of embryo implantation and therefore successful conception.

# 2. Used data & packages 

* The microarray dataset by Xie et al. contained three datasets for three species: bovine, mouse and human. We decided to focus our project on the human dataset. The data was generated by RNA extraction, amplification, and hybridisation onto Affymetrix microarrays. For early human embryogenesis, they prepared three replicates at 1-cell stage, 2-cell stage, 4-cell stage, 8-cell stage, morula, and blastocyst each. Each replicate includes 95,659 transcripts.
Geo accession to the analysed data set: GSE18290 ([link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18290))

* TRA data set: GTEX RNAseq human dataset provided by Dr. Dinkelacker in 2019 (60,131 transcripts from 53 different tissues.


List of all R and Bioconducter packages used in our project
used packages      |                      |                 |                            |
-------------      | -------------        |-------------    | -------------              |
ggplot2 3.3.6      | AnnotationDbi 1.58.0 | dplyr 1.0.9     | hgu133plus2hsenstcdf 25.0.0|
affy 1.74.0        | vsn 3.64.2.          | tidyverse 1.3.1 | hgu133plus2hsenstprobe     |
VennDiagram 1.7.3  | limma 3.52.2         | hexbin 1.28.2   | org.Hs.eg.d 3.15.0         |
cluster 2.1.3      | factoextra 1.0.7     | go.db 3.15.0    | pheatmap 1.0.12            |
imager 0.42.13     | RColorBrewer 1.1-3   | tinytex 0.40    | clusterProfiler 4.4.4      |
enrichplot 1.16.1  |

# 3. Folder structure

- [Data](https://github.com/datascience-mobi-2022/2022-topic-04-team-03/tree/main/Data): This folder contains .CEL files of the raw data and normalized data.
- [Plots](https://github.com/datascience-mobi-2022/2022-topic-04-team-03/tree/main/Plots):  This folder contains the plots generated in our data analysis
- [Quality-Control](https://github.com/datascience-mobi-2022/2022-topic-04-team-03/tree/main/Quality-Control): This folder contains pictures of the chips and all plots generated during our QC.
- [Tables](https://github.com/datascience-mobi-2022/2022-topic-04-team-03/tree/main/Tables): This folder contains our employed and generated tables.
- [Presentations and report](https://github.com/datascience-mobi-2022/2022-topic-04-team-03/tree/main/Presentations%20and%20report): This folder contains the PowerPoints of our project proposal and our final presenatatin. Additionally it contains the report cache and figuere-latex.
