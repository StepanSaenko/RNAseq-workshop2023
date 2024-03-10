# Gevol RNAseq workshop 2023/24

- FileZilla or another equivalent tool to transfer remote files
- WSL, Git bash, PuttY, mobaXterm or another terminal with ssh capability
- Rstudio
- Server: ssh -X -p 6022 user@131.220.238.3

## 1. Install libraries:  
For the purpose of this tutorial, there are some packages, which each of you
will need to install. You can do so using the following commands:  
```
# install bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install(c("BiocFileCache","xml2","AnnotationDbi","XML","httr","KEGGREST","lintr","lattice","Rgraphviz"))
BiocManager::install("biomaRt")
BiocManager::install("topGO")
BiocManager::install("tximport")
BiocManager::install("tidyverse")
BiocManager::install("eulerr")
BiocManager::install("pheatmap")
BiocManager::install("RColorBrewer")
BiocManager::install("DESeq2")
install.packages('CoDiNA')
install.packages("devtools")
install.packages("ncdf4")
install.packages(c("HiClimR","igraph"))

library(devtools)
install_github("deisygysi/wTO")
devtools::install_github('RfastOfficial/Rfast')
```
