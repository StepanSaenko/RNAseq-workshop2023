# Gevol RNAseq workshop 2023/24: Networks, Katja Nowick, 13.03.2024

- FileZilla or another equivalent tool to transfer remote files
- WSL, Git bash, PuttY, mobaXterm or another terminal with ssh capability
- Rstudio
- download all files of the github folder "networks" locally to be used in RStudio


### Install R, Rtools, RStudio and inside R the following packages
```
install.packages('CoDiNA')
install.packages("devtools")
library(devtools)
install_github("deisygysi/wTO")

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




```


