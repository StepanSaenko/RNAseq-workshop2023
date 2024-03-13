# Gevol RNAseq workshop 12./13.03.2024
Eckart stolle (LIB Museum Koenig)
Lars Podsiadlowski (LIB Museum Koenig)
Katja Nowick (FU Berlin)
Joe Colgan (U Mainz)

## Day1: RNAseq, read processing, mapping, read counts, DEG, GO term enrichment (E. Stolle, L. Podsiadlowski, J. Colgan)
## Day2a: co-expression Networks, comparisons of networks (K. Nowick)
## Day2b: differential splicing (J. Colgan)


- Terminal which can do ssh, e.g. WSL, Git bash, PuttY, mobaXterm 
- FileZilla or another equivalent tool to transfer remote files (if you are not using trynsfer via ssh (eg rsync)
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

install.packages("gplots")
```
