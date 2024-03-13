- RStudio
- R packages
- ```
  if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install("maser")
BiocManager::install("Rgraphviz")
BiocManager::install("lintr")
install.packages("knitr")
```
