#RNAseq analysis with kallisto and DeSeq2

### install R packages
```
#install bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install(c("BiocFileCache","xml2","AnnotationDbi","XML","httr","KEGGREST","lintr","lattice","Rgraphviz"))
BiocManager::install("biomaRt")
BiocManager::install("topGO")
install.packages('tidyverse', dependencies=TRUE)
install.packages(c("tximport","eulerr","readr","ggplot2","dplyr","tidyr", dependecies=TRUE))
install.packages('CoDiNA')
install.packages("devtools")
install.packages("ncdf4")
#if you run Linux you can install it also via: sudo apt-get install r-cran-ncdf4
install.packages(c("HiClimR","igraph"))

library(devtools)
install_github("deisygysi/wTO")
devtools::install_github('RfastOfficial/Rfast')
#update RcppArmadillo

BiocManager::install("tximport")
BiocManager::install("tidyverse")
BiocManager::install("eulerr")
BiocManager::install("pheatmap")
BiocManager::install("RColorBrewer")
BiocManager::install("DESeq2")
```


### start analysis of transcript counts from kallisto
```
# load packages
library(tximport)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)

#check you working DIR
getwd()

# get list of PATHs/Files (one output file per Sample

# empty list
paths <- list()

# which folder to scan
paths$count_output   <- "2024-03-06.kallisto.output"

# find the kallisto outputfiles (for salmon outputs change .tsv to .sf)
paths$count_files_relative <- grep(x = list.files(paths$count_output,
                                                     recursive = TRUE),
                                      pattern = ".tsv",
                                      value   = TRUE)

# full PATHs for each file to read later
paths$count_files <- file.path(paths$count_output,
                                  paths$count_files_relative)

# Automatically extract file names from salmon/kallisto outputs to use for colnames
# names(paths$salmon_files) <- paste0("sample", 1:17) #standard sample names sample1 ... sample17
names(paths$salmon_files) <- gsub(paths$count_files_relative,
                                    pattern = "/.*",
                                    replacement = "")

# Extract sample names and put into df for tximport
samples     <- data.frame(treatment = names(paths$count_files))

# load sample information, info about treatments or other groupings
samples_information <- read.table(file = "sample.data.short.tsv",
                                  header = FALSE,
                                  col.names = c("sample_name",
                                                "age",
                                                "sample",
                                                "treatment",
                                                "caste",
                                                "tissue"),
                                  row.names = 1)

# format some columns/data
samples_information$tissue <- as.factor(samples_information$tissue)
samples_information$caste <- as.factor(samples_information$caste)
samples_information$age <- as.factor(samples_information$age)
samples_information$sample_id <- row.names(samples_information)
print(samples_information)
#samples <- data.frame(treatment = samples_information[4])

# load transcripts to gene mapping
transcript_to_gene <- read.table(file = "transcript.gene.map.txt",
                                 col.names = c("transcript",
                                               "locus"))

# import data
txi_counts <- tximport(paths$count_files,
                       type    = "kallisto",
                       tx2gene = transcript_to_gene,
                       countsFromAbundance = "no")

#check
head(txi_counts$counts)

#create Deseq2 objects
print(samples_information)

deseq_txi1 <- DESeqDataSetFromTximport(txi     = txi_counts,
                                      colData = samples_information,
                                      design  = ~age)

deseq_txi2 <- DESeqDataSetFromTximport(txi     = txi_counts,
                                      colData = samples_information,
                                      design  = ~age + tissue)

deseq_txi3 <- DESeqDataSetFromTximport(txi     = txi_counts,
                                      colData = samples_information,
                                      design  = ~age + caste + tissue)

deseq_txi4 <- DESeqDataSetFromTximport(txi     = txi_counts,
                                      colData = samples_information,
                                      design  = ~age + caste)
# Error in checkFullRank(modelMatrix) :   the model matrix is not full rank, so the model cannot be fit as specified.
#  One or more variables or interaction terms in the design formula are linear
#  combinations of the others and must be removed.
# but age is not ressolved by caste, queens are 7 and drones are 15d old


deseq_txi <- DESeqDataSetFromTximport(txi     = txi_counts,
                                      colData = samples_information,
                                      design  = ~caste + tissue)
```

### analysis
```
## how many genes did we detect
nrow(counts(deseq_txi))
#12123
```

```
## Filter out lowly expressed genes:
keep      <- rowSums(counts(deseq_txi)) >= 10
deseq_txi <- deseq_txi[keep, ]
print(nrow(deseq_txi))
#11980
```

```
## Create output 'results' directory:
dir.create(path = "results")
saveRDS(object = deseq_txi,
        file = "results/deseq_txi.rds")
```

#create Deseq2 object
```
deseq_object  <- DESeq(deseq_txi)
```
## 7. Extract normalised gene-level counts:
For examining relationships among samples and tissues, we normalise our read counts, which makes samples more comparable. The form of normalisation is called 'variance stablising transformation' and is implemented using the function 'vst()' from the DESeq2 R package.
```
## Normalise counts:
vsd <- vst(deseq_object,
           blind = FALSE)
```

## 8. Principal component analysis:
Using the normalised counts stored in the 'vsd' object, we can perform a principal component analysis (PCA). PCAs help to visualise the overall similarities and differences among samples allowing for identifying outlier samples or batch effects (e.g., whether individuals from the same colony are more similar to each other than to other ants).  

To perform a PCA, we run the function 'plotPCA()', which performs both the PCA but also plots the output of the analysis. By default, it outputs the first two principal components, which each explain a certain amount of variance in our dataset. 
```
## PCA
pcaplot <- plotPCA(vsd,
        intgroup = c("tissue",
                     "caste")) +
  theme_bw()

jpeg(file="PCA.jpeg")
pcaplot
dev.off()
```

```
## heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$submorph,
                                    vsd$tissue,
                                    sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9,
                                          "Blues")))(255)  

heatmap <- pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

jpeg(file="Heatmap.jpeg")
heatmap
dev.off()
```

## save subsets of the data

```
gonads_samples  <- colnames(counts(deseq_txi)) %in% 
  subset(x = samples_information,
                tissue == "gonad")$sample_id
deseq_txi_gonads <- deseq_txi[, gonads_samples]

## Check if the number subsetted is correct:
print(nrow(deseq_txi_gonads))

## Save as an RDS object:
saveRDS(object = deseq_txi_gonads,
        file = "results/deseq_txi_gonads.rds")


brain_samples  <- colnames(counts(deseq_txi)) %in% 
  subset(x = samples_information,
                tissue == "brains")$sample_id
deseq_txi_brains <- deseq_txi[, brain_samples]

## Check if the number subsetted is correct:
print(nrow(deseq_txi_brains))

## Save:
saveRDS(object = deseq_txi_brains,
        file = "results/deseq_txi_brains.rds")


queen_samples  <- colnames(counts(deseq_txi)) %in% 
  subset(x = samples_information,
                caste == "queen")$sample_id
deseq_txi_queens <- deseq_txi[, queen_samples]

## Check if the number subsetted is correct:
print(nrow(deseq_txi_queens))

## Save:
saveRDS(object = deseq_txi_queens,
        file = "results/deseq_txi_queens.rds")


drone_samples  <- colnames(counts(deseq_txi)) %in% 
  subset(x = samples_information,
                caste == "drone")$sample_id
deseq_txi_drones <- deseq_txi[, drone_samples]

## Check if the number subsetted is correct:
print(nrow(deseq_txi_drones))

## Save:
saveRDS(object = deseq_txi_drones,
        file = "results/deseq_txi_drones.rds")

```

## further analysis specific to each subset
```

# load the Deseq Object we created earlier
deseq_txi_brains <- readRDS(file = "results/deseq_txi_brains.rds")

# modify the design
design(deseq_txi_brains) <- ~caste + tissue

deseq_txi_brains$caste <- relevel(deseq_txi_brains$caste,
                                       ref = "Q")

deseq_object <- DESeq(deseq_txi_brains,
                      test = "LRT",
                      reduced = ~caste)




```
















