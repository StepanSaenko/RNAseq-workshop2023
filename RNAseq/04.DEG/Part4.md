# Part 4. overlap
The purpose of this script is to examine overlap in genes that are differentially expressed between queens and
drones in the gonads and the brains.
The script loads RDS objects generated from previous scripts:
- deseq_object_results_sig_gonads.rds
- deseq_object_results_sig_brains.rds
The scripts examines overlap between differentially expressed genes in both tissues and generates a Euler plot
to visualise genes that share similar expression profiles between nurses and foragers.

## 1. Load libraries
```
library(tidyverse)
library(DESeq2)
library(eulerr)
```

## 2. Load datasets:
We load our outputs for each tissue from DESeq2, which contain genes (rownames), baseMean (mean
expression of a gene across all samples), log2FoldChange (expression difference between nurses and foragers),
as well as adjusted p value (padj).
```
deseq_object_results_sig_gonads <- readRDS(
file = "results/deseq_object_results_sig_gonads.rds")
deseq_object_results_sig_brains <- readRDS(
file = "results/deseq_object_results_sig_brains.rds")
```

## 3. Extract significantly differentially expressed genes:
As mentioned above, the row names of our DESeq2 output contain the gene name, which we can extract
using the function ‘rownames()’ and store as an object for each individual tissue.
```
gonads_degs <- rownames(deseq_object_results_sig_gonads)
brains_degs <- rownames(deseq_object_results_sig_brains)
```

## 4. Combine input genes lists to examine overlap:
Next, we use the ‘euler()’ function from the eulerr R package, which creates a list using the DEGs found in each tissue.
```
euler_lists <- euler(combinations = list("Gonads" = gonads_degs,
"Brains" = brains_degs))
```


## 5. Generate plot:
We can visualise the overlap between DEGs using the ‘plot()’ function.
```
euler_all_plot <- plot(euler_lists,
quantities = TRUE,
edges = TRUE,
labels = list(fontsize = 8))

## Plot to console:
euler_all_plot
```

The last thing to do is save the image to file:
```
ggsave("results/euler_all_plot.png",
plot = euler_all_plot,
height = 5,
width = 5)

```

## 6. Examine changes per caste:
Next, we want to look at specific overlap within a caste/sex.

#### 1) Gonads
##### a) drones:
For drones, we subset the output of DESeq2 to obtain genes that are both significantly expressed (padj < 0.05) and have a bias in terms of expression in drones (log2FoldChange > 0).
Once we subset the data we want, we extract the gene names and store as an object.
```
gonads_degs_drones <- subset(x = deseq_object_results_sig_gonads,
padj < 0.05 & log2FoldChange > 0)
gonads_degs_drones_list <- row.names(gonads_degs_drones)
```

##### b) queens:
Same with queens, we subset the output of DESeq2 to obtain genes that are both significantly expressed
(padj < 0.05) and have a bias in terms of expression in queens (log2FoldChange < 0).
Once we subset the data we want, we extract the gene names and store as an object.
```
gonads_degs_queens <- subset(x = deseq_object_results_sig_gonads,
padj < 0.05 & log2FoldChange < 0)
gonads_degs_queens_list <- row.names(gonads_degs_queens)
```

#### 2) Brains:
task: repeat the above for brains
##### a) drones:
##### b) queens:


##### Check whether caste-biased genes are conserved across tissues.
```
euler_drones_lists <- euler(combinations = list("gonads - drones" = gonads_degs_drones_list,
"Brains - drones" = brains_degs_drones_list))
euler_queens_lists <- euler(combinations = list("gonads - queens" = gonads_degs_queens_list,
"Brains - queens" = brains_degs_queens_list))
```
##### Visualise overlap for genes with drone-biased expression:
```
euler_drones_plot <- plot(euler_drones_lists,
quantities = TRUE,
edges = TRUE,
labels = list(fontsize = 8))
## Print to console:
euler_drones_plot
```

## Save plot to output:
```
ggsave("results/euler_drones_plot.png",
plot = euler_drones_plot,
height = 5,
width = 5)
```

##### TASK:
Visualise overlap for genes with queen-biased expression:
