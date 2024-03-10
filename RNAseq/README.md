#RNAseq analysis with kallisto and DeSeq2
authors: Eckart Stolle (LIB Bonn) & Joe Colgan (U Mainz)
Transcriptomic analysis of Apis mellifera

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


## 2. Load libraries:
Next, we load packages (also called _libraries_) that can be used by R to analyse your data.  
These are packages that contain *functions* written by other researchers in R 
that allow us to perform tasks without having to create or write functions ourselves.  

For this first script, we will load five packages that help with formatting data, 
and visualising plots.  

To load the required packages, we use the 'library()' function:
```
# load packages
library(tximport)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
```

## 3. Load data:
Once we have libaries loaded, we can load our data. 
Our input data comes from the software tool `salmon` or `kallisto`, which performs transcriptome-based
mappings and quantifies transcript expression (this has already been performed) for
each sample. The main output of `salmon/kallisto` that we will use for our transcriptomic
analysis is the `quant.sf` (salmon) or `abundace.tsv` (kallisto) file with a file being generated for each sample.  
Each `quant.sf` file contains information on:  
- transcript id  
- transcript length   
- effective length  
- transcripts per million (TPM)
- number of reads

Kallisto's `abundance.tsv` file has a similar structure
- target_id
- length
- eff_length
- est_counts
- tpm

For the present analysis, TPM values will be used.  

We can load the data in using the following steps - it may look a bit complicated
but this is the most complicated part of the tutorial!  

```
## check you working DIR
getwd()

## get list of PATHs/Files (one output file per Sample
## Set paths to folders containing output files from salmon/kallisto
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
# names(paths$count_files) <- paste0("sample", 1:17) #standard sample names sample1 ... sample17
names(paths$count_files) <- gsub(paths$count_files_relative,
                                    pattern = "/.*",
                                    replacement = "")

# Extract sample names and put into df for tximport
samples     <- data.frame(treatment = names(paths$count_files))
```

## 4. Load sample information data:  
For our analysis, researchers that studied the honeybees in our experiment, collected information about the samples, which is stored in a file `sample.info.short.tsv`.  This file contains information on whether an individual sample is a queen or drone (in the table referred to as `caste`, although this is not really correct as we are talking about different sexes), which is important as we want to ask whether they differ in gene expression. In addition, the file also contains information on whether the sample originated from brain or gonad tissue while also what age the bees had. 
```
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
```

## 5. Generate gene-level counts:
To perform our analysis, first, we read in a file `transcript.gene.map.txt` and store it as an object 'transcript_to_gene'. 
The file looks like the example below, where column1 is the transcript ID and column 2 is the gene/locus ID. We can extract such information from the GTF file
```
XR_001705491.2  LOC551580
XR_001705490.2  LOC551580
```

```
# load transcripts to gene mapping
transcript_to_gene <- read.table(file = "transcript.gene.map.txt",
                                 col.names = c("transcript",
                                               "locus"))
```

Next, we use the function 'tximport()' from the R package 'tximport' to match our transcripts to individual genes. This produces summarised counts for each gene. If we did not do this, our analysis would compare individual transcripts (RNA) but we want to compare expression at the gene level. Therefore, for each gene, we sum the counts for their corresponding transcripts.    
```
# import data
txi_counts <- tximport(paths$count_files,
                       type    = "kallisto",
                       tx2gene = transcript_to_gene,
                       countsFromAbundance = "no")

head(txi_counts$counts)
```

Next, we construct a DESeq2 object which stores:  
- Gene-level count data  
- Sample information  
- Design  
With the 'design' argument, we can specify variables that we want DESeq2 to be aware for. For example, this information includes:  
- 'age': the age that the bee had
- 'caste': whether the individual is a queen or a drone.  
- 'tissue': whether the tissue was from brains or gonads.  

We then create the object using the function 'DESeqDataSetFromTximport'.  

```
#create Deseq object(s)
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
```
We observe an Error here:
```
Error in checkFullRank(modelMatrix) :   the model matrix is not full rank, so the model cannot be fit as specified.
One or more variables or interaction terms in the design formula are linear
combinations of the others and must be removed.
```

Look at the sample information again: `print(samples_information)`. Age is not ressolved by caste, queens are all 7 and drones are all 15 days old. Hence we cannot discriminate age from caste/sex here. We will ignore the age since development is faster in queens than drones and hence likely the sampled age reflect similar developmental stages. However, we cannot test this assumption here.

```
deseq_txi <- DESeqDataSetFromTximport(txi     = txi_counts,
                                      colData = samples_information,
                                      design  = ~caste + tissue)
```

Using this DESeq2 object, we can check how many genes our species has. This can be performed using the 'counts()' function. We can then _nest_ this function inside of another called 'nrow()', a base R function, which will count the *n*umber of *row*s in our dataframe.  

```
## how many genes does our psecies have?
nrow(counts(deseq_txi))
```

How many genes does our species have?  

Now, a potential issue with any transcriptomic analysis is lowly expressed genes or genes that are not expressed in any samples. These are not biologically informative to answering our main question and therefore, we can remove them.  

In the next steps, we remove genes where the total number of reads across all samples is 10, which is quite low. Therefore, for a gene to be kept in our analysis, it must have more than 10 reads across our 17 samples.

```
## Filter out lowly expressed genes:
keep      <- rowSums(counts(deseq_txi)) >= 10
deseq_txi <- deseq_txi[keep, ]
print(nrow(deseq_txi))
```

Last, for this section, we save our DESeq2 object as we can use it again in a later script. To save our object, we use the funtion 'saveRDS()', which is a base R function that saves our object to an R Data Serialisation (RDS) file. These files can store objects and can later be loaded into other R scripts and used there. 
```
## Create output 'results' directory:
dir.create(path = "results")
saveRDS(object = deseq_txi,
        file = "results/deseq_txi.rds")
```

## 6. Generate a DESeq object:
To allow for running statistical tests, we run the DESeq2 which implements generalised linear models to examine relationships between variables we measured and the expression of each gene. However, to perform more basic analyses, such as clustering-based analyses, we need to first run the function 'DESeq()', which creates a DESeq2-formatted object.    
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
Using the normalised counts stored in the 'vsd' object, we can perform a principal component analysis (PCA). PCAs help to visualise the overall similarities and differences among samples allowing for identifying outlier samples or batch effects (e.g., whether individuals from the same colony or sequencing run are more similar to each other than caste or tissues).  

To perform a PCA, we run the function 'plotPCA()', which performs both the PCA but also plots the output of the analysis. By default, it outputs the first two principal components, which each explain a certain amount of variance in our dataset. 
```
## PCA
pcaplot <- plotPCA(vsd,
        intgroup = c("tissue",
                     "caste")) +
  theme_bw()

jpeg(file="results/PCA.jpeg")
pcaplot
dev.off()
```

## 9. Examine expresison profiles between tissues:  
Another use of the transformed data is sample clustering. Here, we apply the 'dist()' function to calculate sample-to-sample distances. The greater the differences in gene expression between two samples, the greater the distance.
```
sampleDists <- dist(t(assay(vsd)))
```
A heatmap of this distance matrix gives us an overview of similarities and dissimilarities between samples.
```
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$submorph,
                                    vsd$tissue,
                                    sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9,
                                          "Blues")))(255)  
```
Lastly, we plot and save the heatmap:
```
heatmap <- pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

jpeg(file="results/Heatmap.jpeg")
heatmap
dev.off()
```

Based on the profiles above, we can see that the two tissues (brains, gonads) differ considerably 
in terms of overall gene expression (more than the caste/sex). Therefore, to understand differences between 
queens and drones, we should analyse each tissue independently.  

Therefore, we can subset the two tissues:

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
```


# Part 2. tissue-specific analysis - brains
## 1. Load libraries:
As with the first script, the first thing we want to do is load libraries that contain functions that we will use to:
- Perform a principal component analysis  
- Identify differentially expressed genes between queen and drones  
- Visualise differentially expressed genes  

```
library(tidyverse)
library(DESeq2)
```

## 2. Load data:
Here, we will load an object from our previous script 'deseq2_analysis_part1_tissues.Rmd' that contains the DESeq2 object, which includes information on gene-level counts, sample information, as well as the design.  

We can load the RDS object using the function 'readRDS()' and store the information as an object.  

```
deseq_txi_brains <- readRDS(file = "results/deseq_txi_brains.rds")
```

## 3. Perform a differential expression analysis:  
We will use a statistical test called a likelihood ratio test (LRT).
To perform an LRT, we create two models - a full model containing all of our
variables, as well as a reduced model, which contains the variable that we may not specifically be interested in but we want the model to control for. In this case, we have the following:
- full model: age + caste
- reduced model: caste
- BUT: age is not variable among the castes, hence we can go only for the rduced model here

The first thing that we want to change is the 'design' as we no longer have 
'tissue' as an explanatory variable. We can readjust the 'design' of our object using the 'design()' function. 

In addition, we want to set 'queens' as our reference, which we can do 
using the 'relevel()' function.  

```
design(deseq_txi_brains) <- ~caste

deseq_txi_brains$caste <- relevel(deseq_txi_brains$caste,
                                       ref = "queen")

deseq_object <- DESeq(deseq_txi_brains,
                      test = "LRT",
                      reduced = ~1)
```

## 4. Perform a principal component analysis:  

Similar to the tissue-based analysis, we want to first examine how samples cluster together for our tissue of interest. For example, we want to check whether queens are more similar to each other and drones more similar to each other, respectively.

First, like in the last script, we use the 'vst()' function, which performs variance stabilising transformation, a form of normalisation of count data.

```
vsd_brains <- vst(deseq_object,
                  blind = FALSE)
```

Next, using the normalised counts as input, we perform a principal component analysis and plot the results. Both steps are performed using the function 'plotPCA()', which takes our normalised count data ('vsd_brains') as input. We can also tell the function, which information from our samples that we want to visualise, using the argument 'intgroup'. Here, we provide 'caste', which will indicate on the plot, which samples are queens and which are drones.

Once we plot, we can then save the plot as a '.PNG' file.   

```
pca_plot <- plotPCA(vsd_brains,
        intgroup = c("caste")) +
  theme_bw()

## Save pca_plot to file:
ggsave(filename = "results/pca_plot_brains.png",
       plot = pca_plot,
       height = 5,
       width = 5)
```

## 5. Complete differential expression analysis:  
Finding out which genes are differentially expressed:
Using the 'resultsNames()' function, we can see that DESeq2 provides results on differential expression for caste , or (if more variable would have been used such as age then both variables (age and caste), which 
were included in our experimental design. As we are only interested in differences between queens and drones (i.e., caste), we can ignore the other comparisons.

```
resultsNames(deseq_object)
```

To extract results from our DESeq2 object, we use the function 'results()' and store all comparisons performed by DESeq2 in an object called 'res'.  

```
res <- results(deseq_object)
```

We can also specifically extract results relevant to answering our scientific question by using the argument 'name' inside of the function 'results'.  

```
deseq_object_results <-  results(deseq_object,
                                 name = "caste_drone_vs_queen")
```

To identify genes that are significantly differentially expressed, we use the 'subset()' function, which allows us to subset a dataframe based on a condition. Here, we ask to subset the column 'padj', which contains information on the adjusted p-value, which we use to identify the signifiance of a gene. For example, if a gene has an adjusted p value < 0.05, we consider that it significantly differs in terms of expression between queens and drones.  

Below, we use the 'dim()' function, which counts the number of rows and columns that a dataframe has.  

In addition, we save this object containing information on significantly differentially expressed genes as we will use it again in a later script.  

```
deseq_object_results_sig <- subset(x = deseq_object_results,
                                   padj < 0.05)
dim(deseq_object_results_sig)

## Let us save this file and we can use it again in a later script:
saveRDS(object = deseq_object_results_sig,
        file = "results/deseq_object_results_sig_brains.rds")
```

## 6. Identifying direction of expression:  
We now know how many genes are differentially expressed in total but how many have elevated expression in queens and how many in drones?

To answer this question, we can use another metric calculated by DESeq2, which is log2FoldChange. Log2FoldChange is a value that indicates how much the expression of gene (or transcript) differs between two groups. The value is reported on a logarithmic scale to base 2 meaning a log2FoldChange = 1 means the expression of a gene is twice as high in one compared to an another. In addition, whether this value is positive or negative can help identify the direction of expression.

For example, in our analysis, queens are the reference level meaning that genes with elevated log2FoldChange are higher in drones compared to queens while a negative value means expression is higher in queens compared to drones.

To identify these patterns, we can use the 'subset()' function again but instead, here, we subset based on the 'log2FoldChange' column. To count the number of genes, we can count the number of rows in the subsetted dataframe using the function 'nrow()'.

```
nrow(subset(x = deseq_object_results_sig,
            log2FoldChange > 0))

nrow(subset(x = deseq_object_results_sig,
            log2FoldChange < 0))
```

## 7. Visualise differentially expressed genes:  
Lastly, we want to visualise our differentially expressed genes. 
Here, we want to generate a Volcano plot which is a scatterplot where we plot log2FoldChange on the x-axis (so it helps see which genes have a more nurse (negative values) or more forager-biased expression (positive values)). On the y-axis, we want to plot our significance value so we can see how much our genes of interest differ in expression between queens and drones. Therefore, we want to plot adjusted p-value, which we use for defining whether a gene is significantly differently expressed or not. 

However, we know that the adjusted p-values range from 0 - 1 and values can get really, really, really small so instead, we change to scale by using the '-log10()' function which makes the values easier to read when plotting them. 

```
## Convert the output of DESeq2 to a dataframe using the 'as.data.frame()' function:  
deseq_object_results_df<- as.data.frame(deseq_object_results)
```

For plotting, we will use ggplot2, an R package with functions that allow 
for making lots of different type of plots, such as:  
- scatterplots
- histograms
- density plots
- volcano plots
- etc.
It is one of the most commonly used tools and you will see many plots in 
scientific publications that come from using the functions within ggplot2.  

So how do we make a ggplot plot:

First, we need a base - what happens when we run the following command?

```
ggplot(data = deseq_object_results_df)
```
Here, ggplot creates the background / canvas for our plot - next, we need to add data points.  

Next, we add the geom_point() function, which is specificially used for plotting points like
in a scatterplot. Within the 'geom_point()' function, we also use the 'aes()' (aesthetic) 
subfunction, where we can specify what we want plotted on the x and y axes.  

Here, below, we want to plot log2FoldChange on the x-axis and the measure of
significance (adjusted p value' padj) on the y-axis.  

```
ggplot(data = deseq_object_results_df) +
  geom_point(aes(x = log2FoldChange,
                 y = padj))
```

That is a lot of data - we can make the dots more transparent using the 
'alpha' argument of the 'geom_point()' function. By default, the alpha value is 
1, which makes the point opaque (you can't see through it) while an alpha value
of 0 makes the point invisible. Here, we will set it 0.5 so the points are
more transparent.  

```
ggplot(data = deseq_object_results_df) +
  geom_point(aes(x = log2FoldChange,
                 y = padj),
             alpha = 0.5)
```

We can also indicate which genes are identified as significantly differentially
expressed from our DESeq2 analysis. To do this, we can add a second geom_point()
function, which adds a layer on top of the original points. 

We can also add the argument 'colour' to set the colour of specific points. Here, 
red dots will indicate points that are significantly differentially expressed
between nurses and foragers.  

```
ggplot(data = deseq_object_results_df) +
  geom_point(aes(x = log2FoldChange,
                 y = padj),
             alpha = 0.5) +
  geom_point(data = subset(x = deseq_object_results_df,
                           padj < 0.05),
             aes(x = log2FoldChange,
                 y = padj),
             colour = "red")
```

Lastly, the y-axis can be a bit difficult to read - to make it clearer, we can
change the scale of the y-axis values using the function '-log10()', which 
transforms the data onto the logarithmic scale. 

```
ggplot(data = deseq_object_results_df) +
  xlim(-20, 10) +
  geom_point(aes(x = log2FoldChange,
                 y = -log10(padj)), #Add -log10() here
             alpha = 0.5) +
  geom_point(data = subset(x = deseq_object_results_df,
                           padj < 0.05),
             aes(x = log2FoldChange,
                 y = -log10(padj)), #Add -log10() here
             colour = "red")
```

Cool - now it looks like a volcano ... plot!    

Lastly, we can store our final plot as an object and save the output to file.  

```
volcano_plot <- ggplot(data = deseq_object_results_df) +
  geom_point(aes(x = log2FoldChange,
                 y = -log10(padj)),
             alpha = 0.5) +
  geom_point(data = subset(x = deseq_object_results_df,
                           padj < 0.05),
             aes(x = log2FoldChange,
                 y = -log10(padj)),
             colour = "red")

## Save image:
ggsave(filename = "results/volcano_plot_brains.png",
       plot = volcano_plot,
       height = 5,
       width = 5)
```


## 8. optional if data are available: Examine expression of candidate genes from the qPCR-based analysis:  
The last thing we want to do with the data are check the expression of genes that
we examined in the qPCR-based practical.  
Here, we have a list of genes from the NCBI gene database for _vitellogenin_ genes
in the _Apis_ genome assembly.  

- Vitellogenin  (n copies in the _Apis_ genome)  
LOCxxxx1  
LOCxxxx2  
LOCxxxx3  
LOCxxxx4  
LOCxxxx5  

Store these gene names in an object called 'vitellogenin_genes' and then subset
our dataframe containing differentially expresssed genes - the output of which
 is stored in a new object called 'vitellogenin_degs'.  

```
vitellogenin_genes <- c("LOCxxxx1",
                        "LOCxxxx2",
                        "LOCxxxx3",
                        "LOCxxxx4",
                        "LOCxxxx5")

## Subset:
vitellogenin_degs <- subset(deseq_object_results_df,
          row.names(deseq_object_results_df) %in% vitellogenin_genes)

## Are any genes differentially expressed?
vitellogenin_degs
```

Similarly, we can do the same with _PGRP_ genes.  

- PGRPs (twenty copies in the _Apis_ genome):   
LOCxxxx6  
LOCxxxx7  
LOCxxxx8   

Store these gene names in an object called 'pgrp_genes' and then subset
our dataframe containing differentially expresssed genes - the output of which
 is stored in a new object called 'pgrp_degs'.  

```{r, message = FALSE, results = 'hide', warning = FALSE}
pgrp_genes <- c("LOCxxxx6",
                "LOCxxxx7",
                "LOCxxxx8")

## Subset:
pgrp_degs <- subset(deseq_object_results_df,
          row.names(deseq_object_results_df) %in% pgrp_genes)

## Are any genes differentially expressed?
pgrp_degs
```


# Part 3. tissue-specific analysis - gonads

Task: repreat the above analysis for gonads


# Part 4. overlap
The purpose of this script is to examine overlap in genes that are differentially expressed between nurses and
foragers in the antennae and the brains.
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
deseq_object_results_sig_antennae <- readRDS(
file = "results/deseq_object_results_sig_gonads.rds")
deseq_object_results_sig_brains <- readRDS(
file = "results/deseq_object_results_sig_brains.rds")
```

## 3. Extract significantly differentially expressed genes:
As mentioned above, the row names of our DESeq2 output contain the gene name, which we can extract
using the function ‘rownames()’ and store as an object for each individual tissue.
```
antennae_degs <- rownames(deseq_object_results_sig_gonads)
brains_degs <- rownames(deseq_object_results_sig_brains)
```

## 4. Combine input genes lists to examine overlap:
Next, we use the ‘euler()’ function from the eulerr R package, which creates a list using the DEGs found in each tissue.
```
euler_lists <- euler(combinations = list("Antennae" = antennae_degs,
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
```

6. Examine changes per caste:
Next, we want to look at specific overlap within a submorph.

1) Gonads
a) drones: For foragedronesrs, we subset the output of DESeq2 to obtain genes that are both significantly expressed (padj < 0.05) and have a bias in terms of expression in drones (log2FoldChange > 0).
Once we subset the data we want, we extract the gene names and store as an object.
```
gonads_degs_f <- subset(x = deseq_object_results_sig_gonads,
padj < 0.05 & log2FoldChange > 0) gonads_degs_f_list <- row.names(gonads_degs_f)
```

b) queens:
Same with nurses, we subset the output of DESeq2 to obtain genes that are both significantly expressed
(padj < 0.05) and have a bias in terms of expression in queens (log2FoldChange < 0).
Once we subset the data we want, we extract the gene names and store as an object.
```
gonads_degs_n <- subset(x = deseq_object_results_sig_gonads,
padj < 0.05 &
log2FoldChange < 0)
gonads_degs_n_list <- row.names(gonads_degs_n)
```

2) Brains:
task: repeat the above for brains




Check whether caste-biased genes are conserved across tissues.
```
euler_f_lists <- euler(combinations = list("gonads - drones" = gonads_degs_f_list,
"Brains - F" = brains_degs_f_list))
euler_n_lists <- euler(combinations = list("gonads - queens" = gonads_degs_n_list,
"Brains - N" = brains_degs_n_list))
```
Visualise overlap for genes with drone-biased expression:
```
euler_f_plot <- plot(euler_f_lists,
quantities = TRUE,
edges = TRUE,
labels = list(fontsize = 8))
## Print to console:
euler_f_plot
```

## Save plot to output:
```
ggsave("results/euler_drone_plot.png",
plot = euler_f_plot,
height = 5,
width = 5)
```

TASK:
Visualise overlap for genes with queen-biased expression:
