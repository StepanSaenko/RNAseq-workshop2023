# Part1
## 1. install packages
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
                tissue == "brain")$sample_id
deseq_txi_brains <- deseq_txi[, brain_samples]

## Check if the number subsetted is correct:
print(nrow(deseq_txi_brains))

## Save:
saveRDS(object = deseq_txi_brains,
        file = "results/deseq_txi_brains.rds")
```
