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

#inspect the results

head(results(deseq_object, tidy=TRUE))

summary(res)

res <- res[order(res$padj),]
head(res)
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
between queens and drones.  

```
ggplot(data = deseq_object_results_df) +
  geom_point(aes(x = log2FoldChange,
                 y = padj),
             alpha = 0.5) +
  geom_point(data = subset(x = deseq_object_results_df,
                           padj < 0.01),
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
                           padj < 0.01),
             aes(x = log2FoldChange,
                 y = -log10(padj)), #Add -log10() here
             colour = "red")
```

Cool - now it looks like a volcano ... plot!    

Lastly, we can store our final plot as an object and save the output to file.
We will colour all datapoint red which have an adjusted p-value of <0.01 and log2foldchange>2

```
volcano_plot <- ggplot(data = deseq_object_results_df) +
  geom_point(aes(x = log2FoldChange,
                 y = -log10(padj)),
             alpha = 0.5) +
  geom_point(data = subset(x = deseq_object_results_df,
                           padj < 0.01 & abs(log2FoldChange)>2),
             aes(x = log2FoldChange,
                 y = -log10(padj)),
             colour = "red")


## Save image:
ggsave(filename = "results/volcano_plot_brains.png",
       plot = volcano_plot,
       height = 5,
       width = 5)
```
## 8. extract and save the list of genes with sign. diff expression and/or log2fold cange
```
genes <- subset(x = deseq_object_results_df, padj < 0.01 & abs(log2FoldChange)>2)
write.table(genes, file='results/genes.brains.tsv', quote=FALSE, sep='\t', col.names = NA)

```

## 9. optional if data are available: Examine expression of candidate genes from the qPCR-based analysis:  
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
