## This part will create a heatmap with single genes as rows (all samples as columns)
needed libraries (dplyr for the filterig step, gplots includes dendrogramtoplot / heatmap2):
```
library(dplyr)
library(gplots)
library(RColorBrewer)
```

creating a dataframe from abundances of all samples
e.g. if you have an external table with normalised readcounts, selected for e.g. p>0.05 you can load it with:
input <- read.delim('abundance_table.txt', sep='\t', header=TRUE)

We will use our data from previous steps:
- Part1 step 5 created "txi_counts"
- Part2, step 8 created "genes" (with p>0.01 and at least 2-fold change)
  
thus, we create an abundance table fro those genes that are significantly differentially expressed
```
gene_names <- rownames(genes)
DF = data.frame(txi_counts$counts)
DF.selected <- DF %>% filter(row.names(DF) %in% gene_names)

input <- DF.selected
mat_input <- data.matrix(input)
rownames(mat_input) <- rownames(input)
```
because read count vary a lot between genes, log2 transformation does compensation of extreme values
```
    linput <- log2(mat_input)
    scale <- "none"
    linput[linput=="-Inf"] <- 0

    srtCol <- 30
    rlabs <- NULL
    clabs <- NULL
    label_margins <- c(8,8)

    dendrogramtoplot <- "both"
        reorder_cols <- TRUE   
        reorder_rows <- TRUE
        layout_matrix <- rbind(c(4,3), c(2,1))
        key_margins <- list(mar=c(4,0.5,2,1))
        lheight <- c(1, 5)
        lwidth <- c(1,3)
    hclust_fun <- function(x) hclust(x, method='complete')
        dist_fun <- function(x) dist(x, method='euclidean')

ncolors <- 11
    colused <- brewer.pal(ncolors, "BrBG")
```
if you dont want to save the figure as file you can omit the png() and dev.off() lines
```
png(file='results/heatmap_genewise.png', res=240, height=2000, width=2000)

heatmap.2(linput, dendrogram=dendrogramtoplot, Colv=reorder_cols, Rowv=reorder_rows,
    distfun=dist_fun, hclustfun=hclust_fun, scale = scale, labRow = rlabs, labCol = clabs,
    col=colused, trace="none", density.info = "none", margins=label_margins,
    main = '', cexCol=0.8, cexRow=0.8, srtCol=srtCol,
    keysize=3, key.xlab='', key.title='', key.par=key_margins,
    lmat=layout_matrix, lhei=lheight, lwid=lwidth)
 
dev.off()
```
