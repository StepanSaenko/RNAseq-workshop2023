# this part will create a heatmap with single genes as rows (all samples as columns)
```

library(dplyr)
library(gplots)
library(RColorBrewer)

´´´
### creating a dataframe from abundances of all samples
### thereby selecting only those genes that are differentially expressed between queen / worker brains
```
gene_names <- rownames(genes)
DF = data.frame(txi_counts$counts)
DF.selected <- DF %>% filter(row.names(DF) %in% gene_names)
´´´

### input <- read.delim('abundance_table.txt', sep='\t', header=TRUE)  ### table with normalised readcounts / DESEq2 output, genes diff. expressen, p >0.01 o.ä.
```
input <- DF.selected
mat_input <- data.matrix(input)
rownames(mat_input) <- rownames(input)
´´´
### because read count vary a lot between genes, log2 transformation does compensation of extreme values
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

png(file='heatmap.png', res=240, height=2000, width=2000)  
heatmap.2(linput, dendrogram=dendrogramtoplot, Colv=reorder_cols, Rowv=reorder_rows,
    distfun=dist_fun, hclustfun=hclust_fun, scale = scale, labRow = rlabs, labCol = clabs,
    col=colused, trace="none", density.info = "none", margins=label_margins,
    main = '', cexCol=0.8, cexRow=0.8, srtCol=srtCol,
    keysize=3, key.xlab='', key.title='', key.par=key_margins,
    lmat=layout_matrix, lhei=lheight, lwid=lwidth)
 
dev.off()
´´´
