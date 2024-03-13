#!/usr/bin/env Rscript
##############################################################################
# Author: Joe Colgan                   Program: go_enrichment_analysis_fisher.R
# modified by Eckart Stolle, 2024
# Date: 16/08/2017
##############################################################################

setwd('/scratch/rnaseq2023/estolle/go')

# Load libraries; install from scratch if needed
libraries <- c("topGO",
               "lintr",
               "lattice",
		"Rgraphviz")
for (lib in libraries) {
    if (require(package = lib, character.only = TRUE)) {
        print("Successful")
    } else {
        print("Installing")
        source("https://bioconductor.org/biocLite.R")
        avebiocLite(pkgs = lib)
        library(lib, character.only = TRUE)
    }
}

## Step One: Input files, define output and objects for running topGO: 
# Load in genelist and database files.  

## Step One:
## Define input:
input <- "/scratch/rnaseq2023/estolle/go/Apis.mellifera.GOterms.tsv"

## GO annotations
gene_to_go_mapping_file <- input

## genes overlapping with regions of significantly low nucleotide diversity:
gene_input <- "gonads.txt"
gene_list <- scan(file = paste("genes_of_interest/",
                               gene_input,
                               sep = ""),
                  as.character())

## Define node size:
node_size <- 50

output_directory <- paste("results.topGO.gonads/go_term_output_",
                          node_size,
                          "_fisher",
                          sep = "")
if (file.exists(output_directory)) {
  stop("The output directory:", output_directory, ", already exists",
       "Let's avoid overwriting")
} else {
  dir.create(output_directory,
             recursive = TRUE)
}

## Read in input file:
## Read in GO annotations: 
gene_to_go_mapping <- readMappings(file = gene_to_go_mapping_file)
gene_universe <- names(gene_to_go_mapping)

genes_of_interest <- as.character(gene_list)

genelist <- factor(as.integer(gene_universe %in% genes_of_interest))
names(genelist) <- gene_universe

## Steps Two and Three: Create topGO Object & run tests for GO term enrichment
# We create a topGO object for each GO term
# We perform two statistical tests:
# 1. A ks test using the topGO 'weight01' algorithm
# 2. A Fisher's exact test using the topGO 'weight01' algoritm
#We combine the output of each test. 
#We filter out enriched terms.
#We do this for each of the three GO categories (ie. Biological process, Molcular Function, Cellular Component):

for (go_category in c("BP", "MF", "CC")) {
  # STEP TWO
  ## Build the GOdata object in topGO
  my_go_data <- new("topGOdata",
                    description = paste("GOtest", go_category, sep = "_"),
                    ontology    = go_category,
                    allGenes    = genelist,
                    gene2GO     = gene_to_go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = node_size) # Modify to reduce/increase stringency.
  # STEP THREE
  ## Calculate ks test using 'weight01' algorithm:
  result_weight_ks <- runTest(object    = my_go_data,
                              algorithm = "weight01",
                              statistic = "ks")
  ## Calculate fisher exact test using 'weight01' algorithm:
  result_weight_fisher <- runTest(object    = my_go_data,
                                  algorithm = "weight01",
                                  statistic = "fisher")
  ## Combine results from statistical tests:
  result_weight_output <- GenTable(object       = my_go_data,
                                   weight_ks     = result_weight_ks,
                                   weight_fisher = result_weight_fisher,
                                   orderBy       = "weight_fisher",
                                   topNodes      = length(score(result_weight_fisher)))
  ## Correct ks test for multiple testing:
  result_weight_output$weight_ks <- as.numeric(result_weight_output$weight_ks)
  result_weight_output$weight_fisher <- as.numeric(result_weight_output$weight_fisher)
  result_weight_output$weight_ks_adjusted <- p.adjust(p = result_weight_output$weight_ks,
                                                      method = c("BH"))
  result_weight_output$weight_fisher_adjusted <- p.adjust(p = result_weight_output$weight_fisher,
                                                          method = c("BH"))
  ## Write to output:
  write.table(x         = result_weight_output,
              file      = file.path(output_directory,
                                    paste(go_category,
                                          "raw.tsv",
                                          sep = "_")),
              row.names = FALSE,
              sep       = "\t",
              quote = FALSE)

   ## write sign terms
result_weight_output2 <- result_weight_output[which(result_weight_output$weight_fisher<=0.05),]
#0.001
#result_weight_output3[which(result_weight_output$weight_fisher_adjusted<=0.05),]
#save first top 50 ontolgies sorted by adjusted pvalues
  write.table(x         = result_weight_output2[1:50,],
              file      = file.path(output_directory,
                                    paste(go_category,
                                          "top50.tsv",
                                          sep = "_")),
              row.names = FALSE,
              sep       = "\t",
              quote = FALSE)



# PLOT the GO hierarchy plot: the enriched GO terms are colored in yellow/red according to significance level
#BiocManager::install("Rgraphviz")
pdf(file=file.path(output_directory,paste(go_category,"topGOPlot_fullnames.pdf",sep = "_")), height=12, width=12, paper='special', pointsize=18)
#showSigOfNodes(my_go_data, score(result_weight_output$weight_fisher), useInfo = "none", sigForAll=FALSE, firstSigNodes=2,.NO.CHAR=50))
#dev.off()

myterms = result_weight_output$GO.ID
mygenes = genesInTerm(my_go_data, myterms)

var=c()
for (i in 1:length(myterms))
{
myterm=myterms[i]
mygenesforterm= mygenes[myterm][[1]]
mygenesforterm=paste(mygenesforterm, collapse=',')
var[i]=paste("GOTerm",myterm,"genes-",mygenesforterm)
}
  write.table(x         = var,
              file      = file.path(output_directory,
                                    paste(go_category,
                                          "genetoGOmapping.txt",
                                          sep = "_")),
              row.names = FALSE,
              sep       = "\t",
              quote = FALSE)

}
