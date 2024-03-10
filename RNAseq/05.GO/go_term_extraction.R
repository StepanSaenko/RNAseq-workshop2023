#!/usr/bin/env Rscript
#
#title: "extract go terms from ensembl metazoa"
#output: extract_ensembl.html
#author: Joe Colgan (joscolgan)


## Introduction:
#This script takes a tab-delimited file containing two columns:
#- Column 1 - Query (QuerySpecies) protein ID
#- Column 2 - Ensembl Reference Species eg Apis (RefSpecies) protein ID for predicted homolog
#The script uses biomaRt to obtained Gene Ontology terms assigned to homologs in Apis mellifera.
#The script matches Euglossa homologs to Gene Ontology assignaed to Apis homologs and writes a dataframe to file.
#The script performs the same analysis using Drosophila melanogaster.

#1. Load libraries:

require(biomaRt)
require(dplyr)
require(tidyr)
#inputfile = "data/EuglossaDilemma_ApisMellifera_1v1_orthologs.sorted.tsv"
inputfile = "XXXXXX"    #replace with path to file (eg on commandline via sed)


#2. Load ensembl database:


## replace ##XXXX to uncomment your required version
##AMEL ensembl <- useEnsemblGenomes(biomart = "metazoa_mart", dataset = "amellifera_eg_gene")
##BTER ensembl <- useEnsemblGenomes(biomart = "metazoa_mart", dataset = "bterrestris_eg_gene")
##BIMP ensembl <- useEnsemblGenomes(biomart = "metazoa_mart", dataset = "bimpatiens_eg_gene")
##NVIT ensembl <- useEnsemblGenomes(biomart = "metazoa_mart", dataset = "nvitripennis_eg_gene")                            


#3. Extract Gene Ontology terms directly assigned to RefSpecies (eg Apis):


## Extract gene coordinates:
ref_ncbi <- getBM(attributes = c('ensembl_gene_id',
                                  	'ensembl_peptide_id',
                                   	'entrezgene_id',
                                   	'chromosome_name',
                                   	'start_position',
                                   	'end_position',
                                  	'go_id',
                                  	'name_1006'),
                    				 mart = ensembl)



#4.Read in table containing list of QuerySpecies and RefSpecies homologs: 
 


query_data <- read.table(file = inputfile, header = FALSE,
                            col.names = c("query_prot", "ref_prot"))

## Remove extension from protein ID:
query_data$ref_prot <- gsub(pattern = "[.]1", replacement = "",
                                query_data$ref_prot)


#5. Subset ensembl gene id, GO id and GO term description:


## Match RefSpecies terms in both dataframes and add associated QuerySpecies protein ids to original dataframe:
## Match Drosophila gene id in both dataframes and add RefSpecies gene ID to dataframe containing GO terms:

ref_ncbi$query_prot <- query_data[match(ref_ncbi$ensembl_peptide_id,
                                               query_data$ref_prot), ]$query_prot

query_ref_go_terms <- ref_ncbi %>%
        select(c(query_prot, go_id))

## Remove blank lines:
query_ref_go_terms <- subset(query_ref_go_terms,
                            query_prot != "")

query_ref_go_terms <- subset(query_ref_go_terms,
                            go_id != "")


#6. Write to file:


## Create output directory:
# dir.create(path = "results")

## Write to file:
##AMEL write.table(query_ref_go_terms, file = paste(inputfile,'.GO.Amel.raw.txt', sep = ''), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
##BTER write.table(query_ref_go_terms, file = paste(inputfile,'.GO.Bter.raw.txt', sep = ''), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
##BIMP write.table(query_ref_go_terms, file = paste(inputfile,'.GO.Bimp.raw.txt', sep = ''), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
##NVIT write.table(query_ref_go_terms, file = paste(inputfile,'.GO.Nvit.raw.txt', sep = ''), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


#7. To get RefSpecies genes assigned with GO terms for Drosophila orthologs:


dros_ensembl <- useEnsemblGenomes(biomart = "metazoa_mart", dataset = "dmelanogaster_eg_gene")


#8. Extract Gene Ontology terms for Drosophila as well as the :


## Check attributes:
listAttributes(dros_ensembl)[521:530, ]

## Extract gene coordinates:
## Biomart doesnt allow to directly access GO terms and homolog information unfortunately:

dros_go_terms <- getBM(attributes = c('ensembl_gene_id', 'go_id', 'name_1006'), mart = dros_ensembl)

## Independetly extract Drosophila gene ID and RefSpecies homolog:
##AMEL dros_ref_homologs <- getBM(attributes = c('ensembl_gene_id', 'amellifera_eg_homolog_ensembl_gene', 'amellifera_eg_homolog_ensembl_peptide'), mart = dros_ensembl)
##BTER dros_ref_homologs <- getBM(attributes = c('ensembl_gene_id', 'bterrestris_eg_homolog_ensembl_gene', 'bterrestris_eg_homolog_ensembl_peptide'), mart = dros_ensembl)
##BIMP dros_ref_homologs <- getBM(attributes = c('ensembl_gene_id', 'bimpatiens_eg_homolog_ensembl_gene', 'bimpatiens_eg_homolog_ensembl_peptide'), mart = dros_ensembl)
##NVIT dros_ref_homologs <- getBM(attributes = c('ensembl_gene_id', 'nvitripennis_eg_homolog_ensembl_gene', 'nvitripennis_eg_homolog_ensembl_peptide'),  mart = dros_ensembl)

## Match Drosophila gene id in both dataframes and add Apis mellifera gene ID to dataframe containing GO terms:
##AMEL dros_go_terms$amellifera_eg_homolog_ensembl_peptide <- dros_ref_homologs[match(dros_go_terms$ensembl_gene_id, dros_ref_homologs$ensembl_gene_id), ]$amellifera_eg_homolog_ensembl_peptide
##BTER dros_go_terms$bterrestris_eg_homolog_ensembl_peptide <- dros_ref_homologs[match(dros_go_terms$ensembl_gene_id, dros_ref_homologs$ensembl_gene_id), ]$bterrestris_eg_homolog_ensembl_peptide
##BIMP dros_go_terms$bimpatiens_eg_homolog_ensembl_peptide <- dros_ref_homologs[match(dros_go_terms$ensembl_gene_id, dros_ref_homologs$ensembl_gene_id), ]$bimpatiens_eg_homolog_ensembl_peptide
##NVIT dros_go_terms$nvitripennis_eg_homolog_ensembl_peptide <- dros_ref_homologs[match(dros_go_terms$ensembl_gene_id, dros_ref_homologs$ensembl_gene_id), ]$nvitripennis_eg_homolog_ensembl_peptide


##AMEL dros_go_terms$query_prot <- query_data[match(dros_go_terms$amellifera_eg_homolog_ensembl_peptide, query_data$ref_prot), ]$query_prot
##BTER dros_go_terms$query_prot <- query_data[match(dros_go_terms$bterrestris_eg_homolog_ensembl_peptide, query_data$ref_prot), ]$query_prot
##BIMP dros_go_terms$query_prot <- query_data[match(dros_go_terms$bimpatiens_eg_homolog_ensembl_peptide, query_data$ref_prot), ]$query_prot
##NVIT dros_go_terms$query_prot <- query_data[match(dros_go_terms$nvitripennis_eg_homolog_ensembl_peptide,  query_data$ref_prot), ]$query_prot



## Rearrange order of Gene Ontology terms:
query_dros_go_terms <- dros_go_terms %>%
        select(c(query_prot, go_id))

## Remove blank lines:
query_dros_go_terms <- subset(query_dros_go_terms,
                             query_prot != "")
query_dros_go_terms <- subset(query_dros_go_terms,
                             go_id != "")


#9. Write to file:

## Create output directory:
# dir.create(path = "results")

## Write to file:
##AMEL write.table(query_dros_go_terms, file = paste(inputfile,'.GO.Amel-Dmel.raw.txt', sep = ''), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
##BTER write.table(query_dros_go_terms, file = paste(inputfile,'.GO.Bter-Dmel.raw.txt', sep = ''), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
##BIMP write.table(query_dros_go_terms, file = paste(inputfile,'.GO.Bimp-Dmel.raw.txt', sep = ''), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
##NVIT write.table(query_dros_go_terms, file = paste(inputfile,'.GO.Nvit-Dmel.raw.txt', sep = ''), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
