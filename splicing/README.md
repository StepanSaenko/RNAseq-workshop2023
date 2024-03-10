# Differential splicing analysis  
## 1. Introduction  
This tutorial examines differential or alternative splicing patterns between male and queen
honeybees (_Apis mellifera_).  
The original data for the tutorial was obtained from the NCBI BioProject database [PRJNA689223](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA689223).  
The study that described the generation of the data can be found [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9156628/).  

In brief, the dataset consists of brains and gonads sampled from honeybee queens and males (drones).  

The primary aims of the tutorial is to:
- Investigate alternatively spliced transcripts between queens and drones  
- Investigate overlap and differences in splicing patterns across two tissues  
- Examine overlap between alternatively spliced and genes that differ in expression amplitude between queens and drones    

__Expected learning outcomes:__  
By the end of the tutorial, you should be:
- Familiar with different splicing events  
- Criteria by which to identify and determine presence of alternative splicing  
- Differentiate between splicing and differences in terms of amplitude   

## 2. Background to data processing for differential splicing analysis  
For each sample, data were downloaded and extracted from NCBI using the [sratoolkit](https://github.com/ncbi/sra-tools).  
The data were quality assessed using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) before each sample
was aligned using [STAR](https://github.com/alexdobin/STAR) against the latest available reference genome assembly, which obtained from [Ensembl Metazoa](https://metazoa.ensembl.org/Apis_mellifera/Info/Index). Using the alignment BAM files, we then ran [rMATS](https://rnaseq-mats.sourceforge.io/download.html).  

## 3. Analysing the output data from rMATs  

