## GO term enrichment

### Steps
1. obtain a set of gene IDs you want to analyze
2. determine 1:1 orthologs between your species and a species available in ENSEMBL/BioMart (eg Apis mellifera, Bombus terrestris, Bombus impatiens, Nasionia vitripennis, Drosophila melanogaster ...). You can do this with the list of genes of interest or all genes of an annotation in that species. We can use broccoli for this task
3. fetch GO terms for the ENSEMBL/BioMART species and it's Drosophila orthologs
4. create genelist subsets
5. enrichment (TopGO)
6. visualitation (GOfigure)

### ATTENTION: for GO Figure we need to try if we quickly can get cona running
export PATH="/opt/anaconda3-2023.07-2-Linux-x86_64/bin:$PATH
conda init
conda activate /home/ek/.conda/envs/R1



## 3. create GO Term lists via BioMaRT
```
#install.packages('https://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_2.3.4.tar.gz', repos = NULL)

REFERENCEFOLDER="/scratch/rnaseq2023/Apis.mellifera.reference.genome"
REF="GCF_003254395.2_Amel_HAv3.1_genomic.fna"
GTF="GCF_003254395.2_Amel_HAv3.1_genomic.gtf"
GTF="GCF_003254395.2_Amel_HAv3.1_genomic.gff"
CDS="cds_from_genomic.fna"
TRANSCRIPTS="rna.fna"
PROTEIN="protein.faa"

###################################################
#make folder for GO term analyses and copy our results files over
mkdir go
cp results/genes.brains.tsv results/genes.gonads.tsv go/
cd go

###################################################
# make a list of Protein accessions from Amel Reference protein.faa
cat $REFERENCEFOLDER/$PROTEIN | grep ">" | tr -d '>' | cut -d" " -f1 > Apis.mellifera.proteinnames.txt

cat $REFERENCEFOLDER/$GTF | awk ' $3 == "transcript" ' | cut -f9 | cut -d";" -f1-6 | sed "s/ID=rna-//g" | sed "s/Parent=gene-//g" | sed "s/Dbxref=GeneID://g" | sed "s/gbkey=//g" | sed "s/gene=//g" | sed "s/,Genbank:/;/g" | sed "s/,BEEBASE:.*;Name=/;/g" | sed "s/Name=//g" | tr ';' '\t' > Apis.mellifera.transcript-LOC-geneID-data.txt
cat Apis.mellifera.transcript-LOC-geneID-data.txt | cut -f1 > Apis.mellifera.transcriptIDsfromGTF.txt

cat $REFERENCEFOLDER/$TRANSCRIPTS | grep ">" | tr -d '>' | cut -d" " -f1 > Apis.mellifera.transcriptnames.txt


###################################################
### GO term database
#1. Download GO terms from Enembl BioMart 
#https://metazoa.ensembl.org/biomart/martview/656f87cd558a5f254f02906f9621a2b2 
#awk '$3 > 0' mart_export\(1\).txt | tail -n +2 | cut -f 1,3 | sort | uniq | grep 'LOC' >  go_term_database_input.txt
#Rscript go_term_conversion.R go_term_database_input.txt go_term_database_output.txt

#2. GO terms from NCBI
cat $REFERENCEFOLDER/GCF_003254395.2_Amel_HAv3.1_gene_ontology.gaf | grep -v "^!" | cut -f3,5 > Apis.mellifera.ncbi.gene-ontologies.tsv
Rscript go_term_conversion.R Apis.mellifera.ncbi.gene-ontologies.tsv Apis.mellifera.ncbi.gene-ontologies.db.txt


###################################################

# create inputfile (query-protein_id ref_prot_id)
paste Apis.mellifera.proteinnames.txt Apis.mellifera.proteinnames.txt > Apis.mellifera.proteinnames2.txt

# fetch a generalized R script and set it to be used for Apis mellifera GO terms
INPUTFILE="$PWD/Apis.mellifera.proteinnames2.txt"

cat /home/ek/scripts/go_term_extraction.R | sed -r "s#XXXXXX#$INPUTFILE#g" | sed "s/##AMEL //g" > go_term_extraction.Amel.R
chmod 755 go_term_extraction.Amel.R

#Run R script
Rscript go_term_extraction.Amel.R

# if you get this error: "Error in `collect()`: ! Failed to collect lazy table" its related to an incompatibility between dbplyr and BiocFileCache, upgrade to newer BiocManager: BiocManager::install(version = "3.18")
# if your R is too old, you can either use conda to make an env wih newer R or try downgrade dbplyr: install.packages('https://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_2.3.4.tar.gz', repos = NULL)

#output when running with protein list as input (NCBI protein IDs)
head $INPUTFILE.GO.Amel-Dmel.raw.txt
#XP_016768898.1  GO:0046872
#XP_016768898.1  GO:0016020
#XP_016768898.1  GO:0016485

# Poblem: Protein IDs are not mapped to LOCUS IDs (genes)
# Above, we made a GO TErm database from Ensembl WebDOWNLOAD or NCBI REFSEQ (FTP): Apis.mellifera.ncbi.gene-ontologies.db.txt


# optional: run more ortholog lists or retreieve more GO term lists

# then merge possible multiple GO term lists
cat *.GO.*.raw.txt | grep -v "Uncharacterized" > go.raw.tsv

cp /home/ek/scripts/go_term_conversion.R .
chmod 755 go_term_conversion.R
Rscript ./go_term_conversion.R go.raw.tsv go.converted.tsv

head go.converted.tsv
head go.converted.tsv
#NP_001010975.1 GO:0005743,GO:0055085
#NP_001011563.1 GO:0005576,GO:0005615,GO:0019953

cat go.converted.tsv | wc -l
#13897 #more GO terms, but mapping Protein IDs to Locus IDs needs to be done for this to be used
cat Apis.mellifera.ncbi.gene-ontologies.db.txt | wc -l
#7963
cat go_term_database_output.txt | wc -l
#6932

cat Apis.mellifera.ncbi.gene-ontologies.db.txt | sed "s/ /\t/g" > Apis.mellifera.GOterms.tsv


```

## 4. TopGO run TopGo enrichment analysis (change node size to gauge strictness)
```
GOTERMS="$PWD/Apis.mellifera.GOterms.tsv"
mkdir -p genes_of_interest
cat genes.brains.tsv | cut -f1 > genes_of_interest/brains.txt
cat genes.gonads.tsv | cut -f1 > genes_of_interest/gonads.txt
GENELIST1="./genes_of_interest/brains.txt"
GENELIST2="./genes_of_interest/gonads.txt"

grep -w -f $GENELIST1 $GOTERMS | sed "s/ /\t/g" > ./genes_of_interest/brains.go.tsv
grep -w -f $GENELIST2 $GOTERMS | sed "s/ /\t/g" > ./genes_of_interest/gonads.go.tsv


## TopGo
ANALYSIS="brains"
INDIR="$PWD"
GENEUNIVERSEgoTERMS="$GOTERMS"
TESTGENESET="brains.txt"
TESTGENESETDIR="genes_of_interest"
OUTDIR="results.topGO.$ANALYSIS"
mkdir -p $OUTDIR

cat ~/scripts/go_enrichment_analysis_fisher.R | sed "s,XXXXXX,$INDIR,g ; s,YYYYYY,$GENEUNIVERSEgoTERMS,g ; s,ZZZZZZ,$TESTGENESET,g ; s,FFFFFF,$TESTGENESETDIR,g ; s,OOOOOO,$OUTDIR,g" > go_enrichment_analysis_fisher.$ANALYSIS.R

Rscript go_enrichment_analysis_fisher.$ANALYSIS.R

## results
ll $OUTDIR/go_term_output_50_fisher/BP_top50.tsv
ll $OUTDIR/go_term_output_50_fisher/CC_top50.tsv
ll $OUTDIR/go_term_output_50_fisher/MF_top50.tsv
```

## 5. visualize results with [GOfigure](https://gitlab.com/evogenlab/GO-Figure)
```
######
gofigure.py --help
#gofigure.py -i input_file.tsv -o out_directory

### TopGO outputs
ANALYSIS="brains"
INDIR="$PWD"
OUTDIR="results.topGO.$ANALYSIS"


$OUTDIR/go_term_output_50_fisher/BP_top50.tsv
$OUTDIR/go_term_output_50_fisher/CC_top50.tsv
$OUTDIR/go_term_output_50_fisher/MF_top50.tsv

conda activate R1 (contains GOfigure)


GO_CATEGORY1="BP"
GO_CATEGORY2="CC"
GO_CATEGORY2="MF"

#plot single categoy
#cat $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_top50.tsv | cut -f1-7 > $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_top50.1.tsv
#
#N=$(cat $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_top50.tsv | tail -n+2 | wc -l)
#echo $N
#paste <(seq 1 $N) <(cat $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_top50.tsv | tail -n+2) | tee $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_top50.formatted.tsv | cut -f1-7 > $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_top50.formatted.1.tsv
#gofigure.py -j topgo -i $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_top50.formatted.1.tsv -o $OUTDIR


#let gofigure plot all
N=$(cat $OUTDIR/go_term_output_50_fisher/*_top50.tsv | tail -n+2 | wc -l)
echo $N
#K=$(echo $N-1 | bc -l)

paste <(seq 1 $N) <(cat $OUTDIR/go_term_output_50_fisher/*_top50.tsv | tail -n+2) | tee $OUTDIR/go_term_output_50_fisher/all_top50.formatted.tsv | cut -f1-7 > $OUTDIR/go_term_output_50_fisher/all_top50.formatted.1.tsv
gofigure.py -j topgo -i $OUTDIR/go_term_output_50_fisher/all_top50.formatted.1.tsv -o $OUTDIR

gofigure.py -j topgo -i $OUTDIR/go_term_output_50_fisher/all_raw.formatted.1.tsv -o $OUTDIR

