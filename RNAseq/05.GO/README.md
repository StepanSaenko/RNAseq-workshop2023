## GO term enrichment

### Steps
1. obtain a set of gene IDs you want to analyze
2. determine 1:1 orthologs between your species and a species available in ENSEMBL/BioMart (eg Apis mellifera, Bombus terrestris, Bombus impatiens, Nasionia vitripennis, Drosophila melanogaster ...). You can do this with the list of genes of interest or all genes of an annotation in that species. We can use broccoli for this task
3. fetch GO terms for the ENSEMBL/BioMART species and it's Drosophila orthologs
4. create genelist subsets
5. enrichment (TopGO)
6. visualitation (GOfigure)



### install broccoli for orthology inference
```
conda create -n env-broccoli python=3.6 ete3 fasttree diamond
cd ~/progz/
git clone https://github.com/rderelle/Broccoli.git
conda activate env-broccoli
```

### install R and packages
```
# R for using biomaRt

#if Linux doesnt have libcurl installed
#sudo apt install libcurl4-openssl-dev libpng-dev

conda create -n R1
conda activate R1
mamba install -c conda-forge r-base=4.3
mamba install -c conda-forge gcc
mamba install -c conda-forge curl
mamba install -c conda-forge r-curl
mamba install -c conda-forge r-tidyverse
mamba install -c conda-forge r-png
mamba install r-png r-rcurl

#start R
R

#install bioconductor
install.packages('tidyverse', dependencies=TRUE)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

#check which packages are avail. if interested
BiocManager::available()

#install biomart
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("KEGGREST")
BiocManager::install("httr")
BiocManager::install("XML")
BiocManager::install("AnnotationDbi")
BiocManager::install("BiocFileCache")
BiocManager::install("xml2")
BiocManager::install("biomaRt")
BiocManager::install("topGO")
BiocManager::install("lintr")
BiocManager::install("lattice")
BiocManager::install("Rgraphviz")
```

### install GOfigure
```
conda activate R1
mamba install -c conda-forge -c bioconda go-figure

cd /scratch/progz
git clone https://gitlab.com/evogenlab/GO-Figure.git
cd GO-Figure
ln -s $PWD/gofigure.py /usr/local/bin
ln -s $PWD/gofigure.py /usr/local/bin/gofigure

cd Data
wget -P /scratch/progz/GO-Figure/data/ https://purl.obolibrary.org/obo/go/go-basic.obo
wget -P /scratch/progz/GO-Figure/data/ https://purl.obolibrary.org/obo/go.obo
wget -P /scratch/progz/GO-Figure/data/ https://purl.obolibrary.org/obo/go.owl
wget -P /scratch/progz/GO-Figure/data/ https://purl.obolibrary.org/obo/go/extensions/go-plus.owl
wget -P /scratch/progz/GO-Figure/data/ ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gpa.gz
wget -P /scratch/progz/GO-Figure/data/ ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gpi.gz
wget -P /scratch/progz/GO-Figure/data/ ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_gcrp.gaf.gz
wget -P /scratch/progz/GO-Figure/data/ ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_gcrp.gpa.gz
wget -P /scratch/progz/GO-Figure/data/ ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_gcrp.gpi.gz
wget -P /scratch/progz/GO-Figure/data/ ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz

python3 /scratch/progz/GO-Figure/scripts/relations.py /scratch/progz/GO-Figure/data/go.obo > /scratch/progz/GO-Figure/data/relations.tsv
python3 /scratch/progz/GO-Figure/scripts/ics.py /scratch/progz/GO-Figure/data/relations.tsv /scratch/progz/GO-Figure/data/goa_uniprot_all.gaf.gz /scratch/progz/GO-Figure/data/go.obo > /scratch/progz/GO-Figure/data/ic.tsv
```
### install compleasm
```
### install tool
cd ~/progz
wget https://github.com/huangnengCSU/compleasm/releases/download/v0.2.5/compleasm-0.2.5_x64-linux.tar.bz2
tar -jxvf compleasm-0.2.5_x64-linux.tar.bz2
pip3 install pandas
#ln -s $PWD/hmmsearch ~/bin/
ln -s $PWD/miniprot ~/bin/
ln -s $PWD/compleasm.py ~/bin/compleasm
ln -s $PWD/compleasm.py ~/bin/

### get lineages (centrally, do not duplicate
compleasm list --remote | grep -P 'metazoa|hymenoptera|insecta|arthropoda|diptera|coleoptera|lepidoptera|endopterygota'
arthropoda_odb10
diptera_odb10
endopterygota_odb10
hymenoptera_odb10
insecta_odb10
lepidoptera_odb10
metazoa_odb10

mkdir -p /share/pool/ek/database/compleasm_lineages
chown -R estolle:CGI /share/pool/ek/database/compleasm_lineages
cd /share/pool/ek/database/compleasm_lineages
compleasm download metazoa,hymenoptera,insecta,arthropoda,diptera,lepidoptera,endopterygota -L /share/pool/ek/database/compleasm_lineages

### lineages are here -L /share/pool/ek/database/compleasm_lineages

### run example
REF="xxxxx.fa"
OUTDIR="/xxxx/xxx"
CPUs=50
compleasm run -a $REF -o $OUTDIR -t $CPUs -m busco -L /share/pool/ek/database/compleasm_lineages -l hymenoptera
```


## 1. compleasm to generate gene list (this example)
```
REF="/scratch/ek/genomes/Bombus.terrestris/Bombus_terrestris-GCA_910591885.1-iyBomTerr1.2-RefSeq/GCF_910591885.1_iyBomTerr1.2_genomic.fna"
GFF="/scratch/ek/genomes/Bombus.terrestris/Bombus_terrestris-GCA_910591885.1-iyBomTerr1.2-RefSeq/GCF_910591885.1_iyBomTerr1.2_genomic.gff"

OUTDIR="$REF.compleasm"
CPUs=50
compleasm run -a $REF -o $OUTDIR -t $CPUs -m busco -L /scratch/database/compleasm_lineages/ -l hymenoptera
rm -rf $OUTDIR/hymenoptera_odb10/hmmer_output
pigz $OUTDIR/hymenoptera_odb10/miniprot_output.gff
pigz $OUTDIR/hymenoptera_odb10/translated_protein.fasta
pigz $OUTDIR/hymenoptera_odb10/gene_marker.fasta

COMPLEASMPROTEIN="$OUTDIR/hymenoptera_odb10/gene_marker.fasta"
COMPLEASMTABLE="$OUTDIR/hymenoptera_odb10/full_table_busco_format.tsv"
APISPROTEIN="/scratch/ek/genomes/Apis.mellifera/Apis.mellifera.HAv3.1.GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_protein.faa"
```

## 2. orthology
```
# use env
conda activate /scratch/progz/conda_envs/orthology
python /scratch/progz/conda_envs/orthology/bin/broccoli.py -h

#establish orthology relationship vs Apis
mkdir -p orthologs.apis.broccoli
cd orthologs.apis.broccoli

#cp input files and rename
cp $APISPROTEIN $COMPLEASMPROTEIN.gz .
pigz -d *.gz
mv gene_marker.fasta BombusTerrestris.fasta
mv GCF_003254395.2_Amel_HAv3.1_protein.faa ApisMellifera.fasta

# orthology # 9min
time python /scratch/progz/conda_envs/orthology/bin/broccoli.py -dir ./ -threads 20 -not_same_sp

# check output, merge/format
QUERYSPECIES="BombusTerrestris"
REFSPECIES="ApisMellifera"
cat dir_step4/orthologous_pairs.txt | head
341at7399_6     XP_026302264.1
1329at7399_9    XP_026302264.1
11770at7399     XP_026295193.1
13022at7399     XP_026297287.1

cat dir_step4/orthologous_pairs.txt | cut -f1 | sort | wc -l
14570
cat dir_step4/orthologous_pairs.txt | cut -f1 | sort | uniq | wc -l
8829
cat dir_step4/orthologous_pairs.txt | cut -f1 | sort | uniq -c | sort -k1,1nr | head 
     51 352at7399
     50 0at7399_5
     42 3174at7399_2

cat dir_step4/orthologous_pairs.txt | cut -f1 | grep "^XP_" | wc -l
4689
cat dir_step4/orthologous_pairs.txt | grep "^XP_" | head
XP_001120319.3  21434at7399
XP_006561939.1  19701at7399
XP_003249953.1  19701at7399


cat dir_step4/orthologous_pairs.txt | grep "^XP_" > $REFSPECIES.$QUERYSPECIES.1v1_orthologs.tsv
cat dir_step4/orthologous_pairs.txt | grep -v "^XP_" > $QUERYSPECIES.$REFSPECIES.1v1_orthologs.tsv

#take input orthology list and place each species on only 1 side
sort -k1,1 <(cat <(cat dir_step4/orthologous_pairs.txt | grep -v "^XP_") <(cat dir_step4/orthologous_pairs.txt | grep "^XP_" | tabtk cut -r -f2,1)) > $QUERYSPECIES.$REFSPECIES.1v1_orthologs.all.tsv

cat $QUERYSPECIES.$REFSPECIES.1v1_orthologs.all.tsv | wc -l
# 14570
#BombusTerrestris.ApisMellifera.1v1_orthologs.all.tsv

cat $QUERYSPECIES.$REFSPECIES.1v1_orthologs.all.tsv | tabtk cut -r -f1,2 | sort -k1,1 > $REFSPECIES.$QUERYSPECIES.1v1_orthologs.all.tsv
#ApisMellifera.BombusTerrestris.1v1_orthologs.all.tsv
```

## 3. create GO Term lists via BioMaRT
```
cd /scratch/ek/comparative/2023phylogenyB
cd orthologs.apis.broccoli
QUERYSPECIES="BombusTerrestris"
REFSPECIES="ApisMellifera"

# use general Rscript and modify accodrding to the reference species (available: BTER, BIMP, AMEL, NVIT)
INPUTFILE="$QUERYSPECIES.$REFSPECIES.1v1_orthologs.all.tsv"
cat ~/scripts/go_term_extraction.R | sed -r "s#XXXXXX#$INPUTFILE#g" | sed "s/##AMEL //g" > go_term_extraction.Amel.R
chmod 755 go_term_extraction.Amel.R

#Run R script
Rscript go_term_extraction.Amel.R

#output 
ll $INPUTFILE.GO.Amel-Dmel.raw.txt

#optional run again using Bter,Bimp,Nvit
# need to run orthology for this
cat ~/scripts/go_term_extraction.R | sed -r "s#XXXXXX#$INPUTFILE#g" | sed "s/##BTER //g" > go_term_extraction.Bter.R
chmod 755 go_term_extraction.Bter.R
Rscript go_term_extraction.Bter.R
ll $INPUTFILE.GO.Bter-Dmel.raw.txt
cat ~/scripts/go_term_extraction.R | sed -r "s#XXXXXX#$INPUTFILE#g" | sed "s/##NVIT //g" > go_term_extraction.Nvit.R
chmod 755 go_term_extraction.Nvit.R
Rscript go_term_extraction.Nvit.R
ll $INPUTFILE.GO.Nvit-Dmel.raw.txt
cat ~/scripts/go_term_extraction.R | sed -r "s#XXXXXX#$INPUTFILE#g" | sed "s/##BIMP //g" > go_term_extraction.Bimp.R
chmod 755 go_term_extraction.Bimp.R
Rscript go_term_extraction.Bimp.R
ll $INPUTFILE.GO.Bimp-Dmel.raw.txt


# merge GO term lists (from above and/or additional sources)
cat *.GO.*.raw.txt | grep -v "Uncharacterized" > go.raw.tsv
Rscript ~/scripts/go_term_conversion.R go.raw.tsv go.converted.tsv
cat go.converted.tsv | wc -l
#3744

cat go.converted.tsv > Busco.Apis.GOterms.tsv
```

## 4. collect files
```
cd /scratch/ek/comparative/2023phylogenyB

mkdir -p ./08.hyphy.relax.GO/genes_of_interest

GOTERMS="$PWD/orthologs.apis.broccoli/Busco.Apis.GOterms.tsv"
GENES1="./07.hyphy.relax.completed.genenames.nochange.txt"
GENES2="./07.hyphy.relax.completed.genenames.intensification.txt"
GENES3="./07.hyphy.relax.completed.genenames.relaxation.txt"

cp $GOTERMS ./08.hyphy.relax.GO
sed -i "s/ /\t/g" ./08.hyphy.relax.GO/Busco.Apis.GOterms.tsv
cp $GENES1 $GENES2 $GENES3 ./08.hyphy.relax.GO/genes_of_interest
grep -w -f $GENES1 $GOTERMS | sed "s/ /\t/g" > ./08.hyphy.relax.GO/genes_of_interest/nohange.go.tsv
grep -w -f $GENES2 $GOTERMS | sed "s/ /\t/g" > ./08.hyphy.relax.GO/genes_of_interest/intesified.go.tsv
grep -w -f $GENES3 $GOTERMS | sed "s/ /\t/g" > ./08.hyphy.relax.GO/genes_of_interest/relaxed.go.tsv
```

## 5. run TopGo enrichment analysis (change node size to gauge strictness)

### Activate env
```
conda activate R1
cd /scratch/ek/comparative/2023phylogenyB/08.hyphy.relax.GO
```

### 1. intensification
```
ANALYSIS="Meliponini-VS-Apidae.intensification"
INDIR="/scratch/ek/comparative/2023phylogenyB/08.hyphy.relax.GO"
GENEUNIVERSEgoTERMS="Busco.Apis.GOterms.tsv"
TESTGENESET="07.hyphy.relax.completed.genenames.intensification.txt"
TESTGENESETDIR="genes_of_interest"
OUTDIR="results.$ANALYSIS"
mkdir -p $OUTDIR
cat ~/scripts/go_enrichment_analysis_fisher.R | sed "s,XXXXXX,$INDIR,g ; s,YYYYYY,$GENEUNIVERSEgoTERMS,g ; s,ZZZZZZ,$TESTGENESET,g ; s,FFFFFF,$TESTGENESETDIR,g ; s,OOOOOO,$OUTDIR,g" > go_enrichment_analysis_fisher.$ANALYSIS.R
Rscript go_enrichment_analysis_fisher.$ANALYSIS.R
```

### 2. no change
```
ANALYSIS="Meliponini-VS-Apidae.nochange"
INDIR="/scratch/ek/comparative/2023phylogenyB/08.hyphy.relax.GO"
GENEUNIVERSEgoTERMS="Busco.Apis.GOterms.tsv"
TESTGENESET="07.hyphy.relax.completed.genenames.nochange.txt"
TESTGENESETDIR="genes_of_interest"
OUTDIR="results.$ANALYSIS"
mkdir -p $OUTDIR
cat ~/scripts/go_enrichment_analysis_fisher.R | sed "s,XXXXXX,$INDIR,g ; s,YYYYYY,$GENEUNIVERSEgoTERMS,g ; s,ZZZZZZ,$TESTGENESET,g ; s,FFFFFF,$TESTGENESETDIR,g ; s,OOOOOO,$OUTDIR,g" > go_enrichment_analysis_fisher.$ANALYSIS.R
Rscript go_enrichment_analysis_fisher.$ANALYSIS.R
```

### 3. relaxation
```
ANALYSIS="Meliponini-VS-Apidae.relaxation"
INDIR="/scratch/ek/comparative/2023phylogenyB/08.hyphy.relax.GO"
GENEUNIVERSEgoTERMS="Busco.Apis.GOterms.tsv"
TESTGENESET="07.hyphy.relax.completed.genenames.relaxation.txt"
TESTGENESETDIR="genes_of_interest"
OUTDIR="results.$ANALYSIS"
mkdir -p $OUTDIR
cat ~/scripts/go_enrichment_analysis_fisher.R | sed "s,XXXXXX,$INDIR,g ; s,YYYYYY,$GENEUNIVERSEgoTERMS,g ; s,ZZZZZZ,$TESTGENESET,g ; s,FFFFFF,$TESTGENESETDIR,g ; s,OOOOOO,$OUTDIR,g" > go_enrichment_analysis_fisher.$ANALYSIS.R
Rscript go_enrichment_analysis_fisher.$ANALYSIS.R
```

### results
```
$OUTFOLDER/go_term_output_50_fisher/BP_raw.tsv
$OUTFOLDER/go_term_output_50_fisher/CC_raw.tsv
$OUTFOLDER/go_term_output_50_fisher/MF_raw.tsv
```

## 5. visualize results with [GOfigure](https://gitlab.com/evogenlab/GO-Figure)
```
######
gofigure.py --help
gofigure.py -i input_file.tsv -o out_directory

### TopGO outputs
ANALYSIS="Meliponini-VS-Apidae.intensification"
OUTDIR="results.$ANALYSIS"
ANALYSIS="Meliponini-VS-Apidae.nochange"
OUTDIR="results.$ANALYSIS"
ANALYSIS="Meliponini-VS-Apidae.relaxation"
OUTDIR="results.$ANALYSIS"

$OUTFOLDER/go_term_output_50_fisher/BP_raw.tsv
$OUTFOLDER/go_term_output_50_fisher/CC_raw.tsv
$OUTFOLDER/go_term_output_50_fisher/MF_raw.tsv

/scratch/ek/comparative/2023phylogenyB/08.hyphy.relax.GO/results.Meliponini-VS-Apidae.intensification/go_term_output_50_fisher/BP_raw.tsv
/scratch/ek/comparative/2023phylogenyB/08.hyphy.relax.GO/results.Meliponini-VS-Apidae.intensification/go_term_output_50_fisher/CC_raw.tsv
/scratch/ek/comparative/2023phylogenyB/08.hyphy.relax.GO/results.Meliponini-VS-Apidae.intensification/go_term_output_50_fisher/MF_raw.tsv

conda activate R1
cd /scratch/ek/comparative/2023phylogenyB/08.hyphy.relax.GO

ANALYSIS="Meliponini-VS-Apidae.intensification"
OUTDIR="results.$ANALYSIS"
ll $OUTDIR/go_term_output_50_fisher/*_raw.tsv
GO_CATEGORY1="BP"
GO_CATEGORY2="CC"
GO_CATEGORY2="MF"

cat $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.tsv | cut -f1-7 > $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.1.tsv

N=$(cat $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.tsv | tail -n+2 | wc -l)
echo $N
#K=$(echo $N-1 | bc -l)
#550
paste <(seq 1 $N) <(cat $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.tsv | tail -n+2) | tee $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.formatted.tsv | cut -f1-7 > $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.formatted.1.tsv
gofigure.py -j topgo -i $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.formatted.1.tsv -o $OUTDIR

N=$(cat $OUTDIR/go_term_output_50_fisher/*_raw.tsv | tail -n+2 | wc -l)
echo $N
#K=$(echo $N-1 | bc -l)
#550
paste <(seq 1 $N) <(cat $OUTDIR/go_term_output_50_fisher/*_raw.tsv | tail -n+2) | tee $OUTDIR/go_term_output_50_fisher/all_raw.formatted.tsv | cut -f1-7 > $OUTDIR/go_term_output_50_fisher/all_raw.formatted.1.tsv
gofigure.py -j topgo -i $OUTDIR/go_term_output_50_fisher/all_raw.formatted.1.tsv -o $OUTDIR


ANALYSIS="Meliponini-VS-Apidae.nochange"
OUTDIR="results.$ANALYSIS"
ll $OUTDIR/go_term_output_50_fisher/*_raw.tsv
GO_CATEGORY1="BP"
GO_CATEGORY2="CC"
GO_CATEGORY2="MF"

cat $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.tsv | cut -f1-7 > $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.1.tsv

N=$(cat $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.tsv | tail -n+2 | wc -l)
echo $N
#K=$(echo $N-1 | bc -l)
#550
paste <(seq 1 $N) <(cat $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.tsv | tail -n+2) | tee $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.formatted.tsv | cut -f1-7 > $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.formatted.1.tsv
gofigure.py -j topgo -i $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.formatted.1.tsv -o $OUTDIR

N=$(cat $OUTDIR/go_term_output_50_fisher/*_raw.tsv | tail -n+2 | wc -l)
echo $N
#K=$(echo $N-1 | bc -l)
#550
paste <(seq 1 $N) <(cat $OUTDIR/go_term_output_50_fisher/*_raw.tsv | tail -n+2) | tee $OUTDIR/go_term_output_50_fisher/all_raw.formatted.tsv | cut -f1-7 > $OUTDIR/go_term_output_50_fisher/all_raw.formatted.1.tsv
gofigure.py -j topgo -i $OUTDIR/go_term_output_50_fisher/all_raw.formatted.1.tsv -o $OUTDIR


ANALYSIS="Meliponini-VS-Apidae.relaxation"
OUTDIR="results.$ANALYSIS"
ll $OUTDIR/go_term_output_50_fisher/*_raw.tsv
GO_CATEGORY1="BP"
GO_CATEGORY2="CC"
GO_CATEGORY2="MF"

cat $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.tsv | cut -f1-7 > $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.1.tsv

N=$(cat $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.tsv | tail -n+2 | wc -l)
echo $N
#K=$(echo $N-1 | bc -l)
#550
paste <(seq 1 $N) <(cat $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.tsv | tail -n+2) | tee $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.formatted.tsv | cut -f1-7 > $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.formatted.1.tsv
gofigure.py -j topgo -i $OUTDIR/go_term_output_50_fisher/${GO_CATEGORY1}_raw.formatted.1.tsv -o $OUTDIR

N=$(cat $OUTDIR/go_term_output_50_fisher/*_raw.tsv | tail -n+2 | wc -l)
echo $N
#K=$(echo $N-1 | bc -l)
#550
paste <(seq 1 $N) <(cat $OUTDIR/go_term_output_50_fisher/*_raw.tsv | tail -n+2) | tee $OUTDIR/go_term_output_50_fisher/all_raw.formatted.tsv | cut -f1-7 > $OUTDIR/go_term_output_50_fisher/all_raw.formatted.1.tsv
gofigure.py -j topgo -i $OUTDIR/go_term_output_50_fisher/all_raw.formatted.1.tsv -o $OUTDIR

