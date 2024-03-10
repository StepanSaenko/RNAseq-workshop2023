## GO term enrichment

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


### compleasm
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

