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
conda create -n R1
conda activate R1
mamba install -c conda-forge r-base=4.3
mamba install -c conda-forge gcc
mamba install -c conda-forge curl
mamba install -c conda-forge r-curl
mamba install -c conda-forge r-tidyverse
mamba install -c conda-forge r-png
mamba install r-png r-rcurl
mamba install -c conda-forge -c bioconda go-figure

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


cd ~/progz
git clone https://gitlab.com/evogenlab/GO-Figure.git
cd GO-Figure
ln -s $PWD/gofigure.py ~/bin
ln -s $PWD/gofigure.py ~/bin/gofigure

cd Data
wget -P ~/progz/GO-Figure/data/ https://purl.obolibrary.org/obo/go/go-basic.obo
wget -P ~/progz/GO-Figure/data/ https://purl.obolibrary.org/obo/go.obo
wget -P ~/progz/GO-Figure/data/ https://purl.obolibrary.org/obo/go.owl
wget -P ~/progz/GO-Figure/data/ https://purl.obolibrary.org/obo/go/extensions/go-plus.owl
wget -P ~/progz/GO-Figure/data/ ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gpa.gz
wget -P ~/progz/GO-Figure/data/ ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gpi.gz
wget -P ~/progz/GO-Figure/data/ ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_gcrp.gaf.gz
wget -P ~/progz/GO-Figure/data/ ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_gcrp.gpa.gz
wget -P ~/progz/GO-Figure/data/ ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_gcrp.gpi.gz
wget -P ~/progz/GO-Figure/data/ ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz

python3 ~/progz/GO-Figure/scripts/relations.py ~/progz/GO-Figure/data/go.obo > ~/progz/GO-Figure/data/relations.tsv
python3 ~/progz/GO-Figure/scripts/ics.py ~/progz/GO-Figure/data/relations.tsv ~/progz/GO-Figure/data/goa_uniprot_all.gaf.gz ~/progz/GO-Figure/data/go.obo > ~/progz/GO-Figure/data/ic.tsv
```
