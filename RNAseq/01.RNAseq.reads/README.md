## RNAseq read processing

### Slater 2022
- Haploid and Sexual Selection Shape the Rate of Evolution of Genes across the Honey Bee (Apis mellifera L.) Genome
- PRJNA689223
- [DOI](https://doi.org/10.1093/gbe/evac063)
- [PAPERLINK](https://doi.org/10.1093/gbe/evac063](https://academic.oup.com/gbe/article/14/6/evac063/6584681)
- [NCBI](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=4&WebEnv=MCID_65ece0921f4e9f73e3a78866&o=acc_s%3Aa)

## 1. prepare the data download
```
sudo mkdir -p /scratch/rnaseq2023/raw.data.apis
cd /scratch/rnaseq2023/raw.data.apis

ffq 10.1093/gbe/evac063
#doesnt work
```

##### List of files/samples
```
#List of files/samples
cat PRJNA689223.txt
SRR13361446
SRR13361447
SRR13361448
SRR13361449
SRR13361450
SRR13361451
SRR13361452
SRR13361453
SRR13361454
SRR13361455
SRR13361456
SRR13361457
SRR13361458
SRR13361459
SRR13361460
SRR13361461
SRR13361462
```
## 2. Download the data
```
cat PRJNA689223.txt | parallel -k -j2 "ffq --ncbi {} | jq -r '.[] | .url' | xargs curl -O"
```
## 3. transform into fastq
```
cat PRJNA689223.txt | parallel -k -j2 "echo {}; fastq-dump --outdir . --skip-technical --readids --dumpbase --split-3 --read-filter pass --clip {}.lite.1"
```

## 4. install STAR and kallisto
```
cd ~/progz
git clone https://github.com/pachterlab/kallisto
cd kallisto
mkdir build
cd build
cmake ..
make -j5
ln -s $PWD/src/kallisto ~/bin

mkdir -p ~/bin
mkdir -p ~/progz; cd ~/progz
wget https://github.com/alexdobin/STAR/archive/2.7.11a.tar.gz
tar -xzf 2.7.11a.tar.gz
cd STAR-2.7.11a/bin/Linux_x86_64_static/
ln -s $PWD/STAR ~/bin/
ln -s $PWD/STAR ~/bin/STAR_exec
ln -s $PWD/STARlong ~/bin/
```

## 5. Reference Genome and Annotation
```
#Apis mellifera reference
REFERENCEFOLDER="/scratch/rnaseq2023/Apis.mellifera.reference.genome"
REF="GCF_003254395.2_Amel_HAv3.1_genomic.fna"
GTF="GCF_003254395.2_Amel_HAv3.1_genomic.gtf"
GFF="GCF_003254395.2_Amel_HAv3.1_genomic.gff"
```








