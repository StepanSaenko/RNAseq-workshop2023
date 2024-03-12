## Readcounts with Kallisto

## 1. Reference
```
#Apis mellifera reference
REFERENCEFOLDER="/scratch/rnaseq2023/Apis.mellifera.reference.genome"
REF="GCF_003254395.2_Amel_HAv3.1_genomic.fna"
GTF="GCF_003254395.2_Amel_HAv3.1_genomic.gtf"
CDS="cds_from_genomic.fna"
TRANSCRIPTS="rna.fna"
GFF="GCF_003254395.2_Amel_HAv3.1_genomic.gff"

ll $REFERENCEFOLDER/$GTF
/scratch/rnaseq2023/Apis.mellifera.reference.genome/GCF_003254395.2_Amel_HAv3.1_genomic.gtf
```

## 2. make genomic index
```
cd /scratch/rnaseq2023
#Generate an index
kallisto index -i $REFERENCEFOLDER/amel.transcripts $REFERENCEFOLDER/$TRANSCRIPTS


```

## 3. Kallisto counts for all samples
```
cd /scratch/rnaseq2023/USER
INPUTFOLDER="/scratch/rnaseq2023/raw.data.apis"
cp $INPUTFOLDER/PRJNA689223.txt ./
OUTFOLDER="../kallisto.output"
#mkdir -p $OUTFOLDER
CPUs=5

# find your sample
ls *.fastq
SAMPLE="XXXXXXXXXXX"

# Kallisto counts for all samples
## dry-run
cat PRJNA689223.txt | grep "$SAMPLE" | parallel -k -j1 "echo ${OUTFOLDER}/{} && mkdir -p ${OUTFOLDER}/{}"

## kallisto without the --genomebam option
cat PRJNA689223.txt | grep "$SAMPLE" | parallel -k -j1 "echo {}; kallisto quant -i ${REFERENCEFOLDER}/amel.transcripts -l 130 -s 20 -t $CPUs -o ${OUTFOLDER}/{} --gtf ${REFERENCEFOLDER}/${GTF} --plaintext ${INPUTFOLDER}/{}.lite.1_pass_1.fastq ${INPUTFOLDER}/{}.lite.1_pass_2.fastq"

```

## 4. alternatively, salmon can be used
```
OUTFOLDER="./2024-03-06.salmon.output"
mkdir -p $OUTFOLDER
CPUs=50
salmon index --index ${REFERENCEFOLDER}/amel.transcripts.salmon --transcripts $REFERENCEFOLDER/$TRANSCRIPTS
salmon quant --help-reads
CPUs=30

## transcripts to genes table example
transcript_id   gene_id
ENST00000456328.2   ENSG00000223972.5
ENST00000461467.1   ENSG00000237613.2
cat $REFERENCEFOLDER/$GTF | awk '$3 == "transcript" ' | cut -f9 | cut -d";" -f1,2 | tr -d ";" | sed "s/ transcript_id /\t/g" | sed "s/gene_id //g" | tr -d "\"" > $REFERENCEFOLDER/$GTF.gene.transcript.map.txt
cat $REFERENCEFOLDER/$GTF.gene.transcript.map.txt | tabtk cut -r -f2,1 > transcript.gene.map.txt

#run salmon
cat PRJNA689223.txt | parallel -k -j2 "echo {}; time salmon quant --index ./amel.transcripts.salmon -l A -p ${CPUs} --validateMappings --mates1 ${INPUTFOLDER}/{}.lite.1_pass_1.fastq --mates2 ${INPUTFOLDER}/{}.lite.1_pass_2.fastq --output ${OUTFOLDER}/{}/ --geneMap transcript.gene.map.txt"

#merge
mkdir -p 2024-03-06.salmon.merge
salmon quantmerge --output 2024-03-06.salmon.merge/PRJNA689223.quantmerge.txt --quants ${OUTFOLDER}/*
salmon quantmerge --output 2024-03-06.salmon.merge/PRJNA689223.quantmerge.genes.txt --genes --quants ${OUTFOLDER}/*

```







