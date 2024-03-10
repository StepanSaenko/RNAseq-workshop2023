# Mapping RNAseq reads with STAR

## 1. define reference and annotation
```
#Apis mellifera reference
DIR="/scratch/rnaseq2023"
REFDIR="$DIR/Apis.mellifera.reference.genome"
INPUT_FASTA="GCF_003254395.2_Amel_HAv3.1_genomic.fna"
INPUT_GTF="GCF_003254395.2_Amel_HAv3.1_genomic.gtf"
```

## 2. Mapping
```
READDIR="$DIR/raw.data.apis"
CPUs=80
OVERHANG="100"
## readlength -1
ls -1 $READDIR/*.lite.1_pass_1.fastq | rev | cut -d "/" -f1 | cut -d "." -f4- | rev | parallel -k -j 1 "echo {}"
ls -1 $READDIR/*.lite.1_pass_1.fastq | rev | cut -d "/" -f1 | cut -d "." -f4- | rev | parallel -k -j 1 "echo {}; STAR --genomeDir $REFDIR --runThreadN $CPUs --readFilesIn $READDIR/{}.lite.1_pass_1.fastq $READDIR/{}.lite.1_pass_2.fastq --outSAMtype BAM Unsorted --outBAMsortingThreadN 15 --alignEndsType EndToEnd --quantMode GeneCounts --outFileNamePrefix $PWD/results/{}.; samtools sort --output-fmt BAM --threads 15 -l 5 -T {}.tmp --write-index --verbosity 1 -o $PWD/results/{}.Aligned.out.CoordSrted.bam $PWD/results/{}.Aligned.out.bam"
```

## 3. Qualimap
If you run QualiMap in parallel for many samples, make sure to create a different tmp-folder for each sample; e.g., ./tmp/SRR3091420_1_chr6Aligned
#Make sure to give to the output folder the name corresponding to a running sample; e.g., ./QC/SRR3091420_1_chr6; otherwise output files will be overwritten
#QualiMap sorts BAM files by read names. To speed up this part of the program execution, you can use samtools to sort the BAM files in parallel and using multiple CPUs and then to give to QualiMap a BAM file sorted by read names and provide an option –sorted
```
mkdir -p results/qualimap
ls -lh results
ls -1 $READDIR/*.lite.1_pass_1.fastq | rev | cut -d "/" -f1 | cut -d "." -f4- | rev | parallel -k -j 1 "echo {}; qualimap rnaseq -outdir ./results/qualimap/{} -a proportional -bam ./results/{}.Aligned.out.CoordSrted.bam -p strand-specific-reverse -gtf $REFDIR/$INPUT_GTF --java-mem-size=20G"
```






