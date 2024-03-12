## Background to steps taken in preparation for practical:  
1. Alignment of RNA-seq data  
Each sample were aligned using the short-read RNA-seq aligner STAR.  

The steps for analysis included:
- Indexing reference genome assembly, which was performed using STAR.   

```
STAR --runThreadN 20 \
     --runMode genomeGenerate \
     --genomeDir ./database/ \
     --genomeFastaFiles "$input_fasta" \
     --sjdbGTFfile "$input_gtf" \
     --genomeSAindexNbases 12 \
     --sjdbOverhang "$overhang"
```

- Alignment of samples against the STAR-indexed reference genome assembly.  

```
STAR  \
--genomeDir ./database/ \
--runThreadN 20 \
--readFilesCommand gunzip -c \
--readFilesIn "$input_pair_1" "$input_pair_2" \
--outSAMtype BAM Unsorted \
--outBAMsortingThreadN 2 \
--quantMode GeneCounts \
--outFileNamePrefix ./results/"$output".  
```
