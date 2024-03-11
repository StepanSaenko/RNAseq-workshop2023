## GO term enrichment

### Steps
1. obtain a set of gene IDs you want to analyze
2. determine 1:1 orthologs between your species and a species available in ENSEMBL/BioMart (eg Apis mellifera, Bombus terrestris, Bombus impatiens, Nasionia vitripennis, Drosophila melanogaster ...). You can do this with the list of genes of interest or all genes of an annotation in that species. We can use broccoli for this task
3. fetch GO terms for the ENSEMBL/BioMART species and it's Drosophila orthologs
4. create genelist subsets
5. enrichment (TopGO)
6. visualitation (GOfigure)





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

