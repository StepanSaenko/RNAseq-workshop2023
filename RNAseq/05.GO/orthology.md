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
