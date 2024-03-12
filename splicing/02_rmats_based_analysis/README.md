## Splicing event identification and quantification using rMATs  
For the present tutorial, rMATs was used for identifying and quantifying 
different splicing events between honeybee sexes.  

Due to time-constraints, we will provide and analysis the output files 
from rMATs, which was run seperately for brains and gonads, respectively. 

The following commands were used:  
python src/rmats_turbo_v4_1_2/rmats.py \
--b1 bam_list_male_gonads.txt  \
--b2 bam_list_queen_gonads.txt  \
--readLength 150  \
--variable-read-length  \
--gtf /scratch/rnaseq2023/Apis.mellifera.reference.genome/GCF_003254395.2_Amel_HAv3.1_genomic.gtf  \
--nthread 20  \
--od results_gonad_comp/  \
--tmp tmp_gonads  
