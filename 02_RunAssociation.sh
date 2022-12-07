#!/bin/bash

module load plink/1.90b3
module load gcta/1.91.1

export LSB_JOB_REPORT_MAIL=N

function write_Rscript()
{
	local outDir=/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov
	local chr=$1
	local pheno=$2
	local cell=$3

cat << EOF | bsub

export LSB_JOB_REPORT_MAIL=N
#BSUB -L /bin/bash
#BSUB -o STDfiles/fastgwa_${chr}_${pheno}_${cell}.out
#BSUB -e STDfiles/fastgwa_${chr}_${pheno}_${cell}.err
#BSUB -J fastgwa
#BSUB -M 8000
#BSUB -q rerunnable
#BSUB -R 'rusage[mem=8000]'
#BSUB -R 'select[hname!=cn001]'
#BSUB -R 'select[hname!=cn002]'
#BSUB -R 'select[hname!=cn003]'
#BSUB -R 'select[hname!=cn004]'
#BSUB -R 'select[hname!=cn005]'

/data/srlab/sasgari/Software/gcta_1.93.2beta/gcta64 --bfile  ../SecondRound_JatinFeatures/PlinkFiles/V2F_chr${chr}.nodup.qc \
--pheno ${outDir}/${cell}/${cell}.pheno \
--mpheno ${pheno} \
--covar ${outDir}/${cell}/${cell}.categorical.cov \
--qcovar ${outDir}/${cell}/${cell}.quantitative.cov \
--remove ../SecondRound_JatinFeatures/PlinkFiles/ExcludeFromAssociation.ID \
--maf 0.05 \
--fastGWA-lr \
--threads 10 \
--out ${outDir}/${cell}/${cell}_chr${chr}_pheno${pheno}.LR

gzip -f ${outDir}/${cell}/${cell}_chr${chr}_pheno${pheno}.LR.fastGWA

EOF
}


cat Chr_Pheno_Cell.list | while read chr pheno cell 
	do
	write_Rscript $chr $pheno $cell
	done
