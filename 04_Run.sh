#!/bin/bash
module load R/3.4.0

function write_Rscript()
{
        local pheno=$1
        local cell=$2
	module load R/3.4.0

cat << EOF | bsub

export LSB_JOB_REPORT_MAIL=N
#BSUB -L /bin/bash
#BSUB -o STDfiles/ForManh_${pheno}_${cell}.out
#BSUB -e STDfiles/ForManh_${pheno}_${cell}.err
#BSUB -J DownSamp
#BSUB -M 8000
#BSUB -q rerunnable
#BSUB -R 'rusage[mem=6000]'
#BSUB -R 'select[hname!=cn001]'
#BSUB -R 'select[hname!=cn002]'
#BSUB -R 'select[hname!=cn003]'
#BSUB -R 'select[hname!=cn004]'
#BSUB -R 'select[hname!=cn005]'

/apps/lib-osver/R/3.4.0/bin/Rscript 04_GetdataForManhattanPlot.R $pheno $cell

#gzip -f /data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/${cell}/${cell}_pheno${pheno}.DownSampled.tmp
 
EOF
}

cat Pheno_Cell.list | while read pheno cell 
        do
        write_Rscript $pheno $cell
        done

