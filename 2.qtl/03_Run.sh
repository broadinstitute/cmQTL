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
#BSUB -o STDfiles/QQ_${pheno}_${cell}.out
#BSUB -e STDfiles/QQ_${pheno}_${cell}.err
#BSUB -J QQ
#BSUB -M 12000
#BSUB -q rerunnable
#BSUB -R 'rusage[mem=12000]'
#BSUB -R 'select[hname!=cn001]'
#BSUB -R 'select[hname!=cn002]'
#BSUB -R 'select[hname!=cn003]'
#BSUB -R 'select[hname!=cn004]'
#BSUB -R 'select[hname!=cn005]'

/apps/lib-osver/R/3.4.0/bin/Rscript 03_GetLambda_QQ.R $pheno $cell
 
EOF
}

cat Pheno_Cell.list | awk '{if ($1 < 4) print $0}' | while read pheno cell 
        do
        write_Rscript $pheno $cell
        done
