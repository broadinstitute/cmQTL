PROJECT_NAME=2018_06_05_GWAS_in_a_dish

BATCH_ID=PILOT_1

BUCKET=imaging-platform

MAXPROCS=3 # m4.xlarge has 4 cores; keep 1 free

mkdir -p ~/efs/${PROJECT_NAME}/workspace/

cd ~/efs/${PROJECT_NAME}/workspace/

mkdir -p log/${BATCH_ID}

PIPELINE_SET=cellpainting_ipsc_20x_phenix_with_bf_bin1

PLATES=$(readlink -f ~/efs/${PROJECT_NAME}/workspace/scratch/${BATCH_ID}/plates_to_process.txt)
