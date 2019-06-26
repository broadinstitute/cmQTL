#!/bin/bash
#
# cmQTL PILOT processing pipeline
# Gregory Way, 2019 (adapted by Shantanu Singh)
#
# We aggregate two different kinds of profiles from distinct single cell populations
#
# 1) Cells that are formed in colonies
# 2) Isolated Cells
#
# The thought is that biological and genomic properties may be differential in each
# population. This pipeline will generate single cell profiles based on the following
# rules:
#
# 1) Colonies = Number of Neighboring Cells greater than or equal to 4
# 2) Isolated = Number of Neighboring Cells equal to 0
#
# This means that we are removing all cells that have between 1 and 3 neighboring cells
#
# Instructions to generate cell painting profiles for the cmQTL experiments
# Pipeline generated using the profiling handbook:
#     https://cytomining.github.io/profiling-handbook/

############################
# Step 1 - Configure the Environment
############################
# Step 1.1: Setup a virtual machine on AWS
#
# Follow the setup instructions in `generate-pilot-profiles.sh`
# This pipeline will be run in the same VM

# Step 1.2: Define Variables
PROJECT_NAME=2018_06_05_cmQTL
BATCH_ID=2019_06_10_Batch3
BUCKET=imaging-platform
MAXPROCS=3 # m4.xlarge has 4 cores; keep 1 free
mkdir -p ~/efs/${PROJECT_NAME}/workspace/
cd ~/efs/${PROJECT_NAME}/workspace/
mkdir -p log/${BATCH_ID}
PLATES=$(readlink -f ~/efs/${PROJECT_NAME}/workspace/scratch/${BATCH_ID}/plates_to_process.txt)

############################
# Step 2 - Generate Single Cell Profiles
############################
# I customized the aggregate script select each class of single cells
# Note that this still generates per well profiles!
# Except it will select only certain cells (either in colony or isolated)

# Step 2.1 - Aggregate per well profiles for each single cell category
SC_TYPES=( "colony" "isolated" )
for sc_type in "${SC_TYPES[@]}"
do
  parallel \
    --no-run-if-empty \
    --eta \
    --joblog ../../log/${BATCH_ID}/aggregate_${sc_type}_test.log \
    --results ../../log/${BATCH_ID}/aggregate_${sc_type}_test \
    --files \
    --keep-order \
    ./aggregate.R \
    --sqlite_file /home/ubuntu/bucket/projects/2018_06_05_cmQTL/workspace/backend/${BATCH_ID}/{1}/{1}.sqlite \
    --output ../../backend/${BATCH_ID}/{1}/{1}_${sc_type}.csv \
    --sc_type $sc_type :::: ${PLATES}
done

# Step 2.2 - Annotate the single cell types
for sc_type in "${SC_TYPES[@]}"
do
  parallel \
    --no-run-if-empty \
    --eta \
    --joblog ../../log/${BATCH_ID}/annotate_${sc_type}.log \
    --results ../../log/${BATCH_ID}/annotate_${sc_type} \
    --files \
    --keep-order \
    ./annotate.R \
    --batch_id ${BATCH_ID} \
    --plate_id {1}_${sc_type} :::: ${PLATES}
done

# Step 2.3 - Normalize single cell types
for sc_type in "${SC_TYPES[@]}"
do
  parallel \
    --no-run-if-empty \
    --eta \
    --joblog ../../log/${BATCH_ID}/normalize_${sc_type}.log \
    --results ../../log/${BATCH_ID}/normalize_${sc_type} \
    --files \
    --keep-order \
    ./normalize.R \
    --batch_id ${BATCH_ID} \
    --plate_id {1}_${sc_type} \
    --subset \"Metadata_Well != \'\'\'dummy\'\'\'\" :::: ${PLATES}
done

# Step 2.4 - Apply variable selection defined in `generate-pilot-profiles.sh`
for sc_type in "${SC_TYPES[@]}"
do
  parallel \
    --no-run-if-empty \
    --eta \
    --joblog ../../log/${BATCH_ID}/select_${sc_type}.log \
    --results ../../log/${BATCH_ID}/select_${sc_type} \
    --files \
    --keep-order \
    ./select.R \
    --batch_id ${BATCH_ID} \
    --plate_id {1}_${sc_type} \
    --filters variance_threshold,correlation_threshold,manual :::: ${PLATES}
done
