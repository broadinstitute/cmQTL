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
# NOTE: We needed at 64Gb instance (r5.2xlarge) for this to complete successfully
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
# Step 2 - Generate single cell profiles
############################
# I customized the aggregate script select each class of single cells
# Note that this still generates per well profiles!
# Except it will select only certain cells (either in colony or isolated)

# Step 2.1 - Aggregate per well profiles for each single cell category

cd ~/efs/${PROJECT_NAME}/workspace/software/cytominer_scripts

mkdir -p ~/ebs_tmp/2018_06_05_cmQTL/workspace/backend/${BATCH_ID}

parallel \
  --no-run-if-empty \
  --eta \
  aws s3 sync  \
  s3://imaging-platform/projects/2018_06_05_cmQTL/workspace/backend/${BATCH_ID}/{1}/ \
  ~/ebs_tmp/2018_06_05_cmQTL/workspace/backend/${BATCH_ID}/{1}/ \
  :::: ${PLATES}

    
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
    --sqlite_file ~/ebs_tmp/2018_06_05_cmQTL/workspace/backend/${BATCH_ID}/{1}/{1}.sqlite \
    --output ../../backend/${BATCH_ID}/{1}/{1}_${sc_type}.csv \
    --sc_type $sc_type :::: ${PLATES}
done

# check rows
for sc_type in "${SC_TYPES[@]}"
do
  parallel \
    --no-run-if-empty \
    --keep-order \
    wc -l ../../backend/${BATCH_ID}/{1}/{1}_${sc_type}.csv :::: ${PLATES}
done  

############################
# Step 2.2 - Annotate the single cell types
############################

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

# check rows
for sc_type in "${SC_TYPES[@]}"
do
  parallel \
    --no-run-if-empty \
    --keep-order \
    wc -l ../../backend/${BATCH_ID}/{1}/{1}_${sc_type}_augmented.csv :::: ${PLATES}
done  

############################
# Step 2.3 - Normalize single cell types
############################

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

# check rows
for sc_type in "${SC_TYPES[@]}"
do
  parallel \
    --no-run-if-empty \
    --keep-order \
    wc -l ../../backend/${BATCH_ID}/{1}/{1}_${sc_type}_normalized.csv :::: ${PLATES}
done

############################
# Step 2.4 - Apply variable selection 
############################

# Save previously generated variable selection files that were generated
# based on the regular (all cells) profiles

mv \
  ../../parameters/${BATCH_ID}/variable_selection/correlation_threshold.txt \
  ../../parameters/${BATCH_ID}/variable_selection/correlation_threshold_all.txt
  
mv \
  ../../parameters/${BATCH_ID}/variable_selection/variance_threshold.txt \
  ../../parameters/${BATCH_ID}/variable_selection/variance_threshold_all.txt

# Sample normalized and unnormalized data
# Normalized
for sc_type in "${SC_TYPES[@]}"
do
  ./sample.R \
    --batch_id ${BATCH_ID} \
    --pattern "_${sc_type}_normalized.csv$" \
    --replicates 1 \
    --output ../../parameters/${BATCH_ID}/sample/${BATCH_ID}_${sc_type}_normalized_sample.feather
done

# Unnormalized
for sc_type in "${SC_TYPES[@]}"
do
  ./sample.R \
    --batch_id ${BATCH_ID} \
    --pattern "_${sc_type}_augmented.csv$" \
    --replicates 1 \
    --output ../../parameters/${BATCH_ID}/sample/${BATCH_ID}_${sc_type}_augmented_sample.feather
done  
  
# Using the sampled feather files, perform a series of three variable selection steps

# Remove variables that have high correlations with other variables
for sc_type in "${SC_TYPES[@]}"
do
  ./preselect.R \
    --batch_id ${BATCH_ID} \
    --input ../../parameters/${BATCH_ID}/sample/${BATCH_ID}_${sc_type}_normalized_sample.feather \
    --operations correlation_threshold
    
  cp \
    ../../parameters/${BATCH_ID}/variable_selection/correlation_threshold.txt \
    ../../parameters/${BATCH_ID}/variable_selection/correlation_threshold_${sc_type}.txt

done

# Remove variables that have low variance
for sc_type in "${SC_TYPES[@]}"
do
  ./preselect.R \
    --batch_id ${BATCH_ID} \
    --input ../../parameters/${BATCH_ID}/sample/${BATCH_ID}_${sc_type}_augmented_sample.feather \
    --operations variance_threshold
    
  cp \
    ../../parameters/${BATCH_ID}/variable_selection/variance_threshold.txt \
    ../../parameters/${BATCH_ID}/variable_selection/variance_threshold_${sc_type}.txt
  
done

# Apply the variable selection steps to the profiles
# Note - This creates the _${sc_type}_normalized_variable_selected.csv files in `backend`

for sc_type in "${SC_TYPES[@]}"
do
  cp \
    ../../parameters/${BATCH_ID}/variable_selection/variance_threshold_${sc_type}.txt \
    ../../parameters/${BATCH_ID}/variable_selection/variance_threshold.txt

  cp \
    ../../parameters/${BATCH_ID}/variable_selection/correlation_threshold_${sc_type}.txt \
    ../../parameters/${BATCH_ID}/variable_selection/correlation_threshold.txt

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

  rm \
    ../../parameters/${BATCH_ID}/variable_selection/correlation_threshold.txt
  
  rm \
    ../../parameters/${BATCH_ID}/variable_selection/variance_threshold.txt

done

parallel \
  --no-run-if-empty \
  --keep-order \
  wc -l ../../backend/${BATCH_ID}/{1}/{1}_${sc_type}_normalized_variable_selected.csv :::: ${PLATES}

############################
# Step 6 - Audit
############################

# Adapt this step 
# https://cytomining.github.io/profiling-handbook/create-profiles.html#audit
# i.e.:

PLATE_MAPS=../../scratch/${BATCH_ID}/plate_maps.txt

csvcut -c Plate_Map_Name \
  ../../metadata/${BATCH_ID}/barcode_platemap.csv | \
  tail -n +2|sort|uniq > \
  ${PLATE_MAPS}
  
mkdir -p ../../audit/${BATCH_ID}/

for sc_type in "${SC_TYPES[@]}"
do
  parallel \
    --no-run-if-empty \
    --eta \
    --joblog ../../log/${BATCH_ID}/audit_${sc_type}.log \
    --results ../../log/${BATCH_ID}/audit_${sc_type} \
    --files \
    --keep-order \
    ./audit.R \
    -b ${BATCH_ID} \
    -m {1} \
    -f _${sc_type}_normalized_variable_selected.csv \
    -o ../../audit/${BATCH_ID}/{1}_audit_${sc_type}.csv \
    -l ../../audit/${BATCH_ID}/{1}_audit_${sc_type}_detailed.csv \
    -p Metadata_Plate_Map_Name,Metadata_line_ID,Metadata_plating_density :::: ${PLATE_MAPS}
done

############################
# Step 6 - Convert to other formats
############################

# Adapt this step 
# https://cytomining.github.io/profiling-handbook/create-profiles.html#convert-to-other-formats
# i.e.:

for sc_type in "${SC_TYPES[@]}"
do
  parallel \
    --no-run-if-empty \
    --eta \
    --joblog ../../log/${BATCH_ID}/csv2gct_backend.log \
    --results ../../log/${BATCH_ID}/csv2gct_backend \
    --files \
    --keep-order \
    ./csv2gct.R \
    ../../backend/${BATCH_ID}/{1}/{1}_{2}.csv \
    -o ../../backend/${BATCH_ID}/{1}/{1}_{2}.gct :::: ${PLATES} ::: ${sc_type}_augmented ${sc_type}_normalized ${sc_type}_normalized_variable_selected
done

############################
# Step 7 - Upload data
############################

# Adapt this step 
# https://cytomining.github.io/profiling-handbook/create-profiles.html#upload-data
# i.e.:

parallel \
  aws s3 sync \
  ../../{1}/${BATCH_ID}/ \
  s3://${BUCKET}/projects/${PROJECT_NAME}/workspace/{1}/${BATCH_ID}/ ::: audit backend log scratch

PROJECT_NAME=2018_06_05_cmQTL
BATCH_ID=2019_06_10_Batch3
BUCKET=imaging-platform

cd ~/work/projects/${PROJECT_NAME}/workspace

aws s3 sync \
  s3://${BUCKET}/projects/${PROJECT_NAME}/workspace/audit/${BATCH_ID}/ \
  audit/${BATCH_ID}/
  
aws s3 sync --exclude "*.sqlite" \
  s3://${BUCKET}/projects/${PROJECT_NAME}/workspace/audit/${BATCH_ID}/ \
  backend/${BATCH_ID}/
