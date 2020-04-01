#!/bin/bash
#
# cmQTL processing pipeline
# Gregory Way, 2019 (adapted by Shantanu Singh)
#
# Here, we profile all the 480 cell lines plates that were plated under the conditions we found optimal in `0.pilot-determine-conditions`
#
# Instructions to generate cell painting profiles for the cmQTL experiments
# Pipeline generated using the profiling handbook:
#     https://cytomining.github.io/profiling-handbook/

############################
# Step 1 - Configure the Environment
############################
# Step 1.1: Setup a virtual machine on AWS
#
#     * Launch an ec2 instance on AWS
#     * AMI: cytomining/images/hvm-ssd/cytominer-ubuntu-trusty-14.04-amd64-server-1529668435
#     * Instance Type: m4.xlarge
#     * Network: vpc-35149752 - Subnet: subnet-55625c0d (cellpainting_1b)
#     * IAM role: `s3-imaging-platform-role
#     * Add New Volume (if necessary): `EBS` with 110 GiB
#     * No Tags
#     * Select Existing Security Group: `SSH_HTTP`
#     * Review and Launch
#     * ssh -i <USER>.pem ubuntu@<Public DNS IPv4>
#     * Inside AWS terminal: `aws configure` and input security credentials
#
# See https://cytomining.github.io/profiling-handbook/configure-environment.html#set-up-a-virtual-machine
# for more details

# Step 1.2: Define Variables
PROJECT_NAME=2018_06_05_cmQTL
BATCH_ID=2019_06_10_Batch3
BUCKET=imaging-platform
MAXPROCS=3 # m4.xlarge has 4 cores; keep 1 free
mkdir -p ~/efs/${PROJECT_NAME}/workspace/
cd ~/efs/${PROJECT_NAME}/workspace/
mkdir -p log/${BATCH_ID}

# the scripts were run verbatim for these two additional batches with appropriate modifications to 
# SAMPLE_PLATE_ID and contents of $PLATES
#BATCH_ID=2019_08_15_Batch4 
#BATCH_ID=2019_09_06_Batch5
#BATCH_ID=2020_03_05_Batch6

# Step 1.3 - Create an EBS temp directory for creating the backend
mkdir ~/ebs_tmp

############################
# Step 2 - Configure Tools to Process Images
############################
cd ~/efs/${PROJECT_NAME}/workspace/
mkdir -p software
cd software

# check whether the repos below already exist
# If they do, run `git status` to make sure that the files haven't been modified
# If they don't exist then clone as below
git clone git@github.com:broadinstitute/cellpainting_scripts.git
git clone git@github.com:broadinstitute/pe2loaddata.git
git clone git@github.com:broadinstitute/cytominer_scripts.git
git clone git@github.com:CellProfiler/Distributed-CellProfiler.git

# Copy config_files/config.yml to cellpainting_scripts/ (overwrite the existing file)
      
# Follow these steps verbatim
# https://cytomining.github.io/profiling-handbook/configure-tools-to-process-images.html#setup-distributed-cellprofiler
# https://cytomining.github.io/profiling-handbook/setup-pipelines-and-images.html#get-cellprofiler-pipelines

# Specify pipeline set
PIPELINE_SET=cellpainting_ipsc_20x_phenix_with_bf_bin1_cp3/

# Follow these steps verbatim
# https://cytomining.github.io/profiling-handbook/setup-pipelines-and-images.html#prepare-images

# Create list of plates
mkdir -p scratch/${BATCH_ID}/
PLATES=$(readlink -f ~/efs/${PROJECT_NAME}/workspace/scratch/${BATCH_ID}/plates_to_process.txt)
echo "cmqtlpl1.5-31-2019-mt cmqtlpl261-2019-mt"|tr " " "\n" > ${PLATES}

# Follow these steps verbatim
# https://cytomining.github.io/profiling-handbook/setup-pipelines-and-images.html#create-loaddata-csvs
# Note: the cellpainting_scripts/config.ini file had not been set up correctly, and pe2loaddata got stuck
# so I did it by hand (i.e. `parallel --dry-run` and then run each command)

# Follow these steps verbatim
# https://cytomining.github.io/profiling-handbook/setup-jobs.html#illumination-correction
# https://cytomining.github.io/profiling-handbook/setup-jobs.html#analysis
# but change the docker image for both 
# `--cp_docker_image cellprofiler/cellprofiler:3.1.8`
# This errored but the config files got created without a hitch. See https://github.com/broadinstitute/cmQTL/issues/14#issuecomment-505551405

# Copy the `dcp_config_files` directory to `cellpainting_scripts`
cp ../2018_06_05_cmQTL/1.profile-cell-lines/dcp_config_files/* ../cellpainting_scripts/dcp_config_files/

# Follow these steps verbatim 
# https://cytomining.github.io/profiling-handbook/run-jobs.html#run-illum-dcp
# https://cytomining.github.io/profiling-handbook/run-jobs.html#dcp
# https://cytomining.github.io/profiling-handbook/create-profiles.html#create-database-backend

############################
# Step 3 - Annotate
############################
# NOTE - The annotate step creates `augmented` profiles in the `backend` folder
# `augmented` profiles represent aggregated per-well data annotated with metadata
# Follow this step with some modifications (see below)
# https://cytomining.github.io/profiling-handbook/create-profiles.html#annotate

# Retrieve metadata information
aws s3 sync s3://${BUCKET}/projects/${PROJECT_NAME}/workspace/metadata/${BATCH_ID}/ ~/efs/${PROJECT_NAME}/workspace/metadata/${BATCH_ID}/

# Use cytominer_scripts to run annotation
cd  ~/efs/${PROJECT_NAME}/workspace/software/cytominer_scripts

parallel \
  --no-run-if-empty \
  --eta \
  --joblog ../../log/${BATCH_ID}/annotate.log \
  --results ../../log/${BATCH_ID}/annotate \
  --files \
  --keep-order \
  ./annotate.R \
  --batch_id ${BATCH_ID} \
  --plate_id {1} :::: ${PLATES}

############################
# Step 4 - Normalize
############################
# Note - The normalize step creates `normalized` profiles in the `backend` folder
# The step z-scores each feature using all wells (i.e. use all "non-dummy" wells)

# Follow this step with some modifications (see below)
# https://cytomining.github.io/profiling-handbook/create-profiles.html#normalize

parallel \
  --no-run-if-empty \
  --eta \
  --joblog ../../log/${BATCH_ID}/normalize.log \
  --results ../../log/${BATCH_ID}/normalize \
  --files \
  --keep-order \
  ./normalize.R \
  --batch_id ${BATCH_ID} \
  --plate_id {1} \
  --subset \"Metadata_Well != \'\'\'dummy\'\'\'\" :::: ${PLATES}

############################
# Step 5 - Variable Selection
############################
# Follow this step with some modifications (no replicate correlation step; see below)
# https://cytomining.github.io/profiling-handbook/create-profiles.html#select-variables

# Note - Variable selection uses both normalized and unnormalized data
mkdir -p ../../parameters/${BATCH_ID}/sample/

# In the sampling steps below, we specify replicates = 1 because there is only 1 
# replicate plate per platemap. Currently, sampling replicates within a plate is
# not supported in cytominer_scripts

# Step 5.0 - Sample normalized and unnormalized data
# Normalized
./sample.R \
  --batch_id ${BATCH_ID} \
  --pattern "_normalized.csv$" \
  --replicates 1 \
  --output ../../parameters/${BATCH_ID}/sample/${BATCH_ID}_normalized_sample.feather

# Unnormalized
./sample.R \
  --batch_id ${BATCH_ID} \
  --pattern "_augmented.csv$" \
  --replicates 1 \
  --output ../../parameters/${BATCH_ID}/sample/${BATCH_ID}_augmented_sample.feather

# Using the sampled feather files, perform a series of three variable selection steps
# Step 5.1 - Remove variables that have high correlations with other variables
./preselect.R \
  --batch_id ${BATCH_ID} \
  --input ../../parameters/${BATCH_ID}/sample/${BATCH_ID}_normalized_sample.feather \
  --operations correlation_threshold

# Step 5.2 - Remove variables that have low variance
./preselect.R \
  --batch_id ${BATCH_ID} \
  --input ../../parameters/${BATCH_ID}/sample/${BATCH_ID}_augmented_sample.feather \
  --operations variance_threshold

# Step 5.3 - Remove features known to be noisy
SAMPLE_PLATE_ID='cmqtlpl1.5-31-2019-mt'
echo "variable" > ../../parameters/${BATCH_ID}/variable_selection/manual.txt

head -1 \
  ../../backend/${BATCH_ID}/${SAMPLE_PLATE_ID}/${SAMPLE_PLATE_ID}.csv \
  |tr "," "\n"|grep -v Meta|grep -E -v 'Granularity_14|Granularity_15|Granularity_16|Manders|RWC' >> \
  ../../parameters/${BATCH_ID}/variable_selection/manual.txt

# Step 5.4 - Apply the variable selection steps to the profiles
# Note - This creates the _normalized_variable_selected.csv files in `backend`
parallel \
  --no-run-if-empty \
  --eta \
  --joblog ../../log/${BATCH_ID}/select.log \
  --results ../../log/${BATCH_ID}/select \
  --files \
  --keep-order \
  ./select.R \
  --batch_id ${BATCH_ID} \
  --plate_id {1} \
  --filters variance_threshold,correlation_threshold,manual :::: ${PLATES}

############################
# Compute cell counts
############################

# If jumping in directly to this step, define BATCH_ID and PLATES
# BATCH_ID=2019_06_10_Batch3
# BATCH_ID=2019_08_15_Batch4 
# BATCH_ID=2019_09_06_Batch5
# PLATES=$(readlink -f ~/efs/${PROJECT_NAME}/workspace/scratch/${BATCH_ID}/plates_to_process.txt)

cd ~/efs/${PROJECT_NAME}/workspace/software/cytominer_scripts

mkdir -p ~/ebs_tmp/2018_06_05_cmQTL/workspace/backend/${BATCH_ID}

parallel \
  --no-run-if-empty \
  --eta \
  aws s3 sync \
  --exclude '"*"' \
  --include '"*.sqlite"' \
  s3://imaging-platform/projects/2018_06_05_cmQTL/workspace/backend/${BATCH_ID}/{1}/ \
  ~/ebs_tmp/2018_06_05_cmQTL/workspace/backend/${BATCH_ID}/{1}/ \
  :::: ${PLATES}

# cell counts
parallel \
  --no-run-if-empty \
  --max-procs 1 \
  --joblog ../../log/${BATCH_ID}/stats.log \
  --results ../../log/${BATCH_ID}/stats \
  --files \
  --keep-order \
  --eta \
  ./stats.R  \
  ~/ebs_tmp/2018_06_05_cmQTL/workspace/backend/${BATCH_ID}/{1}/{1}.sqlite \
  -o ../../backend/${BATCH_ID}/{1}/{1}_count.csv \
  :::: ${PLATES}

# check rows
parallel \
  --no-run-if-empty \
  --keep-order \
  wc -l ../../backend/${BATCH_ID}/{1}/{1}_count.csv :::: \
  ${PLATES}

# remove sqlite
parallel \
  --no-run-if-empty \
  --eta \
  rm ~/ebs_tmp/2018_06_05_cmQTL/workspace/backend/${BATCH_ID}/{1}/{1}.sqlite \
  :::: ${PLATES}

############################
# Step 6 - Audit
############################

# Audit for replicate reproducibility
# https://cytomining.github.io/profiling-handbook/create-profiles.html#audit

PLATE_MAPS=../../scratch/${BATCH_ID}/plate_maps.txt

csvcut -c Plate_Map_Name \
  ../../metadata/${BATCH_ID}/barcode_platemap.csv | \
  tail -n +2|sort|uniq > \
  ${PLATE_MAPS}
  
mkdir -p ../../audit/${BATCH_ID}/

parallel \
  --no-run-if-empty \
  --eta \
  --joblog ../../log/${BATCH_ID}/audit.log \
  --results ../../log/${BATCH_ID}/audit \
  --files \
  --keep-order \
  ./audit.R \
  -b ${BATCH_ID} \
  -m {1} \
  -f _normalized_variable_selected.csv \
  -o ../../audit/${BATCH_ID}/{1}_audit.csv \
  -l ../../audit/${BATCH_ID}/{1}_audit_detailed.csv \
  -p Metadata_Plate_Map_Name,Metadata_line_ID,Metadata_plating_density :::: ${PLATE_MAPS}
  
############################
# Step 6 - Convert to other formats
############################

# Follow only the first part of this step 
# https://cytomining.github.io/profiling-handbook/create-profiles.html#convert-to-other-formats
# i.e.:

parallel \
  --no-run-if-empty \
  --eta \
  --joblog ../../log/${BATCH_ID}/csv2gct_backend.log \
  --results ../../log/${BATCH_ID}/csv2gct_backend \
  --files \
  --keep-order \
  ./csv2gct.R \
  ../../backend/${BATCH_ID}/{1}/{1}_{2}.csv \
  -o ../../backend/${BATCH_ID}/{1}/{1}_{2}.gct :::: ${PLATES} ::: augmented normalized normalized_variable_selected
  
############################
# Step 7 - Upload data
############################

# Follow this step nearly verbatim 
# https://cytomining.github.io/profiling-handbook/create-profiles.html#upload-data

parallel \
  aws s3 sync \
  ../../{1}/${BATCH_ID}/ \
  s3://${BUCKET}/projects/${PROJECT_NAME}/workspace/{1}/${BATCH_ID}/ ::: audit backend batchfiles  load_data_csv log metadata parameters scratch

