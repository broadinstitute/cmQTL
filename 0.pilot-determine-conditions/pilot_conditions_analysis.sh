#!/bin/bash
#
# cmQTL PILOT processing pipeline
# Gregory Way, 2019
#
# Here, we are using two different plates with various conditions to identify optimal
# conditions for cell painting experiments.
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
#     * Network: vpc-35149752 - Subnet: Default (imaging platform terraform)
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
BATCH_ID=2019_05_13_Batch2
BUCKET=imaging-platform
MAXPROCS=3 # m4.xlarge has 4 cores; keep 1 free
mkdir -p ~/efs/${PROJECT_NAME}/workspace/
cd ~/efs/${PROJECT_NAME}/workspace/
mkdir -p log/${BATCH_ID}
PLATES=$(readlink -f ~/efs/${PROJECT_NAME}/workspace/scratch/${BATCH_ID}/plates_to_process.txt)

# Step 1.3 - Create an EBS temp directory for creating the backend
mkdir ~/ebs_tmp

############################
# Step 2 - Configure Tools to Process Images
############################
cd ~/efs/${PROJECT_NAME}/workspace/
mkdir software
cd software
git clone git@github.com:broadinstitute/cytominer_scripts.git

# Note that additional steps are required if processing profiles directly from images
# Beth Cimini has already done this, so no need to process directly from raw

############################
# Step 3 - Annotate
############################
# NOTE - The annotate step creates `augmented` profiles in the `backend` folder
# `augmented` profiles represent aggregated per-well data annotated with metadata
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
# Note - Variable selection uses both normalized and unnormalized data
mkdir -p ../../parameters/${BATCH_ID}/sample/

# Step 5.0 - Sample normalized and unnormalized data
# Normalized
./sample.R \
  --batch_id ${BATCH_ID} \
  --pattern "_normalized.csv$" \
  --replicates 2 \
  --output ../../parameters/${BATCH_ID}/sample/${BATCH_ID}_normalized_sample.feather

# Unnormalized
./sample.R \
  --batch_id ${BATCH_ID} \
  --pattern "_augmented.csv$" \
  --replicates 2 \
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
SAMPLE_PLATE_ID='BR00103267'
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
