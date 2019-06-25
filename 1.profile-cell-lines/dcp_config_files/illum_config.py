# Constants (User configurable)
# See https://github.com/CellProfiler/Distributed-CellProfiler/wiki/Step-1%3A--Configuration

APP_NAME = '2018_06_05_cmQTL_Illum'                # Used to generate derivative names unique to the application.

# DOCKER REGISTRY INFORMATION:
DOCKERHUB_TAG = 'bethcimini/distributed-cellprofiler:1.2.0_3.1.8'

# AWS GENERAL SETTINGS:
AWS_REGION = 'us-east-1'
AWS_PROFILE = 'default'                 # The same profile used by your AWS CLI installation
SSH_KEY_NAME = 'CellProfiler.pem'      # Expected to be in ~/.ssh
AWS_BUCKET = 'imaging-platform'

# EC2 AND ECS INFORMATION:
ECS_CLUSTER = 'Beth'
CLUSTER_MACHINES = 2
TASKS_PER_MACHINE = 1
EBS_VOL_SIZE = 22                       # In GB.  Minimum allowed is 22.

# DOCKER INSTANCE RUNNING ENVIRONMENT:
DOCKER_CORES = 1                        # Number of CellProfiler processes to run inside a docker container
CPU_SHARES = DOCKER_CORES * 1024        # ECS computing units assigned to each docker container (1024 units = 1 core)
MEMORY = 15000                          # Memory assigned to the docker container in MB
SECONDS_TO_START = 5*60                 # Wait before the next CP process is initiated to avoid memory collisions
DOCKER_BASE_SIZE = 20                   # In GB; default is 10.  The amount of hard disk space each docker container uses.
                                        # EBS_VOL_SIZE should be >= DOCKER_BASE_SIZE * TASKS_PER_MACHINE

# SQS QUEUE INFORMATION:
SQS_QUEUE_NAME = APP_NAME + 'Queue'
SQS_MESSAGE_VISIBILITY = 720*60           # Timeout (secs) for messages in flight (average time to be processed)
SQS_DEAD_LETTER_QUEUE = 'arn:aws:sqs:us-east-1:385009899373:DeadMessages_Beth'

# LOG GROUP INFORMATION:
LOG_GROUP_NAME = APP_NAME

# REDUNDANCY CHECKS
CHECK_IF_DONE_BOOL = 'False'  #True or False- should it check if there are a certain number of non-empty files and delete the job if yes?
EXPECTED_NUMBER_FILES = 7    #What is the number of files that trigger skipping a job?
