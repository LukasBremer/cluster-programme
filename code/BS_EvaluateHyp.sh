#!/bin/bash
# This script will be executed by bash on the backend HPC computers
# It can also run on Sirona and there should be no difference
# call it on sirona by bash SCRIPTNAME.sh and submit to cluster
# queue by submission script or qsub [qsub arguments] SCRIPTNAME.sh.

usage () {
        cat <<HELP_USAGE
        $(basename $0) OutputFolder InputFolder ResultFolder HyperparameterFolder
        args:
        1 NumericalExperimentOutputFolder: In this directory all output from the NumericalExperiment will be written
        2 NumericalExperimentInput: Folder containing Input for the numerical experiment.
        3 Subfolder of ExperimentOutput to store Data and Hyperparameter Configs to take Parameters from.
        4 Number of HyperparameterSets with same Transient Frames.
        5 Model for which prediction is performed.
HELP_USAGE
}


# Argument checking
if [[ $# == 0 ]]
    then
        usage
        echo "No arguments provided. Aborting."
    exit
fi


source ~/.bashrc
# usually for a experiment you want to have a certain environment, for example from anaconda
# or pyenv or virtualenv or something else... I use anaconda but I am not sure this is the best
conda activate # baltasar hatte hier: experiment_environment. Ich hab nur eine

export OMP_NUM_THREADS=1

# Store the date when this script was executed to DATE variable:
DATE=`date +"%Y-%m-%d"`

# I defined the first argument to be the Simulation output folder 
# The Syntax :-ALTERNATIVE set the variable 'NumericalExperimentOutputFolder' to ALTERNATIVE
# if it was not given when called the script (-> $1 is empty)
NumericalExperimentOutputFolder=${1:-SimulationOutput/${DATE}_default}
NumericalExperimentInputFolder=${2:-SimulationInput/${DATE}_default}
# Make use of Cluster Queue Related Variables:
# [Background: a cluster job has a unique identifier, the JOB_ID
# this is stored in the variable name $JOB_ID by the cluster
# it also stores the name of the job given on submission in the
# JOB_NAME variable and -- if an array job (one submission with
# subtasks) is performed -- it stores an identifier of the 
# task in SGE_TASK_ID.

TASK_ID="${SGE_TASK_ID:-1}"
TASK_ID_ZEROED=$(expr $TASK_ID - 1) # start at 0 not at 

NAME="${JOB_NAME:-Reservoir}"
#NAME="${TrainingsData}"
MY_JOB_ID="${JOB_ID:-00000000}"

# SIMID: a unique identifier [note: not 100% unique because the different
# queues can produce similar SIMIDs as the JOB_ID is only a running integer,
# would be better if queue would be represented in number]


# Output Directory handling

if [ -d "$NumericalExperimentOutputFolder" ]; then
    echo "$NumericalExperimentOutputFolder existed. "
    echo "Writing to the same directory."
else
    mkdir -p $NumericalExperimentOutputFolder
    echo "Created $NumericalExperimentOutputFolder."
fi


FILE_NR=$((${TASK_ID_ZEROED} * ${4}))

#Create Config Input Directory
SUB_CONFIG_DIR="/ExperimentOutput_scratch/$5/$3/"
CONFIG_DIR="${NumericalExperimentInputFolder}HyperparameterConfigs/$SUB_CONFIG_DIR/Parameter${FILE_NR}.yml"


if [ -f "$CONFIG_DIR" ]; then
    echo "$CONFIG_DIR exists and will be used."
else
    echo "$CONFIG_DIR does not exist. Exiting. "
    exit
fi

#Resetting Transient Frames for all Hyperparameterfiles [${TASK_ID_ZEROED} * ${4}, ${TASK_ID_ZEROED} * ${4} + 1 ,..., (${TASK_ID_ZEROED} + 1) * ${4}]

#echo "python3 SetTransientFrames.py ${2} ${4} ${FILE_NR} ${4}"
#python3 SetTransientFrames.py ${2} ${4} ${FILE_NR} ${4}

for ((i=0;i<${4};i++))
do
    FILE_NR=$((${i} + ${TASK_ID_ZEROED} * ${4}))
    
    # Create Task OutputFolder
    SIMID="${MY_JOB_ID}_${NAME}_${FILE_NR}"
    SIM_OUTPUT_DIR="${NumericalExperimentOutputFolder}ExperimentOutput_scratch/$5/$3/${SIMID}/"
    if [ -d "$SIM_OUTPUT_DIR" ]; then
        echo "$SIM_OUTPUT_DIR exists. Writing to the same Directory."
    else
        mkdir -p $SIM_OUTPUT_DIR
        echo "Created $SIM_OUTPUT_DIR."
    fi
    # Calling Function
    echo "python3 EvaluateHyperparameter.py $SIM_OUTPUT_DIR $2 $SUB_CONFIG_DIR ${FILE_NR} "
    python3 EvaluateHyperparameter.py 0 $SIM_OUTPUT_DIR $2 $SUB_CONFIG_DIR ${FILE_NR}

done


