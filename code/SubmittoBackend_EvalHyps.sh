#!/bin/bash
# This script is used to submit a backend script.

## Explanation of qsub params also via man qsub
# Use bash as shell 
# -S /bin/bash
# -V # Preserve environment variables
# -cwd # Execute from current working directory
# -j yes # Merge standard output and standard error into one file
# -N ms_fem_complete_run # Name of job
# -o /data.bmp/$USER/q-out/ # path for the output files
# Limit memory usage # maximales memory: anzahl kerne * 6 gb - 0.5 gb (puffer fuer alle faelle):
# -hard -l h_vmem=5.8G
# -t 1-30 # array range for array jobs.
# -tc 20 # num tasks at once
# which queue to use
# -q mvapich2-grannus.q # mvapich2 is the parallel q

usage () {
        cat <<HELP_USAGE
        Submit script to cluster. Usage:
        $(basename $0) meshing_paramater_file parallel
        args:
        1 NumericalExperimentOuputFolder: Output folder for this experiment (pathto bmp /data.bmp/lfleddermann/SimulationData)
        2 NumericalExperimentInput: Input for this experiment. Each task will
        take its input from this folder
        3 Name, OutputSubfolder, InputSubfolder: the experiments name (Name of Folders and Name of submition order)
        4 Name of Model for which prediction is performed.

HELP_USAGE
}

# change this variable to the project you are using
PROJECT_NAME="Lorenz63Reservoir"


# Argument checking
if [[ $# != 4 ]]
    then
        usage
        echo "No arguments provided. Aborting."
    exit
fi


# Handle parallel
# The Backend Script is not written for parallel jobs
parallel=''

if test -z "$parallel"; then
    queue="-q teutates.q" 
    MEMINGB=5.8
else -hard -l h_vmem="${MEMINGB}G" \
        -t $ARRAYSTART-$ARRAYSTOP\
        -tc 20\
        $queue\
    # the -pe option is needed for  parallel computing
    cores=32
    queue="-pe mvapich2-teutates $cores"
    MEMINGB=$(expr $cores \* 4 - 1)
fi


NAME=$3
DATE=`date +"%Y-%m-%d"`
YEAR=`date +"%Y"`
# Folder to save output of each task
OUTPUT_PATH="${1}/ClusterOutput_scratch/${4}/${3}"
EX_OUTPUT_PATH="${1}/ExperimentOutput_scratch/${4}/${3}"


echo "Saving cluster output to $OUTPUT_PATH."
mkdir -p $OUTPUT_PATH
echo "Saving experiment output to $EX_OUTPUT_PATH."
mkdir -p $EX_OUTPUT_PATH

#Checking If config directory exists
CONFIG_DIR="${2}/HyperparameterConfigs/ExperimentOutput_scratch/${4}/${3}/"
if [ -d "$CONFIG_DIR" ]; then
    echo "$CONFIG_DIR exists and will be used."
else
    echo "$CONFIG_DIR does not exist. Exiting. "
    exit
fi
cp -a "${CONFIG_DIR}/MetaParameter.yml" ${EX_OUTPUT_PATH}


ARRAYSTART=1
# Either add number of different experiments by hand # can be found in MetaFile
# or caculate in bash (google helps you with it :) )

TOTALSETS=1247400
SETSPERJOB=50
ARRAYSTOP=$(expr ${TOTALSETS} / ${SETSPERJOB}) #For this purpise it should equal the whole state space/ all quantities independend of the TransientFrameNr


qsub -S "/bin/bash" \
    -V \
    -cwd \
    -j "yes" \
    -N $NAME \
    -o $OUTPUT_PATH\
    -hard -l h_vmem="${MEMINGB}G" \
    -t $ARRAYSTART-$ARRAYSTOP\
    -tc 1500\
    $queue\
    BS_EvaluateHyp.sh $1 $2 $3 ${SETSPERJOB} $4


# BS_EvaluateHyp.sh takes the following inputs
# 1 NumericalExperimentOutputFolder: In this directory all output from the NumericalExperiment will be written
# 2 NumericalExperimentInput: Folder containing Input for the numerical experiment.
# 3 Subfolder of ExperimentOutput to store Data and Hyperparameter Configs to take Parameters from.
# 4 Number of HyperparameterSets with same Transient Frames.
# 5 Model for which prediction is performed.
