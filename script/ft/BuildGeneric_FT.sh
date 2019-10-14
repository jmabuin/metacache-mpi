#!/bin/bash

#USE: BuildGeneric_DB.sh [Dataset] [Num_procs]

DATASET="AFS20"
INPUT_SEQUENCES="/mnt/lustre/scratch/home/usc/ec/jam/Genomica/AFS20/"
PROCESSES_NUMBER=8
TAXONOMY="/mnt/lustre/scratch/home/usc/ec/jam/Genomica/Taxonomy/"
OUTPUT_DIR="/mnt/lustre/scratch/home/usc/ec/jam/Genomica/Databases/AFS20_${PROCESSES_NUMBER}/DB_AFS20.db"
BIN_PATH="/mnt/lustre/scratch/home/usc/ec/jam/Genomica/metacache-mpi/build/metacache_mpi"

if [[ "$#" -eq 1 ]]; then
    DATASET=$1

elif [[ "$#" -eq 2 ]]; then
    DATASET=$1
    PROCESSES_NUMBER=$2
fi

case "$DATASET" in
        AFS20)
            INPUT_SEQUENCES="/mnt/lustre/scratch/home/usc/ec/jam/Genomica/AFS20RS90/AFS/"
            ;;

        AFS20RS90)
            INPUT_SEQUENCES="/mnt/lustre/scratch/home/usc/ec/jam/Genomica/AFS20RS90/"
            ;;

        AFS31)
            INPUT_SEQUENCES="/mnt/lustre/scratch/home/usc/ec/jam/Genomica/AFS31RS90/AFS/"
            ;;

        AFS31RS90)
            INPUT_SEQUENCES="/mnt/lustre/scratch/home/usc/ec/jam/Genomica/AFS31RS90/"
            ;;

        *)
            INPUT_SEQUENCES="/mnt/lustre/scratch/home/usc/ec/jam/Genomica/AFS20RS90/AFS/"
            ;;

#Set dataset and output dir

OUTPUT_DIR="/mnt/lustre/scratch/home/usc/ec/jam/Genomica/Databases/${DATASET}_${PROCESSES_NUMBER}/DB_${DATASET}.db"

echo "Dataset: ${DATASET}"
echo "Number of processes: $PROCESSES_NUMBER"
echo "Database name: $OUTPUT_DIR"
echo "Input sequences: $INPUT_SEQUENCES"
echo "Taxonomy: $TAXONOMY"

${BIN_PATH} build ${OUTPUT_DIR} ${INPUT_SEQUENCES} -taxonomy ${TAXONOMY} -remove-overpopulated-features
