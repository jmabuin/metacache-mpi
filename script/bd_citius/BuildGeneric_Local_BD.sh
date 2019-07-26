#!/bin/bash

#USE: BuildGeneric_DB.sh [Dataset] [Num_procs]

DATASET="AFS20"
INPUT_SEQUENCES="/home/josemanuel.abuin/Genomics/AFS20/"
PROCESSES_NUMBER=8
TAXONOMY="/home/josemanuel.abuin/Genomics/Taxonomy/"
OUTPUT_DIR="/home/josemanuel.abuin/Genomics/Databases/AFS20_${PROCESSES_NUMBER}/DB_AFS20.db"
NODES_FILE="/home/remoto/josemanuel.abuin/nodefilesmall"


if [ "$#" -eq 1 ]; then
    DATASET=$1

elif [ "$#" -eq 2 ]; then
    DATASET=$1
    PROCESSES_NUMBER=$2
fi


case "$PROCESSES_NUMBER" in
        4)
            NODES_FILE="/home/remoto/josemanuel.abuin/nodefileverysmall"
            ;;

        8)
            NODES_FILE="/home/remoto/josemanuel.abuin/nodefilesmall"
            ;;

        128)
            NODES_FILE="/home/remoto/josemanuel.abuin/nodefilebig"
            ;;

        *)
            NODES_FILE="/home/remoto/josemanuel.abuin/nodefile"
            ;;

esac


#Set dataset and output dir
INPUT_SEQUENCES="/home/josemanuel.abuin/Genomics/${DATASET}/"
OUTPUT_DIR="/home/josemanuel.abuin/Genomics/Databases/${DATASET}_${PROCESSES_NUMBER}/DB_${DATASET}.db"

echo "Dataset: ${DATASET}"
echo "Number of processes: $PROCESSES_NUMBER"
echo "Nodefile: $NODES_FILE"
echo "Database name: $OUTPUT_DIR"
echo "Input sequences: $INPUT_SEQUENCES"
echo "Taxonomy: $TAXONOMY"

mpiexec -mca btl_tcp_if_include em1 -np $PROCESSES_NUMBER -hostfile $NODES_FILE /home/josemanuel.abuin/Genomics/metacache-mpi/build/metacache_mpi build $OUTPUT_DIR $INPUT_SEQUENCES -taxonomy $TAXONOMY -remove-overpopulated-features
