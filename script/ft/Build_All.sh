#!/usr/bin/env bash

DATASETS=( "AFS20", "AFS20RS90", "AFS31", "AFS31RS90")
PROCESSES_NUMBERS=( "8" "16" "32" "64" )

DATABASES_PATH="/mnt/lustre/scratch/home/usc/ec/jam/Genomica/Databases/"

for i in 1 2
do
    for DATASET in ${DATASETS[@]}
    do
        for PROCESSES_NUMBER in "${PROCESSES_NUMBERS[@]}"
        do
            echo "rm -Rf /mnt/lustre/scratch/home/usc/ec/jam/Genomica/Databases/${DATASET}_${PROCESSES_NUMBER}/*"

            rm -Rf /mnt/lustre/scratch/home/usc/ec/jam/Genomica/Databases/${DATASET}_${PROCESSES_NUMBER}/*

            sbatch ./Run_Build_${DATASET}_${PROCESSES_NUMBER}.sh

        done

    done

done
