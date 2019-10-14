#!/usr/bin/env bash

DATASETS=( "AFS20" "AFS20RS90" "AFS31" "AFS31RS90")
PROCESSES_NUMBERS=( "8" "16" "32" "64" "128")

nodes=( "nodo2" "nodo5" "nodo6" "nodo7" "nodo8" "nodo9" "nodo10" "nodo11" "nodo12" "nodo13" "nodo14" "nodo15")

for i in 1 2 3
do
    for DATASET in "${DATASETS[@]}"
    do
        for PROCESSES_NUMBER in "${PROCESSES_NUMBERS[@]}"
        do
            ORDER="rm -Rf /home/josemanuel.abuin/Genomics/Databases/${DATASET}_${PROCESSES_NUMBER}/*"

            for node in "${nodes[@]}"
            do
	              ssh $node "${ORDER}"
            done

           ./BuildGeneric_Local_BD.sh ${DATASET} ${PROCESSES_NUMBER}

        done

    done

done
