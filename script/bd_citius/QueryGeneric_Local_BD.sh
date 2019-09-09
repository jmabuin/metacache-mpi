#!/usr/bin/env bash

#USE: QueryGeneric_DB.sh [Dataset] [Input dataset] [Num_procs] [Num_Threads] [Query limit]
#Example: QueryGeneric_DB.sh AFS20 S12 4 2 1000000

DATASET="AFS20"
INPUT_DATASET="S4"
INPUT_SEQUENCES="/media/glusterfs/genomics/Inputs/fastq/113003602_S4_forward_paired.fq /media/glusterfs/genomics/Inputs/fastq/113003602_S4_reverse_paired.fq"
PROCESSES_NUMBER=8
THREADS_NUMBER=2
DATABASE_NAME="/home/josemanuel.abuin/Genomics/Databases/${DATASET}_${PROCESSES_NUMBER}/DB_${DATASET}.db"
OUTPUT_DIR="/home/josemanuel.abuin/Genomics/Outputs/Output_${DATASET}_${INPUT_DATASET}_${PROCESSES_NUMBER}P_${THREADS_NUMBER}T.txt"
NODES_FILE="/home/remoto/josemanuel.abuin/nodefilesmall"
QUERY_LIMIT=25000

if [[ "$#" -eq 1 ]]; then
    DATASET=$1

elif [[ "$#" -eq 2 ]]; then
    DATASET=$1
    INPUT_DATASET=$2

elif [[ "$#" -eq 3 ]]; then
    DATASET=$1
    INPUT_DATASET=$2
    PROCESSES_NUMBER=$3

elif [[ "$#" -eq 4 ]]; then
    DATASET=$1
    INPUT_DATASET=$2
    PROCESSES_NUMBER=$3
    THREADS_NUMBER=$4

elif [[ "$#" -eq 5 ]]; then
    DATASET=$1
    INPUT_DATASET=$2
    PROCESSES_NUMBER=$3
    THREADS_NUMBER=$4
    QUERY_LIMIT=$5

fi


case "$THREADS_NUMBER" in
        4)
            NODES_FILE="/home/remoto/josemanuel.abuin/nodefilesmall"
            ;;

        8)
            NODES_FILE="/home/remoto/josemanuel.abuin/nodefileverysmall"
            ;;

        16)
            NODES_FILE="/home/remoto/josemanuel.abuin/nodefileverysmall"
            ;;

        32)
            NODES_FILE="/home/remoto/josemanuel.abuin/nodefiletiny"
            ;;

        *)
            NODES_FILE="/home/remoto/josemanuel.abuin/nodefile"
            ;;

esac

case "$INPUT_DATASET" in
        S4)
            INPUT_SEQUENCES="/media/glusterfs/genomics/Inputs/fastq/113003602_S4_forward_paired.fq /media/glusterfs/genomics/Inputs/fastq/113003602_S4_reverse_paired.fq"
            ;;

        S5)
            INPUT_SEQUENCES="/media/glusterfs/genomics/Inputs/fastq/113003610_S5_forward_paired.fq /media/glusterfs/genomics/Inputs/fastq/113003610_S5_reverse_paired.fq"
            ;;

        S6)
            INPUT_SEQUENCES="/media/glusterfs/genomics/Inputs/fastq/113003628_S6_forward_paired.fq /media/glusterfs/genomics/Inputs/fastq/113003628_S6_reverse_paired.fq"
            ;;

        S7)
            INPUT_SEQUENCES="/media/glusterfs/genomics/Inputs/fastq/113003636_S7_forward_paired.fq /media/glusterfs/genomics/Inputs/fastq/113003636_S7_reverse_paired.fq"
            ;;

        S8)
            INPUT_SEQUENCES="/media/glusterfs/genomics/Inputs/fastq/113003644_S8_forward_paired.fq /media/glusterfs/genomics/Inputs/fastq/113003644_S8_reverse_paired.fq"
            ;;

        S9)
            INPUT_SEQUENCES="/media/glusterfs/genomics/Inputs/fastq/113003652_S9_forward_paired.fq /media/glusterfs/genomics/Inputs/fastq/113003652_S9_reverse_paired.fq"
            ;;

        S10)
            INPUT_SEQUENCES="/media/glusterfs/genomics/Inputs/fastq/113003660_S10_forward_paired.fq /media/glusterfs/genomics/Inputs/fastq/113003660_S10_reverse_paired.fq"
            ;;
        S11)
            INPUT_SEQUENCES="/media/glusterfs/genomics/Inputs/fastq/113003678_S11_forward_paired.fq /media/glusterfs/genomics/Inputs/fastq/113003678_S11_reverse_paired.fq"
            ;;
        S12)
            INPUT_SEQUENCES="/media/glusterfs/genomics/Inputs/fastq/113003686_S12_forward_paired.fq /media/glusterfs/genomics/Inputs/fastq/113003686_S12_reverse_paired.fq"
            ;;
        S13)
            INPUT_SEQUENCES="/media/glusterfs/genomics/Inputs/fastq/113003694_S13_forward_paired.fq /media/glusterfs/genomics/Inputs/fastq/113003694_S13_reverse_paired.fq"
            ;;
        S14)
            INPUT_SEQUENCES="/media/glusterfs/genomics/Inputs/fastq/113003701_S14_forward_paired.fq /media/glusterfs/genomics/Inputs/fastq/113003701_S14_reverse_paired.fq"
            ;;
        S15)
            INPUT_SEQUENCES="/media/glusterfs/genomics/Inputs/fastq/113003719_S15_forward_paired.fq /media/glusterfs/genomics/Inputs/fastq/113003719_S15_reverse_paired.fq"
            ;;
        S16)
            INPUT_SEQUENCES="/media/glusterfs/genomics/Inputs/fastq/113003727_S16_forward_paired.fq /media/glusterfs/genomics/Inputs/fastq/113003727_S16_reverse_paired.fq"
            ;;

        KAL)
            INPUT_SEQUENCES="/media/glusterfs/genomics/Inputs/Kal_D/raw_lane2_R1.fq /media/glusterfs/genomics/Inputs/Kal_D/raw_lane2_R2.fq"
            ;;

        *)
            echo $"Usage: $0 {start|stop|restart|condrestart|status}"
            exit 1

esac


#Set dataset and output dir
DATABASE_NAME="/home/josemanuel.abuin/Genomics/Databases/${DATASET}_${PROCESSES_NUMBER}/DB_${DATASET}.db"
OUTPUT_DIR="/home/josemanuel.abuin/Genomics/Outputs/Output_${DATASET}_${INPUT_DATASET}_${PROCESSES_NUMBER}P_${THREADS_NUMBER}T.txt"

echo "Dataset: ${DATASET}"
echo "Number of processes: $PROCESSES_NUMBER"
echo "Number of threads: $THREADS_NUMBER"
echo "Nodefile: $NODES_FILE"
echo "Database name: $DATABASE_NAME"
echo "Output: $OUTPUT_DIR"
echo "Input sequences: $INPUT_SEQUENCES"

mpiexec -mca btl_tcp_if_include em1 -np ${PROCESSES_NUMBER} -hostfile ${NODES_FILE} ./metacache_mpi query $DATABASE_NAME ${INPUT_SEQUENCES} -lowest species -threads ${THREADS_NUMBER} -abundance-per species -pairfiles -out ${OUTPUT_DIR} -maxcand 4 -hitmin 4 -hitdiff 80 -query-limit ${QUERY_LIMIT}
