#!/usr/bin/env bash

BIN_PATH=/mnt/lustre/scratch/home/usc/ec/jam/Genomica/metacache-mpi/build/metacache_mpi
DATABASE_NAME=/mnt/lustre/scratch/home/usc/ec/jam/Genomica/Databases/AFS20_8/DB_ref.db
INPUT_SEQUENCES="/mnt/lustre/scratch/home/usc/ec/jam/Genomica/Inputs/fastq/113003602_S4_forward_paired.fq /mnt/lustre/scratch/home/usc/ec/jam/Genomica/Inputs/fastq/113003602_S4_reverse_paired.fq"
NUM_THREADS=4
OUTPUT_FILE=Saida_AFS20_8_2T_S4

$BIN_PATH query $DATABASE_NAME $INPUT_SEQUENCES -lowest species -threads $NUM_THREADS -abundance-per species -pairfiles -out $OUTPUT_FILE -maxcand 4 -hitmin 4 -hitdiff 80 -query-limit 500000
