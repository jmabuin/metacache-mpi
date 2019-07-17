#!/usr/bin/env bash

BIN_PATH=/mnt/lustre/scratch/home/usc/ec/jam/Genomica/metacache-mpi/build/metacache_mpi
DATABASE_NAME=/mnt/lustre/scratch/home/usc/ec/jam/Genomica/Databases/AFS20_8/DB_ref.db
INPUT_SEQUENCES="/mnt/lustre/scratch/home/usc/ec/jam/Genomica/Inputs/Kal_D/raw_lane2_R1.fq /mnt/lustre/scratch/home/usc/ec/jam/Genomica/Inputs/Kal_D/raw_lane2_R2.fq"
NUM_THREADS=4
OUTPUT_FILE=Saida_AFS20_8_4T_KAL

$BIN_PATH query $DATABASE_NAME $INPUT_SEQUENCES -lowest species -threads $NUM_THREADS -abundance-per species -pairfiles -out $OUTPUT_FILE -maxcand 4 -hitmin 4 -hitdiff 80 -query-limit 1000000
