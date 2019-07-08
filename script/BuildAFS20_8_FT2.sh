#!/usr/bin/env bash

BIN_PATH=/mnt/lustre/scratch/home/usc/ec/jam/Genomica/metacache-mpi/build/metacache_mpi
DATABASE_NAME=/mnt/lustre/scratch/home/usc/ec/jam/Genomica/Databases/AFS20_8/DB_ref.db
INPUT_SEQUENCES=/mnt/lustre/scratch/home/usc/ec/jam/Genomica/AFS20/
TAXONOMY=/mnt/lustre/scratch/home/usc/ec/jam/Genomica/Taxonomy/

$BIN_PATH build $DATABASE_NAME $INPUT_SEQUENCES -taxonomy $TAXONOMY -remove-overpopulated-features
