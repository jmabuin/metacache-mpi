#!/usr/bin/env bash

mpirun -np 4 ./metacache_mpi build DB_AFS20.db /home/remoto/josemanuel.abuin/Datasets/AFS20-Robin/ -taxonomy /home/remoto/josemanuel.abuin/Datasets/NCBI/Refseq-Release90-Taxonomy/ -remove-overpopulated-features