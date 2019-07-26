#!/usr/bin/env bash

mpirun -np 8 ./metacache_mpi query DB_AFS20.db /home/remoto/josemanuel.abuin/Datasets/afs_kalibrator/fastq/113003602_S4_forward_paired.fq /home/remoto/josemanuel.abuin/Datasets/afs_kalibrator/fastq/113003602_S4_reverse_paired.fq -lowest species -threads 2 -abundance-per species -pairfiles -out SaidaS4_New -maxcand 4 -hitmin 4 -hitdiff 80
