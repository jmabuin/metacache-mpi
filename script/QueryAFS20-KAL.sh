#!/usr/bin/env bash

mpirun -np 4 ./metacache_mpi query DB_AFS20.db /home/remoto/josemanuel.abuin/Datasets/afs_kalibrator/Kal_D/raw_lane2_R1.fq /home/remoto/josemanuel.abuin/Datasets/afs_kalibrator/Kal_D/raw_lane2_R2.fq -lowest species -threads 4 -abundance-per species -pairfiles -out SaidaKAL -maxcand 4 -hitmin 4 -hitdiff 80 -query-limit 500000
