#!/usr/bin/env bash

mpirun -np 4 ./metacache_mpi query myRefseq.db /data/Software/Citius/metacache-mpi/data/datasets/test_query/113003694_S13_forward_paired.fq /data/Software/Citius/metacache-mpi/data/datasets/test_query/113003694_S13_reverse_paired.fq -lowest species -threads 2 -abundance-per species -pairfiles -out Saida