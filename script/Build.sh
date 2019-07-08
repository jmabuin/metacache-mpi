#!/usr/bin/env bash

mpirun -np 4 ./metacache_mpi build myRefseq /data/Software/Citius/metacache-mpi/data/datasets/test/ -taxonomy /data/Software/Citius/metacache-mpi/data/taxonomy/ -remove-overpopulated-features