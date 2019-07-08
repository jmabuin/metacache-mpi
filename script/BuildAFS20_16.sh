#!/usr/bin/env bash

mpirun -np 16 -hostfile ~/nodefile /media/glusterfs/genomics/metacache-mpi/build/metacache_mpi build /media/glusterfs/genomics/Databases/AFS20_16/DB_AFS20.db /media/glusterfs/genomics/AFS20/ -taxonomy /media/glusterfs/genomics/Taxonomy/ -remove-overpopulated-features
