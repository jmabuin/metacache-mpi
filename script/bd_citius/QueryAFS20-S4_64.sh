#!/usr/bin/env bash

mpiexec -mca btl_tcp_if_include em1 -np 64 -hostfile ~/nodefile ./metacache_mpi query /media/glusterfs/genomics/Databases/AFS20_64/DB_AFS20.db /media/glusterfs/genomics/Inputs/fastq/113003602_S4_forward_paired.fq /media/glusterfs/genomics/Inputs/fastq/113003602_S4_reverse_paired.fq -lowest species -threads 4 -abundance-per species -pairfiles -out /media/glusterfs/genomics/Outputs/SaidaS4_64_4T -maxcand 4 -hitmin 4 -hitdiff 80 -query-limit 500000
