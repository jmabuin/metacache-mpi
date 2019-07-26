#!/usr/bin/env bash

mpiexec -mca btl_tcp_if_include em1 -np 16 -hostfile ~/nodefile ./metacache_mpi query /media/glusterfs/genomics/Databases/AFS20_16/DB_AFS20.db /media/glusterfs/genomics/Inputs/Kal_D/raw_lane2_R1.fq /media/glusterfs/genomics/Inputs/Kal_D/raw_lane2_R2.fq -lowest species -threads 4 -abundance-per species -pairfiles -out /media/glusterfs/genomics/Outputs/SaidaKAL_D_16_4T -maxcand 4 -hitmin 4 -hitdiff 80 -query-limit 100000
