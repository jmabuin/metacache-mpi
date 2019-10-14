#!/bin/bash
#SBATCH -n 32
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH -p thinnodes
#SBATCH -t 00:30:00
module load openmpi-runtime/2.1.6
srun ./QueryGeneric_FT.sh AFS31 KAL 32 4 3268000