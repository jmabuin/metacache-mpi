#!/bin/bash
#SBATCH -n 8
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH -p thinnodes
#SBATCH -t 00:30:00
module load openmpi-runtime/2.1.6
srun ./QueryGeneric_FT.sh AFS20RS90 KAL 8 4 3268000