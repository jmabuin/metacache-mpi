#!/bin/bash
#SBATCH -n 16
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH -p thinnodes
#SBATCH -t 00:30:00
module load openmpi-runtime/2.1.6
srun ./QueryGeneric_FT.sh AFS20 KAL 16 4 3268000