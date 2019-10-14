#!/bin/bash
#SBATCH -n 32
#SBATCH --ntasks-per-node=1
#SBATCH -c 16
#SBATCH -p thinnodes
#SBATCH -t 00:30:00
module load openmpi-runtime/2.1.6
srun ./QueryGeneric_FT.sh AFS20RS90 KAL 32 16 817000