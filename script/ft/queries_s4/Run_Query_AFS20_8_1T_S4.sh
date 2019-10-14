#!/bin/bash
#SBATCH -n 8
#SBATCH --ntasks-per-node=1
#SBATCH -p thinnodes
#SBATCH -t 00:30:00
module load openmpi-runtime/2.1.6
srun ./QueryGeneric_FT.sh AFS20 KAL 8 1 13072000