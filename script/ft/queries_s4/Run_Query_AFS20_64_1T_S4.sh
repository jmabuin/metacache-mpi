#!/bin/bash
#SBATCH -n 64
#SBATCH --ntasks-per-node=2
#SBATCH -p thinnodes
#SBATCH -t 00:30:00
module load openmpi-runtime/2.1.6
srun ./QueryGeneric_FT.sh AFS20 KAL 64 1 13072000