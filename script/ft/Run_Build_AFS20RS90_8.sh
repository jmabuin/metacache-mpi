#!/bin/bash
#SBATCH -n 8
#SBATCH --ntasks-per-node=2
#SBATCH -p thinnodes
#SBATCH -t 03:10:00
module load openmpi-runtime/2.1.6
srun ./BuildGeneric_FT.sh AFS20RS90 8