#!/bin/bash
#SBATCH -n 64
#SBATCH --ntasks-per-node=2
#SBATCH -p thinnodes
#SBATCH -t 01:20:00
module load openmpi-runtime/2.1.6
srun ./BuildGeneric_FT.sh AFS31 64