#!/bin/bash
#SBATCH -n 16
#SBATCH --ntasks-per-node=2
#SBATCH -p thinnodes
#SBATCH -t 02:10:00
module load openmpi-runtime/2.1.6
srun ./BuildGeneric_FT.sh AFS31RS90 16