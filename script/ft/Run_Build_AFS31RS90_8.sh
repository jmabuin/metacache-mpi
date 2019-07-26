#!/bin/bash
#SBATCH -n 8
#SBATCH --ntasks-per-node=1
#SBATCH -p thinnodes
#SBATCH -t 03:10:00
srun ./BuildGeneric_FT.sh AFS31RS90 8