#!/bin/bash
#SBATCH -n 8
#SBATCH --ntasks-per-node=4
#SBATCH -p thinnodes
#SBATCH -t 01:10:00
srun ./BuildAFS20_8_FT2.sh