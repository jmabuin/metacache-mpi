#!/bin/bash
#SBATCH -n 8
#SBATCH --ntasks-per-node=4
#SBATCH -c 4
#SBATCH -p thinnodes
#SBATCH -t 01:10:00
srun ./QueryAFS20_8_KAL_FT2_4T.sh