#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --partition=general
#SBATCH --job-name=TEST
#SBATCH --output=TEST
#SBATCH --time=01:00:00

cd /home/uqlyengo/h2WGS/splitGRM/

srun ./splitGRM --grm ../16grms_pgen/WGS_GRM --ngrp 3 --out WGS_GRM


