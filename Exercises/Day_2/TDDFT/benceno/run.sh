#!/bin/bash
#SBATCH --job-name=G09
#SBATCH -p labs
#SBATCH --reservation=abinitio
#SBATCH -c 2 #numero de CPUs
#SBATCH --output=example.out
#SBATCH --error=example.err
#SBATCH --mem=8000
#SRUN --export=ALL

ml intel/2019b
ml g09/D01

g09 benceno.com > benceno.log

