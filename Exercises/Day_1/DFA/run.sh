#!/bin/bash
#SBATCH --job-name=G09
#SBATCH --partition=labs
#SBATCH -c 2 #numero de CPUs
#SBATCH --output=example.out
#SBATCH --error=example.err
#SBATCH --mem=1500
#SRUN --export=ALL

ml intel/2019b 
ml g09/D01


g09 Be.gjf > Be.log
g09 Be_ccsd.gjf > Be_ccsd.log


