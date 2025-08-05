#!/bin/bash
#SBATCH --job-name=G09
#SBATCH -p general
#SBATCH -c 1 #numero de CPUs
#SBATCH --output=example.out
#SBATCH --error=example.err
#SBATCH --mem=1000
#SRUN --export=ALL

ml intel/2019b
ml g09/D01

formchk benceno_GS.chk 
formchk benceno_state1.chk 
formchk benceno_FC_NO.chk 

