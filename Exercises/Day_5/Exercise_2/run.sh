#!/bin/bash
#SBATCH -J QE
#SBATCH -p labs
#SBATCH --reservation=abinitio
#SBATCH -n 10
#SBATCH --ntasks-per-node=10  
#SBATCH --mem-per-cpu=1000 
#SBATCH -o %j_%x.out
#SBATCH -e %j_%x.err

ml purge 
ml intel-compilers/2022.0.1  
ml impi/2021.5.0
ml  QuantumESPRESSO/7.2

srun pw.x < opt.inp > opt.out
