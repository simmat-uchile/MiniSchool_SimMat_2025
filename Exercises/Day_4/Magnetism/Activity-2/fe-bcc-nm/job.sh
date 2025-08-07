#!/bin/bash                                                                     
#---------------Script SBATCH - NLHPC ----------------                          
#SBATCH -J foo                                                                  
#SBATCH -p general                                                              
#SBATCH -n 1                                                                   
#SBATCH --ntasks-per-node=1                                                    
#SBATCH -c 1                                                                    
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-user=                     
#SBATCH --mail-type=ALL                                                         
#SBATCH -o foo_%j.out                                                           
#SBATCH -e foo_%j.err                                                           

#-----------------Toolchain---------------------------                          
ml purge
ml intel/2022.00
# ----------------Modulos----------------------------                           
ml  VASP-VTST/6.2.1
# ----------------Comando--------------------------                             
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE

srun vasp_std > out2.log
