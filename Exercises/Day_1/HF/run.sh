#!/bin/bash
#SBATCH --job-name=G09
#SBATCH --partition=general
#SBATCH -c 2 #numero de CPUs
#SBATCH --output=example.out
#SBATCH --error=example.err
#SBATCH --mem=1500
#SRUN --export=ALL

ml intel/2019b 
ml g09/D01


g09 He.gjf > He.log #linea descomentada

formchk He_sto3g.chk He_sto3g.fchk
#formchk He_sto6g.chk He_sto6g.fchk
#formchk He_631g.chk He_631g.fchk
#formchk He_6311g++g3df3pd.chk He_6311g++g3df3pd.fchk
formchk He_augccpv6z.chk He_augccpv6z.fchk

# Crear archivo con los parámetros de la caja
cat > cubegen_box.txt << EOF
-1  -5.0  -5.0  -5.0
300   0.033 0.000  0.000
300   0.000 0.033  0.000
300   0.000 0.000  0.033
EOF

# Name
BASENAME="He_sto3g"

# Generar densidad electrónica
cubegen 1 FDensity=SCF "${BASENAME}.fchk" "${BASENAME}_density.cube" -1 h < cubegen_box.txt

# Crear archivo con los parámetros de la caja
cat > cubegen_box.txt << EOF
-1  -5.0  -5.0  -5.0
300   0.033 0.000  0.000
300   0.000 0.033  0.000
300   0.000 0.000  0.033
EOF

# Name
BASENAME="He_augccpv6z"

# Generar densidad electrónica
cubegen 1 FDensity=SCF "${BASENAME}.fchk" "${BASENAME}_density.cube" -1 h < cubegen_box.txt

cat > cubman.in << EOF
ToXYZ
He_sto3g_density.cube
yes
He_sto3g_density.xyz
yes
EOF

cubman < cubman.in > He_sto3g.out


cat > cubman.in << EOF
ToXYZ
He_augccpv6z_density.cube
yes
He_augccpv6z_density.xyz
yes
EOF

cubman < cubman.in > He_augccpv6z.out

