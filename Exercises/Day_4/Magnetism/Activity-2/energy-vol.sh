#!/bin/bash

# Lista de parámetros de escala
parametros=(1.05 1.06 1.07)

# Archivos que se copian sin modificación
archivos_base=(INCAR POTCAR KPOINTS)

# Verificación de existencia
for archivo in "${archivos_base[@]}"; do
    if [ ! -f "$archivo" ]; then
        echo "❌ ERROR: Archivo $archivo no encontrado en el directorio actual."
        exit 1
    fi
done

# Loop sobre parámetros
for p in "${parametros[@]}"; do
    carpeta="scale_${p}"
    mkdir -p "$carpeta"

    # Crear POSCAR con estructura fija y escala $p
    cat > "$carpeta/POSCAR" <<EOF
Fe
$p
  2.8312852282750209   0.0000000000000000    0.0000000000000002
  0.0000000000000002   2.8312852282750209    0.0000000000000002
  0.0000000000000000   0.0000000000000000    2.8312852282750209
Fe
2
direct
   0.0000000000000000    0.0000000000000000    0.0000000000000000 Fe0+
   0.5000000000000000    0.5000000000000000    0.5000000000000000 Fe0+
EOF

    # Copiar archivos comunes
    for archivo in "${archivos_base[@]}"; do
        cp "$archivo" "$carpeta/"
    done

    # Crear job.sh personalizado
    cat > "$carpeta/job.sh" <<EOF
#!/bin/bash
#---------------Script SBATCH - NLHPC ----------------
#SBATCH -J ev_${p}
#SBATCH -p general
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o ev_${p}_%j.out
#SBATCH -e ev_${p}_%j.err

#-----------------Toolchain---------------------------
ml purge
ml intel/2022.00

# ----------------Modulos----------------------------
ml VASP-VTST/6.2.1

# ----------------Comando----------------------------
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE

srun vasp_std > out2.log
EOF

    # Enviar el trabajo al scheduler
    (cd "$carpeta" && sbatch job.sh)

    echo "📁 Carpeta $carpeta creada y job ev_${p} enviado."
done

