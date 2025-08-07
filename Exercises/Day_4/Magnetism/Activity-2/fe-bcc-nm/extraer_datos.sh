#!/bin/bash

# Archivo de salida
output_file="data-Fe-nm.txt"

# Escribir encabezado en archivo y pantalla
echo -e "Carpeta\t\tVol_final(Å³)\tEnergía (eV)\tMagnetización (μB)" | tee "$output_file"

# Recorrer carpetas scale_*
for dir in scale_*; do
    outcar="$dir/OUTCAR"
    oszicar="$dir/OSZICAR"

    # Volumen final desde OUTCAR
    if [ -f "$outcar" ]; then
        vol_fin=$(grep "volume of cell" "$outcar" | tail -1 | awk '{print $5}')
    else
        vol_fin="---"
    fi

    # Energía y magnetización desde OSZICAR
    if [ -f "$oszicar" ]; then
        ultima_linea=$(tail -1 "$oszicar")
        energia=$(echo "$ultima_linea" | awk '{print $3}')
        energia=$(printf "%.8f" "$energia")
        mag=$(echo "$ultima_linea" | grep -o "mag= *[-0-9.]*" | awk -F= '{print $2}' | awk '{printf "%.3f", $1}')
    else
        energia="---"
        mag="---"
    fi

    # Imprimir en pantalla y guardar en archivo
    echo -e "$dir\t$vol_fin\t\t$energia\t$mag" | tee -a "$output_file"
done

