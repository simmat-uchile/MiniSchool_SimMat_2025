# -*- coding: utf-8 -*-


run  pip install --upgrade pyprocar
import pyprocar
import numpy as np
import matplotlib.pyplot as plt

# 2. Llamar a la función para graficar las bandas
#    pyprocar buscará automáticamente los archivos necesarios en el directorio.
print("Generando el gráfico de la estructura de bandas...")

pyprocar.bandsplot(
    code='qe',                 # Especifica que el cálculo se hizo con Quantum Espresso.
    dirname='./',              # Directorio donde están los archivos de salida (el actual).
    mode='plain',              # Modo de gráfico 'simple', sin proyecciones de color.
    elimit=[-10,18],
    fermi=17.5842,             # Nivel de Fermi en eV para dibujarlo como referencia.
    savefig='./bands1.png',    # Guarda la figura resultante con este nombre y formato.
    show=False                 # No muestra la figura en una ventana emergente (útil para scripts).
)

print("\n¡Gráfico guardado exitosamente como 'bands1.png'!")


