import matplotlib.pyplot as plt

# Valor mínimo de referencia para la energía
E0 = -16.47455800

# Función para leer volumen y energía relativa desde archivo
def leer_datos(nombre_archivo):
    volumenes = []
    energias_rel = []
    with open(nombre_archivo, "r") as f:
        next(f)  # Saltar encabezado
        for line in f:
            if line.strip() == "":
                continue
            parts = line.split()
            volumen = float(parts[1])
            energia = float(parts[2])
            volumenes.append(volumen)
            energias_rel.append(energia - E0)
    return volumenes, energias_rel

# Leer los tres archivos
vol1, ener1 = leer_datos("data-Fe-fm.txt")
vol2, ener2 = leer_datos("data-Fe-nm.txt")
vol3, ener3 = leer_datos("data-Fe-hcp.txt")

# Graficar
plt.figure(figsize=(8,6))

plt.plot(vol1, ener1, 'o-', color='blue', linewidth=2, label='Estructura 1 (Fe FM.txt)')
plt.plot(vol2, ener2, 's--', color='red', linewidth=2, label='Estructura 2 (Fe NM.txt)')
plt.plot(vol3, ener3, 'd-.', color='green', linewidth=2, label='Estructura 3 (Fe-hcp.txt)')

# Estética
plt.title("Comparación: Energía relativa vs Volumen", fontsize=16)
plt.xlabel("Volumen (Å³)", fontsize=14)
plt.ylabel("Energía - E₀ (eV)", fontsize=14)
plt.grid(True)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("energia_vs_volumen_compare.png", dpi=300)
plt.show()
