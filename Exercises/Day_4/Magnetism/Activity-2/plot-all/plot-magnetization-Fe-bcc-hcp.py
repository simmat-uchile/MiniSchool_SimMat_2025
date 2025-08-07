import matplotlib.pyplot as plt

# Función para leer volumen y magnetización por átomo
def leer_magnetizacion(nombre_archivo):
    volumenes = []
    magnetizaciones = []
    with open(nombre_archivo, "r") as f:
        next(f)  # saltar encabezado
        for line in f:
            if line.strip() == "":
                continue
            parts = line.split()
            volumen = float(parts[1])
            mag_total = float(parts[3])
            mag_por_atomo = mag_total / 2.0
            volumenes.append(volumen)
            magnetizaciones.append(mag_por_atomo)
    return volumenes, magnetizaciones

# Leer datos
vol1, mag1 = leer_magnetizacion("data-Fe-fm.txt")
vol2, mag2 = leer_magnetizacion("data-Fe-hcp.txt")

# Graficar
plt.figure(figsize=(8,6))

plt.plot(vol1, mag1, marker='o', linestyle='-', color='blue', linewidth=2, label='Fe-FM bcc')
plt.plot(vol2, mag2, marker='s', linestyle='--', color='red', linewidth=2, label='Fe-NM hcp')

# Estética
plt.title("Magnetización por átomo vs Volumen", fontsize=16)
plt.xlabel("Volumen (Å³)", fontsize=14)
plt.ylabel("Magnetización (μB / átomo)", fontsize=14)
plt.grid(True)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("magnetizacion_vs_volumen.png", dpi=300)
plt.show()
