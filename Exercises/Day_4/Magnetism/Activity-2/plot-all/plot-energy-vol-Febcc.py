import matplotlib.pyplot as plt

# Parámetro de energía mínima (eV)
E0 = -16.47455800

# Listas para almacenar datos
volumenes = []
energias_relativas = []

# Leer data.txt ignorando la cabecera
with open("data-Fe-fm.txt", "r") as f:
    next(f)  # saltar encabezado
    for line in f:
        if line.strip() == "":
            continue
        parts = line.split()
        volumen = float(parts[1])
        energia = float(parts[2])
        volumenes.append(volumen)
        energias_relativas.append(energia - E0)

# Graficar
plt.figure(figsize=(8,6))
plt.plot(volumenes, energias_relativas, marker='o', linestyle='-', color='blue', linewidth=2)

# Estilo
plt.title("Curva Energía vs Volumen", fontsize=16)
plt.xlabel("Volumen (Å³)", fontsize=14)
plt.ylabel("Energía - E₀ (eV)", fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.savefig("curva_energia_volumen.png", dpi=300)
plt.show()

