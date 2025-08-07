import matplotlib.pyplot as plt
import numpy as np

# Nombre del archivo DOS
filename = "TDOS.dat"  # puedes cambiarlo si tu archivo tiene otro nombre

# Leer datos
data = np.loadtxt(filename)
energia = data[:, 0]
dos_total = data[:, 1]

# Graficar DOS
plt.figure(figsize=(7,6))
plt.plot(energia, dos_total, color='purple', linewidth=2)

# Línea vertical en E = 0 (nivel de Fermi)
plt.axvline(x=0.0, color='gray', linestyle='--', linewidth=1)

# Estética
plt.title("Densidad de Estados Total", fontsize=16)
plt.xlabel("Energía (eV)", fontsize=14)
plt.ylabel("DOS (Estados/eV)", fontsize=14)
plt.grid(True)
plt.tight_layout()

# Guardar y mostrar
plt.savefig("dos-Fe-nm.png", dpi=300)
plt.show()
