import matplotlib.pyplot as plt
import numpy as np

# Nombre del archivo de DOS magnética
filename = "TDOS.dat"  # cambia el nombre si es necesario

# Leer datos
data = np.loadtxt(filename)
energia = data[:, 0]
dos_up = data[:, 1]
dos_down = data[:, 2]  # signo negativo para visualizar hacia abajo

# Graficar
plt.figure(figsize=(8,6))
plt.plot(energia, dos_up, label='Spin ↑', color='red', linewidth=2)
plt.plot(energia, dos_down, label='Spin ↓', color='blue', linewidth=2)

# Línea vertical en E=0
plt.axvline(0.0, color='gray', linestyle='--', linewidth=1)

# Ejes y etiquetas
plt.title("Fe (FM)", fontsize=18)
plt.xlabel("Energía (eV)", fontsize=18)
plt.ylabel("DOS (Estados/eV)", fontsize=18)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()

emin = -8                       # Límite inferior del eje X (energía)
emax = 5                        # Límite superior del eje X (energía)
set_xlim = True

if set_xlim:
    plt.xlim(emin, emax)

# Guardar y mostrar
plt.savefig("dos-Fe-fm.png", dpi=300)
plt.show()
