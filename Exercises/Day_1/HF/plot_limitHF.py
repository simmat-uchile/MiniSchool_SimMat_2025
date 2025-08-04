import re
import matplotlib.pyplot as plt

# === Archivo log con todos los cálculos ===
log_file = "He.log"

# === Nombres de las bases ===
basis_names = [
    "STO-3G",
    "STO-6G",
    "6-31G",
    "6-311G",
    "6-311++G(3df,3pd)"
]

# Leer el archivo completo
with open(log_file, 'r') as f:
    text = f.read()

# === 1. Extraer energías HF ===
hf_matches = re.findall(r'\\HF=(-?\d+\.\d+)', text)
hf_energies = [float(e) for e in hf_matches]

# === 2. Extraer los HOMO (primer valor de alpha occ. eigenvalues) ===
pi_matches = re.findall(r'alpha\s+occ\. eigenvalues --\s+(-?\d+\.\d+)', text, re.IGNORECASE)
koopmans_ionization = [-float(e) for e in pi_matches]  # aplicar -1

# === Verificar consistencia ===
if len(hf_energies) != len(basis_names) or len(koopmans_ionization) != len(basis_names):
    print("[ERROR] La cantidad de datos extraídos no coincide con el número de bases.")
    print(f"HF encontrados: {len(hf_energies)}, HOMO encontrados: {len(koopmans_ionization)}, Bases: {len(basis_names)}")
    exit()

# === 3. Graficar HF vs basis ===
plt.figure(figsize=(9, 5))
plt.plot(basis_names, hf_energies, marker='o', linestyle='-')
plt.xlabel("Basis set")
plt.ylabel("Energía total HF (Hartree)")
plt.title("Energía Hartree-Fock vs basis set")
plt.tight_layout()
plt.savefig("hf_vs_basis.png", dpi=300)

# === 4. Graficar Koopmans vs basis con línea experimental ===
plt.figure(figsize=(9, 5))
plt.plot(basis_names, koopmans_ionization, marker='s', linestyle='-', label=r"$- \varepsilon_{\mathrm{HOMO}}$ (Koopmans)")
plt.axhline(0.903569, color='red', linestyle='--', label="Experiment (0.903569 a.u.)")
plt.xlabel("Basis set")
plt.ylabel("Ionization energy (Hartree)")
plt.title("Koopmans' Theorem")
plt.legend()
plt.tight_layout()
plt.savefig("koopmans_vs_basis.png", dpi=300)
