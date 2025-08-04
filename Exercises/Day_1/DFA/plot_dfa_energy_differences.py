import re
import matplotlib.pyplot as plt

# === INPUT LOG FILE ===
logfile = "Be.log"

# === REFERENCE ENERGY (e.g., CCSD) ===
reference_energy = -14.6174605  # Hartree

# === PARSE LOG FILE ===
methods = []
energies = []

with open(logfile, 'r') as f:
    for line in f:
        if "SCF Done:" in line:
            match = re.search(r'SCF Done:\s+E\(([^)]+)\)\s+=\s+(-?\d+\.\d+)', line)
            if match:
                method = match.group(1)
                energy = float(match.group(2))
                methods.append(method)
                energies.append(energy)

# === COMPUTE ENERGY DIFFERENCES ===
delta_E_mHa = [(E - reference_energy) * 1000 for E in energies]

# === PLOT ===
plt.figure()
plt.bar(methods, delta_E_mHa, color="royalblue")
plt.axhline(0, color='black', linestyle='--')
plt.ylabel("ΔE (mHartree)")
plt.title("SCF Energy Deviations from CCSD")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig("scf_energy_deviations_dfa.png", dpi=300)

