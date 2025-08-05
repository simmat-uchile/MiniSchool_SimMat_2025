# =========================================
# Cálculo Hartree-Fock (HF)
# =========================================
from pyscf import gto, scf

mol = gto.M(
    atom = 'H 0 0 0',
    basis = 'cc-PVDZ',
    charge = 0,
    spin = 1
)

mf = scf.UHF(mol)
energia_hf = mf.kernel()

print("\n================== Hartree-Fock ================\n")
print("\tEnergía HF:", energia_hf)


# =========================================
# Cálculo FCI (exacto en base finita)
# =========================================
from pyscf import fci

fci_calc = fci.FCI(mol, mf.mo_coeff)
energia_fci, _ = fci_calc.kernel()

print("\n======== Full Configuration Interactions ========\n")
print("\tEnergía FCI:", energia_fci)
print("\tCorrelación FCI:", energia_fci - energia_hf)

