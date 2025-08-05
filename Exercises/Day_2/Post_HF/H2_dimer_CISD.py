#!/usr/bin/env python


# =========================================
# Cálculo Hartree-Fock (HF)
# =========================================
from pyscf import gto, scf
print("\nSystem: H2\n")
print("\n================== Hartree-Fock ================\n")

mol = gto.M(
    atom = """
            H   0.0 0.0 0.0
            H   0.0 0.0 0.74
           """,
    basis = 'STO-3g',
    charge = 0,
    spin = 0
)

mf = scf.RHF(mol)
energia_hf = mf.kernel()

print("\tEnergía HF:", energia_hf)

# =========================================
# Cálculo CISD
# =========================================
from pyscf import ci
print("\n====================== CISD ====================\n")

cisd_calc = ci.CISD(mf)
e_corr_cisd, _ = cisd_calc.kernel()
e_tot_cisd = energia_hf + e_corr_cisd
print("\tEnergía CISD:", e_tot_cisd)
print("\tCorrelación FCI:", e_corr_cisd)


# =========================================
# Cálculo FCI (exacto en base finita)
# =========================================
from pyscf import fci

fci_calc = fci.FCI(mol, mf.mo_coeff)
energia_fci, _ = fci_calc.kernel()

print("\n======== Full Configuration Interactions ========\n")
print("\tEnergía FCI:", energia_fci)
print("\tCorrelación FCI:", energia_fci - energia_hf)


#--------------------------------------------------------
print("\nSystem: H2-H2\n")
print("\n================== Hartree-Fock ================\n")

mol = gto.M(
    atom = """
            H   0.0 0.0 0.0
            H   0.0 0.0 0.74
            H   0.0 0.0 10.0
            H   0.0 0.0 10.74
           """,
    basis = 'STO-3G',
    charge = 0,
    spin = 0
)

mf = scf.RHF(mol)
energia_hf = mf.kernel()

print("\tEnergía HF:", energia_hf)

# =========================================
# Cálculo CISD
# =========================================
from pyscf import ci
print("\n====================== CISD ====================\n")

cisd_calc = ci.CISD(mf)
e_corr_cisd, _ = cisd_calc.kernel()
e_tot_cisd = energia_hf + e_corr_cisd
print("\tEnergía CISD:", e_tot_cisd)
print("\tCorrelación FCI:", e_corr_cisd)


# =========================================
# Cálculo FCI (exacto en base finita)
# =========================================
from pyscf import fci
print("\n======== Full Configuration Interactions ========\n")

fci_calc = fci.FCI(mol, mf.mo_coeff)
energia_fci, _ = fci_calc.kernel()

print("\tEnergía FCI:", energia_fci)
print("\tCorrelación FCI:", energia_fci - energia_hf)

