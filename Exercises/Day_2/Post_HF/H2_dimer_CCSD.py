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
# Cálculo CCSD
# =========================================
from pyscf import cc
print("\n=================== CCSD ======================\n")

ccsd_calc = cc.CCSD(mf)
E_corr_ccsd, t1, t2 = ccsd_calc.kernel()
print("Energía CCSD:", E_corr_ccsd + energia_hf)
print("Correlación CCSD:", E_corr_ccsd)

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
# Cálculo CCSD
# =========================================
from pyscf import cc
print("\n=================== CCSD ======================\n")

ccsd_calc = cc.CCSD(mf)
E_corr_ccsd, t1, t2 = ccsd_calc.kernel()
print("Energía CCSD:", E_corr_ccsd + energia_hf)
print("Correlación CCSD:", E_corr_ccsd)
