#!/usr/bin/env python

# =========================================
# Cálculo Hartree-Fock (HF)
# =========================================
from pyscf import gto, scf

# Definimos la molécula de H2
mol = gto.M(
    atom = 'H 0 0 0; H 0 0 0.74',
    basis = 'sto-3g'
)

# Ejecutamos SCF (Hartree-Fock restringido)
mf = scf.RHF(mol)
energia_hf = mf.kernel()
print("Energía HF:", energia_hf)

# =========================================
# Cálculo MP2
# =========================================
from pyscf import mp

mp2_calc = mp.MP2(mf)
E_corr_mp2, t2 = mp2_calc.kernel()
print("Energía MP2:", E_corr_mp2 + energia_hf)
print("Correlación MP2:", E_corr_mp2)

# =========================================
# Cálculo CCSD
# =========================================
from pyscf import cc

ccsd_calc = cc.CCSD(mf)
E_corr_ccsd, t1, t2 = ccsd_calc.kernel()
print("Energía CCSD:", E_corr_ccsd + energia_hf)
print("Correlación CCSD:", E_corr_ccsd)

# =========================================
# Correcciones para excitaciones triples (T)
# =========================================

E_correction_T = ccsd_calc.ccsd_t()

# Total energy = HF + CCSD correlation + (T) correction
energia_total_ccsdt = energia_hf + E_corr_ccsd + E_correction_T

print("Corrección (T):", E_correction_T)
print("Energía total CCSD(T):", energia_total_ccsdt)


# =========================================
# Cálculo FCI (exacto en base finita)
# =========================================
from pyscf import fci

fci_calc = fci.FCI(mol, mf.mo_coeff)
energia_fci, _ = fci_calc.kernel()
print("Energía FCI:", energia_fci)
print("Correlación FCI:", energia_fci - energia_hf)

# =========================================
# Orbitales naturales (opcional)
# =========================================
import numpy as np

# Matriz de densidad de 1 partícula a partir de MP2
rdm1 = mp2_calc.make_rdm1()
valores_propios = np.linalg.eigvalsh(rdm1)
print("Ocupaciones de orbitales naturales (MP2):")
print(valores_propios)

# =========================================
# Notas adicionales:
# - Se puede repetir con diferentes bases (cc-pVDZ, cc-pVTZ)
# - Para comparar energías, conviene usar la misma molécula y geometría
# - Si hay tiempo, se puede hacer un escaneo del enlace H-H

