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
# Cálculo manual de energía MP2 
# =========================================
import numpy as np
from pyscf import ao2mo

# Obtener integrales ERI en base MO, notación químico (pq|rs)
nmo = mf.mo_coeff.shape[1]
eri_1d = ao2mo.kernel(mol, mf.mo_coeff)
eri_mo = ao2mo.restore(1, eri_1d, nmo)  # tensor 4D (p,q,r,s)

eps = mf.mo_energy
nocc = mol.nelectron // 2


MP2_corr_E = ? #Completar script
