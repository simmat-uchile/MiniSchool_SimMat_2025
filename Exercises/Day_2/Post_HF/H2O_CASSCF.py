#!/usr/bin/env python

from pyscf import gto, scf, mcscf

mol = gto.M(
    atom   = 'O  0.   0.   0.;  H  0.   0.   0.96;  H  0.   0.76  0.58',
    basis  = 'cc-pVDZ',
    charge = 0,
    spin   = 0,       # closed shell
)
mf = scf.RHF(mol).run()

ncas    = 6       # e.g. 6 orbitals in active space
nelecas = 6       # e.g. 6 electrons in those orbitals

mc = mcscf.CASSCF(mf, ncas, nelecas)
results = mc.kernel()
e_casscf = results[0]

print(f'CASSCF energy: {e_casscf:.8f}  a.u.')
