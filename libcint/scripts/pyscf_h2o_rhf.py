from pyscf import gto, lib, scf
import numpy as np
import scipy

mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5", basis="def2-TZVP").build()

# start of rhf total energy algorithms
coords = mol.atom_coords()
charges = mol.atom_charges()
dist = scipy.spatial.distance.cdist(coords, coords)
np.fill_diagonal(dist, np.inf)
eng_nuc = 0.5 * (charges * charges[:, None] / dist).sum()
print(f"Total energy nuc: {eng_nuc}")

hcore = mol.intor("int1e_kin") + mol.intor("int1e_nuc")
ovlp = mol.intor("int1e_ovlp")
int2e = mol.intor("int2e")
nocc = 5

dm = np.zeros_like(ovlp)
for _ in range(40):
    # hardcoded SCF iterations
    fock = hcore + ((1 * int2e - 0.5 * int2e.swapaxes(1, 2)) * dm).sum(axis=(-1, -2))
    e, c = scipy.linalg.eigh(fock, ovlp)
    dm = 2 * c[:, :nocc] @ c[:, :nocc].T

eng_scratch = 1 * hcore + ((0.5 * int2e - 0.25 * int2e.swapaxes(1, 2)) * dm).sum(axis=(-1, -2))
eng_elec = (dm * eng_scratch).sum()
print(f"Total elec energy: {eng_elec}")
print(f"Total RHF energy: {eng_elec + eng_nuc}")
# end of rhf total energy algorithms

mf = scf.RHF(mol).run()
print(mf.energy_elec()[0], mf.energy_nuc(), mf.energy_tot())
