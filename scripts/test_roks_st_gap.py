"""Regression checks for the portable frozen-study tools."""
import tempfile
from pathlib import Path
import unittest
import roks_st_gap_common as c
from parse_roks_st_gap import donor_energy

class PortableFrozenTests(unittest.TestCase):
    def test_basis_label_is_verified(self):
        c.check_basis('BASIS_SET aug-cc-pVTZ\nBASIS_SET RI_AUX RI-5Z\n')
        with self.assertRaisesRegex(ValueError, 'aug-cc-pVTZ'):
            c.check_basis('BASIS_SET aug-cc-pVQZ\n')

    def test_hybrid_component_sum(self):
        text='''Total energy: -10.0
Exchange-correlation energy: -1.0
GAPW| Exc from hard and soft atomic rho1: -2.0
Hartree energy: 4.0
GAPW| local Eh = 1 center integrals: -3.0
Core Hamiltonian energy: -4.0
Overlap energy of the core charge distribution: 0.0
Self energy of the core charge distribution: -3.0
Hartree-Fock Exchange energy: -1.0
'''
        self.assertEqual(c.energy_components(text,-10)['exact_exchange'],-1)
        with self.assertRaisesRegex(ValueError,'Incomplete'):
            c.energy_components(text.replace('Hartree-Fock Exchange energy: -1.0',''),-10)

    def test_incomplete_donor_rejected(self):
        with tempfile.TemporaryDirectory() as root:
            path=Path(root)/'donor.out'
            path.write_text('SCF run converged in 10 steps\nENERGY| Total FORCE_EVAL ( QS ) energy [hartree] -10')
            with self.assertRaisesRegex(ValueError,'incomplete'):
                donor_energy(path)

    def test_rectangular_metric_rejected(self):
        with self.assertRaisesRegex(ValueError,'no complete square'):
            c.recover_overlap({'spins':[{'C':c.np.ones((3,2))}]})

    def test_density_dump_requires_both_spins(self):
        with tempfile.TemporaryDirectory() as root:
            path=Path(root)/'density.Log'
            path.write_text('DENSITY MATRIX FOR ALPHA SPIN\n 1\n 1 H 1 s 1.0\n')
            with self.assertRaisesRegex(ValueError,'Missing spin'):
                c.density_matrices(path,1)

if __name__=='__main__':
    unittest.main()
