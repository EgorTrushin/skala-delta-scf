"""Run with: python -m unittest discover -s scripts -p 'test_*.py'."""
import tempfile
from pathlib import Path
import unittest

from parse_cp2k import HARTREE_TO_EV, collect, main, parse_output


def log(spin, energy):
    return (f'PROGRAM STARTED AT now\nDFT| Calculation {spin}\n'
            'SCF WAVEFUNCTION OPTIMIZATION\nSCF run converged in 5 steps\n'
            f'ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] {energy}\n'
            'PROGRAM ENDED AT now\n')


class ParserTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def case(self, method='uks', dataset='QUEST1'):
        p = self.root / 'aug-cc-pVTZ' / dataset / method / 'pbe' / 'Water'
        p.mkdir(parents=True)
        for stage, spin, energy in [('ground_state', 'RKS', -10),
                                     ('singlet', method.upper(), -9.8),
                                     ('triplet', method.upper(), -9.9)]:
            (p / (stage + '.out')).write_text(log(spin, energy))
        return p

    def test_uks_and_roks_have_different_singlet_formulas(self):
        self.case('uks'); self.case('roks')
        rows = {r['method']: r for r in collect(self.root)}
        self.assertAlmostEqual(rows['UKS']['singlet_eV'], 0.3 * HARTREE_TO_EV)
        self.assertAlmostEqual(rows['ROKS']['singlet_eV'], 0.2 * HARTREE_TO_EV)
        for r in rows.values():
            self.assertAlmostEqual(r['triplet_eV'], 0.1 * HARTREE_TO_EV)
            self.assertAlmostEqual(r['gap_eV'], r['singlet_eV'] - r['triplet_eV'])

    def test_quest6_uses_triplet_but_exports_singlet_only(self):
        self.case(dataset='QUEST6')
        row, = collect(self.root)
        self.assertAlmostEqual(row['singlet_eV'], 0.3 * HARTREE_TO_EV)
        self.assertEqual((row['triplet_eV'], row['gap_eV']), ('', ''))

    def test_rejects_appended_incomplete_run(self):
        p = self.case() / 'singlet.out'
        with p.open('a') as f:
            f.write('PROGRAM STARTED AT later\nDFT| Calculation UKS\n')
        with self.assertRaises(ValueError):
            parse_output(p, 'UKS')

    def test_rejects_wrong_spin_and_failed_or_truncated_outputs(self):
        p = self.case() / 'singlet.out'
        for content in [log('ROKS', -9.8),
                        log('UKS', -9.8).replace('SCF run converged', 'SCF run NOT converged'),
                        log('UKS', -9.8).replace('PROGRAM ENDED AT now', ''),
                        log('UKS', 'NaN')]:
            with self.subTest(content=content):
                p.write_text(content)
                with self.assertRaises(ValueError):
                    parse_output(p, 'UKS')

    def test_fortran_exponent(self):
        p = self.root / 'test.out'
        p.write_text(log('ROKS', '-1.0D+01'))
        self.assertEqual(parse_output(p, 'ROKS'), -10.0)

    def test_missing_triplet_preserves_existing_csv(self):
        p = self.case(dataset='QUEST6')
        (p / 'triplet.out').unlink()
        output = self.root / 'existing.csv'
        output.write_text('previous results\n')
        with self.assertRaises(SystemExit) as exc:
            main(['--data-dir', str(self.root), '--output', str(output)])
        self.assertEqual(exc.exception.code, 1)
        self.assertEqual(output.read_text(), 'previous results\n')


if __name__ == '__main__':
    unittest.main()
