import gzip
import hashlib
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'audit.py'
spec = importlib.util.spec_from_file_location('audit', SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class CoverageTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.f = Path(self.tmp.name) / 'ref.fasta'
    def write(self, text):
        self.f.write_text(text)
    def test_wrapped_exact_and_coordinates(self):
        self.write('>p\nAAKLH\nGILVEAAA\n')
        z = mod.audit(self.f, ['KLHGILVEA'])['peptides'][0]
        self.assertEqual(z['occurrences'], 1)
        self.assertEqual(z['matches'][0]['start_1based'], 3)
        self.assertEqual(z['matches'][0]['end_1based_inclusive'], 11)
    def test_absence_requires_hash(self):
        self.write('>p\nAAA\n')
        self.assertEqual(mod.audit(self.f, ['CCC'])['peptides'][0]['status'],
                         'no_exact_match_in_unverified_or_ambiguous_reference')
        sha = hashlib.sha256(self.f.read_bytes()).hexdigest()
        self.assertEqual(mod.audit(self.f, ['CCC'], sha)['peptides'][0]['status'],
                         'not_in_verified_reference')
    def test_hash_mismatch(self):
        self.write('>p\nAAA\n')
        with self.assertRaises(ValueError): mod.audit(self.f, ['AAA'], '0'*64)
    def test_truncated_size(self):
        self.write('>p\nAAA\n')
        with self.assertRaises(ValueError): mod.audit(self.f, ['AAA'], expected_bytes=100)
    def test_no_cross_record_or_stop_matching(self):
        self.write('>a\nKLHG\n>b\nILVEA\n>c\nKLHG*ILVEA\n')
        self.assertEqual(mod.audit(self.f, ['KLHGILVEA'])['peptides'][0]['occurrences'], 0)
    def test_overlapping_matches_and_cap(self):
        self.write('>p\nAAAAA\n')
        z = mod.audit(self.f, ['AAA'], max_hits=1)['peptides'][0]
        self.assertEqual(z['occurrences'], 3)
        self.assertTrue(z['matches_truncated'])
    def test_ambiguous_reference(self):
        self.write('>p\nAAXAA\n')
        sha = hashlib.sha256(self.f.read_bytes()).hexdigest()
        self.assertEqual(mod.audit(self.f, ['CCC'], sha)['peptides'][0]['status'],
                         'no_exact_match_in_unverified_or_ambiguous_reference')
    def test_il_distinction(self):
        self.write('>p\nILMILQPQL\n')
        self.assertEqual(mod.audit(self.f, ['LLMILQPQL'])['peptides'][0]['occurrences'], 0)
    def test_malformed_fasta(self):
        for text in ('AAA\n', '>p\n', '>p\nAAA123\n', '>\nAAA\n', ''):
            self.write(text)
            with self.assertRaises(ValueError): mod.audit(self.f, ['AAA'])
    def test_invalid_peptide(self):
        self.write('>p\nAAA\n')
        for p in ('aaa', 'AA[+80]', 'AAX', ''):
            with self.assertRaises(ValueError): mod.audit(self.f, [p])
    def test_gzip(self):
        self.f = self.f.with_suffix('.fasta.gz')
        with gzip.open(self.f, 'wt') as out: out.write('>p\nAAA\n')
        self.assertEqual(mod.audit(self.f, ['AAA'])['peptides'][0]['occurrences'], 1)
    def test_cli_error_is_nonzero_without_json(self):
        self.write('>p\nAAA\n')
        z = subprocess.run([sys.executable, str(SCRIPT), '--fasta', str(self.f),
              '--peptide', 'AAA', '--expected-sha256', '0'*64], capture_output=True)
        self.assertEqual(z.returncode, 2)
        self.assertEqual(z.stdout, b'')

if __name__ == '__main__': unittest.main()
