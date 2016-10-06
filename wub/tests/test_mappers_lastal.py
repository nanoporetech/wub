import unittest
import tempfile
import os
import numpy as np

from wub.mappers import lastal
from wub.util import seq as seq_util
from wub.simulate import seq as sim_seq

error_rate = 0.1
ref_length = 1000


class TestMappersLastal(unittest.TestCase):

    def setUp(self):
        fh_ref = tempfile.NamedTemporaryFile(suffix=".fas", delete=False)
        self.ref_fasta = fh_ref.name
        fh_target = tempfile.NamedTemporaryFile(suffix=".fas", delete=False)
        self.target_fasta = fh_target.name

        self.ref = sim_seq.simulate_sequence(ref_length)
        nr_errors = int(len(self.ref) * error_rate)
        self.target = sim_seq.add_mismatches(self.ref, nr_errors)

        left_flanking = sim_seq.simulate_sequence(50)
        right_flanking = sim_seq.simulate_sequence(50)

        self.ref = left_flanking + self.ref + right_flanking
        self.target = left_flanking + self.target + right_flanking

        fh_ref.write(">ref\n{}\n".format(self.ref))
        fh_ref.flush()
        fh_ref.close()

        fh_target.write(">target\n{}\n".format(self.target))
        fh_target.flush()
        fh_target.close()

    def tearDown(self):
        os.unlink(self.ref_fasta)
        os.unlink(self.target_fasta)

    def test_parse_lastal_identical(self):
        raw = """\
# batch 0
a score=23 EG2=3.8e+06 E=5.2e-13
s Simulomonas 0 23 + 23 ATGCGGGGGATAGGACCATATCT
s tig00000000 0 23 + 23 ATGCGGGGGATAGGACCATATCT
        """
        parsed = lastal.parse_lastal(raw).next()
        acc = seq_util.alignment_stats(parsed.r_aln, parsed.q_aln).accuracy
        self.assertEqual(acc, 1.0)
        self.assertEqual(parsed.score, 23)

    def test_parse_lastal_difference(self):
        raw = """\
# batch 0
a score=23 EG2=3.8e+06 E=5.2e-13
s Simulomonas 0 23 + 23 TTGCGGGGGATAGGACCATATCT
s tig00000000 0 23 + 23 ATGCGGGGGATAGGACCATATCT
        """
        parsed = lastal.parse_lastal(raw).next()
        acc = seq_util.alignment_stats(parsed.r_aln, parsed.q_aln).accuracy
        self.assertAlmostEqual(acc, 0.9565, places=3)

    def test_parse_lastal_zero(self):
        raw = """\
# batch 0
a score=23 EG2=3.8e+06 E=5.2e-13
s Simulomonas 0 23 + 23 CCCTCCCCCCCCCCCTTCCCCAC
s tig00000000 0 23 + 23 ATGCGGGGGATAGGACCATATCT
        """
        parsed = lastal.parse_lastal(raw).next()
        acc = seq_util.alignment_stats(parsed.r_aln, parsed.q_aln).accuracy
        self.assertAlmostEqual(acc, 0.0, places=3)

    def test_lastal_compare(self):
        substs = lastal.compare_genomes_lastal(
            self.ref_fasta, self.target_fasta)['substitutions'][0]
        self.assertEqual(int(ref_length * error_rate), substs)
