import unittest

from wub.mappers import lastal
from wub.util import seq as seq_util

class TestMappersLastal(unittest.TestCase):

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
