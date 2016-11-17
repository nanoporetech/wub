import unittest

from os import path
from wub.bam import stats
from wub.util import seq as seq_util


class TestBamStats(unittest.TestCase):

    """Test BAM statistics functions."""

    def test_error_and_read_stats(self):
        """Test the gathering of error and read statistics."""
        top = path.dirname(__file__)
        ref_fasta = path.join(top, "data/test_bam_stats/stat_ref.fas")
        bam = path.join(top, "data/test_bam_stats/stat_test.bam")
        refs = seq_util.read_seq_records_dict(ref_fasta)
        res = stats.error_and_read_stats(bam, refs, context_sizes=(1, 1), region=None, min_aqual=0)

        # Test evenets:
        self.assertEqual(res['events']['AGA'], {'*': 1, 'G': 2})
        self.assertEqual(res['events']['CGA'], {'-': 1, 'G': 2})
        self.assertEqual(res['events']['ACA'], {'C': 2, 'T': 1})

        # Test indel properties:
        self.assertEqual(res['indel_dists']['insertion_lengths'], {8: 1})
        self.assertEqual(res['indel_dists']['insertion_composition'], {'G': 8})
        self.assertEqual(res['indel_dists']['deletion_lengths'], {9: 1})

        # Test read statistics:
        self.assertEqual(res['read_stats'], {'alignment_lengths': [87], 'mapping_quals': [47], 'unaligned_lengths': [], 'unaligned_quals': [], 'mqfail_alignment_lengths': [], 'mapped': 1, 'unmapped': 0, 'mqfail_aligned_quals': [], 'aligned_quals': [40], 'aligned_lengths': [87]})

    def test_read_stats(self):
        """Test the gathering read statistics."""
        top = path.dirname(__file__)
        bam = path.join(top, "data/test_bam_stats/stat_test.bam")
        res = stats.read_stats(bam, region=None, min_aqual=0)

        self.assertEqual(res, {'alignment_lengths': [87], 'mapping_quals': [47], 'unaligned_lengths': [], 'unaligned_quals': [], 'mqfail_alignment_lengths': [], 'mapped': 1, 'unmapped': 0, 'mqfail_aligned_quals': [], 'aligned_quals': [40], 'aligned_lengths': [87]})

    def test_pileup_stats(self):
        """Test the gathering read statistics."""
        top = path.dirname(__file__)
        bam = path.join(top, "data/test_bam_stats/stat_test.bam")
        res = stats.pileup_stats(bam, region=None)

        self.assertEqual(res, {'coverage': {'seq_0': {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1, 11: 1, 12: 1, 13: 1, 14: 1, 15: 1, 16: 1, 17: 1, 18: 1, 19: 1, 20: 1, 21: 1, 22: 1, 23: 1, 24: 1, 25: 1, 26: 1, 27: 1, 28: 1, 29: 1, 30: 1, 31: 1, 32: 1, 33: 1, 34: 1, 35: 1, 36: 1, 37: 1, 38: 1, 39: 1, 40: 1, 41: 1, 42: 1, 43: 1, 44: 1, 45: 1, 46: 1, 47: 1, 48: 1, 49: 1, 50: 1, 51: 1, 52: 1, 53: 1, 54: 1, 55: 1, 56: 1, 57: 1, 58: 1, 59: 1, 60: 1, 61: 1, 62: 1, 63: 1, 64: 1, 65: 1, 66: 1, 67: 1, 68: 1, 69: 1, 70: 1, 71: 1, 72: 1, 73: 1, 74: 1, 75: 1, 76: 1, 77: 1, 78: 1, 79: 1, 80: 1, 81: 1, 82: 1, 83: 1, 84: 1, 85: 1, 86: 1, 87: 1}}, 'qualities': {'seq_0': {0: [40], 1: [40], 2: [40], 3: [40], 4: [40], 5: [
                         40], 6: [40], 7: [40], 8: [40], 9: [40], 10: [40], 11: [40], 12: [40], 13: [40], 14: [40], 15: [40], 16: [40], 17: [40], 18: [40], 19: [40], 20: [40], 21: [40], 22: [40], 23: [40], 24: [40], 34: [40], 35: [40], 36: [40], 37: [40], 38: [40], 39: [40], 40: [40], 41: [40], 42: [40], 43: [40], 44: [40], 45: [40], 46: [40], 47: [40], 48: [40], 49: [40], 50: [40], 51: [40], 52: [40], 53: [40], 54: [40], 55: [40], 56: [40], 57: [40], 58: [40], 59: [40], 60: [40], 61: [40], 62: [40], 63: [40], 64: [40], 65: [40], 66: [40], 67: [40], 68: [40], 69: [40], 70: [40], 71: [40], 72: [40], 73: [40], 74: [40], 75: [40], 76: [40], 77: [40], 78: [40], 79: [40], 80: [40], 81: [40], 82: [40], 83: [40], 84: [40], 85: [40], 86: [40], 87: [40]}}})
