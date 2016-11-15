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
        self.assertEqual(res['read_stats'],{'unmapped': 0, 'alignment_lengths': [87], 'aligned_quals': [40], 'aligned_lengths': [87], 'mapped': 1})


