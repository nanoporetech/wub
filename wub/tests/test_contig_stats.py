import unittest

from wub.read_stats import contig_stats
import pandas as pd


class TestUtilSeq(unittest.TestCase):

    """Test N50 utility function."""

    def test_N50(self):
        """Test calculation of N50."""
        sequence_lengths = pd.DataFrame({'dummy': [2, 3, 4, 5, 6, 7, 8, 9, 10]}
        self.assertEqual(contig_stats.N50(sequence_lengths, 'dummy'), 8.0)
