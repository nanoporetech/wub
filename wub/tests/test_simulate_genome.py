import unittest

import editdistance
import numpy as np
from wub.simulate import genome as sim_genome
from wub.util import seq as seq_util


class TestSimulateGenome(unittest.TestCase):

    """Test genome simulation utilities."""

    def test_simulate_genome(self):
        """Test genome simulator."""
        record = sim_genome.simulate_genome(number_chromosomes=1, mean_length=1000, gamma_shape=50,
                                            low_truncation=1000, high_truncation=1001, base_frequencies=[0.25] * 4).next()
        self.assertEqual(len(record), 1000)

    def test_simulate_fragment(self):
        """Test fragment simulator."""
        chrom = sim_genome.simulate_genome(number_chromosomes=1, mean_length=1000, gamma_shape=50,
                                           low_truncation=1000, high_truncation=1001, base_frequencies=[0.25] * 4).next()
        frag = sim_genome. simulate_fragment(
            chrom, mean_length=50, gamma_shape=50, low_truncation=50, high_truncation=51, fragment_number=0)
        self.assertEqual(frag.end - frag.start, 50)

    def test_simulate_fragment_edge(self):
        """Test fragment simulator (edge case)."""
        chrom = sim_genome.simulate_genome(number_chromosomes=1, mean_length=1000, gamma_shape=50,
                                           low_truncation=1000, high_truncation=1001, base_frequencies=[0.25] * 4).next()
        frag = sim_genome. simulate_fragment(
            chrom, mean_length=2000, gamma_shape=50, low_truncation=2000, high_truncation=2001, fragment_number=0)
        self.assertEqual(frag.end - frag.start, 1000)
