import unittest

import editdistance
import numpy as np
from wub.simulate import seq as sim_seq
from wub.util import seq as seq_util


class TestSimulateSeq(unittest.TestCase):

    """Test sequence simulation utilities."""

    def test_simulate_sequencing_errors(self):
        """Test function simulating sequencing errors."""
        error_rate = 0.1
        error_weights = {'substitution': 1.0 / 6,
                         'insertion': 1.0 / 6,
                         'deletion': 4.0 / 6}
        sequence = sim_seq.simulate_sequence(5000)
        mutated_record = sim_seq.simulate_sequencing_errors(
            sequence, error_rate, error_weights)
        distance = editdistance.eval(sequence, mutated_record.seq)
        expected_errors = len(sequence) * error_rate
        errors_sd = np.sqrt(len(sequence) * error_rate * (1 - error_rate))
        # Should pass 0.9973 proportion of cases:
        self.assertTrue(expected_errors - errors_sd * 3 < distance < expected_errors +
                        errors_sd * 3, msg="expected: {} realised:{}".format(expected_errors, distance))

    def test_add_errors(self):
        """Test function adding sequencing errors."""
        seq = "ATGCATGCATGC"
        mut_seq = sim_seq.add_errors(seq, 6, 'substitution')
        self.assertSequenceEqual(seq_util.alignment_stats(seq, mut_seq), (12, 6, 0, 0, 0.5))

    def test_compress_raw_cigar_list(self):
        """Test compression of raw cigar lists."""
        cigar_list = [
            (1, 'M'), (1, 'M'), (1, 'M'), (1, 'D'), (1, 'D'), (1, 'M'), (1, 'I'), (1, 'M')]
        compressed = sim_seq.compress_raw_cigar_list(cigar_list)
        expected = [(3, 'M'), (2, 'D'), (1, 'M'), (1, 'I'), (1, 'M')]
        self.assertSequenceEqual(compressed, expected)

    def test_cigar_list_to_string(self):
        """Test formatting of cigar strings."""
        cigar_list = [(3, 'M'), (2, 'D'), (1, 'M'), (1, 'I'), (1, 'M')]
        cigar_string = sim_seq.cigar_list_to_string(cigar_list)
        expected = "3M2D1M1I1M"
        self.assertEqual(cigar_string, expected)
