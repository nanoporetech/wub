import unittest

import editdistance
from wub.simulate import seq as sim_seq
from Bio.SeqRecord import SeqRecord


class TestSimulateSeq(unittest.TestCase):

    """Test sequence utilities."""

    def test_simulate_sequencing_errors(self):
        """Test function simulating sequencing errors."""
        error_rate = 0.1
        error_weights = {'substitution': 1.0 / 6, 'insertion': 1.0 / 6, 'deletion': 4.0 / 6}
        sequence = sim_seq.simulate_sequence(5000)
        mutated_record = sim_seq.simulate_sequencing_errors(sequence, error_rate, error_weights)
        distance = editdistance.eval(sequence, mutated_record.seq)
        expected_errors = len(sequence) * error_rate
        errors_sd = len(sequence) * error_rate * (1 - error_rate)
        self.assertTrue(expected_errors - errors_sd < distance < expected_errors + errors_sd, msg="expected: {} realised:{}".format(expected_errors, distance))

