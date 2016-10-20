import unittest

import tempfile
import os

from wub.wrappers import dnadiff
from wub.util import seq as seq_util
from wub.util import cmd as cmd_util
from wub.simulate import seq as sim_seq

error_rate = 0.1
ref_length = 5000


class TestWrappersDnadiff(unittest.TestCase):

    """Test dnadiff wrapper."""

    def _generate_test_data(self):
        """Generate test data for dnadiff test."""
        fh_ref = tempfile.NamedTemporaryFile(suffix=".fas", delete=False)
        self.ref_fasta = fh_ref.name
        fh_target = tempfile.NamedTemporaryFile(suffix=".fas", delete=False)
        self.target_fasta = fh_target.name

        self.ref = sim_seq.simulate_sequence(ref_length)
        nr_errors = int(len(self.ref) * error_rate)
        self.target = sim_seq.add_errors(self.ref, nr_errors, 'substitution')

        fh_ref.write(">ref\n{}\n".format(self.ref))
        fh_ref.flush()
        fh_ref.close()

        fh_target.write(">target\n{}\n".format(self.target))
        fh_target.flush()
        fh_target.close()

    def _cleanup_test_data(self):
        """Cleanup test dataset."""
        os.unlink(self.ref_fasta)
        os.unlink(self.target_fasta)

    @unittest.skipIf(not cmd_util.find_executable('dnadiff'),
                     "Dnadiff binary not found, skipping integration tests.")
    def test_dnadiff(self):
        """Test dnadiff wrapper."""
        self._generate_test_data()
        res, _, _ = dnadiff.dnadiff(self.ref_fasta, self.target_fasta)
        self.assertAlmostEqual(
            res['Alignments']['1-to-1']['AvgIdentity'].ref, 90.0, places=0)
        self._cleanup_test_data()
