import unittest

from wub.util import seq
from Bio.SeqRecord import SeqRecord

class TestUtilSeq(unittest.TestCase):
    """Test sequence utilities."""

    def test_new_dna_record(self):
        """Test the construction of new DNA SeqRecord."""
        sequence = seq.new_dna_record("ATGC", "test")
        self.assertEqual(type(sequence), SeqRecord)

    def test_mock_qualities(self):
        """Test quality mocking function."""
        sequence = seq.new_dna_record("ATGC", "test")
        mock_qual = 40
        qual_seq = seq.mock_qualities(sequence, mock_qual)
        self.assertSequenceEqual(qual_seq.letter_annotations["phred_quality"], [mock_qual] * len(qual_seq))
