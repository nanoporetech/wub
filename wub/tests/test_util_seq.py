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
        self.assertSequenceEqual(
            qual_seq.letter_annotations["phred_quality"], [mock_qual] * len(qual_seq))

    def test_reverse_complement(self):
        """Test reverse complementing."""
        sequence = "ATGCNXatgcnx-"
        revcomp = "-xngcatXNGCAT"
        self.assertEqual(seq.reverse_complement(sequence), revcomp)

    def test_prob_to_phred(self):
        """Test error probability to phred score conversion."""
        self.assertEqual(seq.prob_to_phred(0.5), 3)

    def test_prob_to_phred_max(self):
        """Test error probability to phred score conversion (very small error)."""
        self.assertEqual(seq.prob_to_phred(1 * 10 ** -10), 93)

    def test_phred_to_prob(self):
        """Test error probability to phred score conversion."""
        self.assertAlmostEqual(seq.phred_to_prob(3), 0.5, places=2)

    def test_mean_qscore(self):
        """Test mean q score calculation (large identical input)."""
        scores = [30] * 5000
        self.assertEqual(seq.mean_qscore(scores), 30)

    def test_mean_qscore(self):
        """Test mean q score calculation."""
        scores = [14, 10]
        self.assertEqual(seq.mean_qscore(scores), 12)

    def test_alignment_stats(self):
        """Test calculation of alignment statistics."""
        seq1 = "ATGCTG-AAAAA"
        seq2 = "TTG-TGCAAAAA"
        self.assertEqual(
            tuple(seq.alignment_stats(seq1, seq2)), (12, 1, 1, 1, 0.75))
