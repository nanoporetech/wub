import unittest

from os import path
from wub.parsers import blastn


class TestBlastnCoordParse(unittest.TestCase):

    def test_nucmer_coord_parse(self):
        """Test blastn outfmt 6 cooridnate parsing."""
        top = path.dirname(__file__)
        coord_file = path.join(top, "data/test_blastn_parse/blastn_test.coords")
        records = blastn.parse_coords(coord_file)
        self.assertEqual(records,  [{'gapopen': 0, 'query_end': 100, 'mismatch': 0, 'ref_end': 300, 'query': 'seq_0', 'identity': 100.0, 'bitscore': 181.0, 'query_start': 1, 'ref_start': 201, 'strand': '+', 'aln_length': 100, 'ref': 'seq_0', 'evalue': 3e-50}, {'gapopen': 0, 'query_end': 100, 'mismatch': 0, 'ref_end': 600, 'query': 'seq_0', 'identity': 100.0, 'bitscore': 181.0, 'query_start': 1, 'ref_start': 501, 'strand': '-', 'aln_length': 100, 'ref': 'seq_0', 'evalue': 3e-50}])
