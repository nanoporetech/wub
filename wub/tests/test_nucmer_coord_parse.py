import unittest

import six
from os import path
from collections import OrderedDict
from Bio import SeqIO
from wub.parsers import mummer


class TestNucmerCoordParse(unittest.TestCase):

    """Test BAM statistics functions."""

    def test_nucmer_coord_parse(self):
        """Test parsing of nucmer coordinate files."""
        top = path.dirname(__file__)
        coord_file = path.join(top, "data/test_nucmer_parse/nucmer_test.coords")
        records = mummer.parse_coords(coord_file)
        self.assertEqual(records, [{'query_len': 49, 'query_end': 49, 'ref_end': 500, 'ref_len': 47, 'query': 'adapter', 'query_start': 1, 'ref_start': 454, 'ref': 'seq_0', 'identity': 90.0}, {'query_len': 48, 'query_end': 50, 'ref_end': 1003, 'ref_len': 48, 'query': 'adapter', 'query_start': 3, 'ref_start': 956, 'ref': 'seq_0', 'identity': 92.0}])
