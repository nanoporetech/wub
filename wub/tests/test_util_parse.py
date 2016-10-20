import unittest

from wub.util import parse
import numpy as np


class TestUtilParse(unittest.TestCase):

    """Test parsing utilities."""

    def test_separated_list_to_floats(self):
        """Test parsing of separated lists."""
        string = "0.1,0.2,0.3"
        parsed = (0.1, 0.2, 0.3)
        self.assertSequenceEqual(parse.separated_list_to_floats(string), parsed)

    def test_args_string_to_dict(self):
        """Test parsing of dictionaries encoded in separated strings."""
        string = "a:0.1,b:0.2,c:0.3"
        parsed = (("a", "0.1"), ("b", "0.2"), ("c", "0.3"))
        self.assertSequenceEqual(parse.args_string_to_dict(string).items(), parsed)

    def test_normalise_array(self):
        """Test array normalization."""
        a = np.array([2, 2, 2, 2])
        a_norm = np.array([0.25, 0.25, 0.25, 0.25], dtype=float)
        self.assertTrue(all(parse.normalise_array(a) == a_norm))
