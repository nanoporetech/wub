# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA, IUPACAmbiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


""" Utilities manipulating biological sequences and formats. Extensions to biopython functionality.
"""

# Shortcut to the list of DNA bases:
bases = sorted(list(IUPACUnambiguousDNA().letters))
ambiguous_bases = sorted(list(IUPACAmbiguousDNA().letters))


def mock_qualities(record, mock_qual):
    """Add mock quality values to SeqRecord object.

    :param record: A SeqRecord object.
    :param mock_qual: Mock quality value used for each base.
    :returns: The record augmented with mock quality values.
    :rtype: object

    """
    rec_copy = record[:]
    rec_copy.letter_annotations["phred_quality"] = [mock_qual] * len(rec_copy)
    return rec_copy

def new_dna_record(sequence, name):
    """Create a new SeqRecord object using IUPACUnambiguousDNA and the specified sequence.

    :param sequence: The sequence.
    :param name: Record identifier.
    :returns: The SeqRecord object.
    :rtype: SeqRecord

    """
    return SeqRecord(Seq(sequence, IUPACUnambiguousDNA), id=name, description="")
