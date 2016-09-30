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
    return SeqRecord(Seq(sequence, IUPACUnambiguousDNA), id=name, description="", name="")

def write_seq_records(records_iterator, output_file, format='fasta'):
    """Write out SeqRecord objects to a file from an iterator in the specified format.

    :param records_iterator: An iterator of SeqRecord objects.
    :param output_file: Output file.
    :param format: Output format (fasta by default).
    :returns: None
    :rtype: object

    """
    with open(output_file, 'w') as output_handle:
        SeqIO.write(records_iterator, output_handle, format)
