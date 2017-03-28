# -*- coding: utf-8 -*-

import sys
from itertools import izip
import itertools
from collections import namedtuple, OrderedDict, Counter
import numpy as np
from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA, IUPACAmbiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import pysam

""" Utilities manipulating biological sequences and formats. Extensions to biopython functionality.
"""

# Reverse complements of bases, taken from dragonet:
comp = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'X': 'X', 'N': 'N',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'x': 'x', 'n': 'n',
    '-': '-'
}

# Shortcut to the list of DNA bases:
bases = sorted(list(IUPACUnambiguousDNA().letters))
ambiguous_bases = sorted(list(IUPACAmbiguousDNA().letters))


def base_complement(k):
    """ Return complement of base.

    Performs the subsitutions: A<=>T, C<=>G, X=>X for both upper and lower
    case. The return value is identical to the argument for all other values.

    :param k: A base.
    :returns: Complement of base.
    :rtype: str

    """
    try:
        return comp[k]
    except KeyError:
        sys.stderr.write(
            "WARNING: No reverse complement for {} found, returning argument.".format(k))
        return k


def reverse_complement(seq):
    """ Return reverse complement of a string (base) sequence.

    :param seq: Input sequence.
    :returns: Reverse complement of input sequence.
    :rtype: str

    """
    if len(seq) == 0:
        return seq
    return reduce(lambda x, y: x + y, map(base_complement, seq[::-1]))


def base_composition(seq):
    """ Return letter counts of a string (base) sequence.

    :param seq: Input sequence.
    :returns: Letter counts.
    :rtype: dict

    """
    return dict(Counter(seq))


def word_composition(seq, size):
    """ Return word counts of a nucleotide sequence.

    :param seq: Input sequence.
    :param size: word length.
    :returns: word counts.
    :rtype: OrderedDict

    """
    composition = OrderedDict()
    for word in itertools.product(*([bases] * size)):
        word = ''.join(word)
        composition[word] = seq.count(word)
    return composition


def gc_content(seq):
    """ Return fraction of GC bases in sequence.

    :param seq: Input sequence.
    :returns: GC content.
    :rtype: float

    """
    return (seq.count('G') + seq.count('g') + seq.count('C') + seq.count('c')) / float(len(seq))


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


def new_dna_record(sequence, name, qualities=None):
    """Create a new SeqRecord object using IUPACUnambiguousDNA and the specified sequence.

    :param sequence: The sequence.
    :param name: Record identifier.
    :param qualities: List of base qualities.
    :returns: The SeqRecord object.
    :rtype: SeqRecord

    """
    seq_record = SeqRecord(
        Seq(sequence, IUPACUnambiguousDNA), id=name, description="", name="")
    if qualities is not None:
        seq_record.letter_annotations["phred_quality"] = qualities
    return seq_record


def write_seq_records(records_iterator, output_object, format='fasta'):
    """Write out SeqRecord objects to a file from an iterator in the specified format.

    :param records_iterator: An iterator of SeqRecord objects.
    :param output_object: Open file object or file name.
    :param format: Output format (fasta by default).
    :returns: None
    :rtype: object

    """
    if type(output_object) == file:
        SeqIO.write(records_iterator, output_object, format)
    else:
        with open(output_object, 'w') as output_handle:
            SeqIO.write(records_iterator, output_handle, format)


def read_seq_records(input_object, format='fasta'):
    """Read SeqRecord objects from a file in the specified format.

    :param input_object: A file object or a file name.
    :param format: Input format (fasta by default).
    :returns: A dictionary with the parsed SeqRecord objects.
    :rtype: generator

    """
    handle = input_object
    if type(handle) != file:
        handle = open(handle, "rU")
    return SeqIO.parse(handle, format)


def read_seq_records_dict(input_object, format='fasta'):
    """Read SeqRecord objects to a dictionary from a file in the specified format.

    :param input_object: A file object or a file name.
    :param format: Input format (fasta by default).
    :returns: An iterator of SeqRecord objects.
    :rtype: dict

    """
    handle = input_object
    if type(handle) != file:
        handle = open(handle, "rU")
    return SeqIO.to_dict(SeqIO.parse(handle, format))


def record_lengths(input_iter):
    """Return lengths of SeqRecord obejcts in the input iterator.

    :param input_iter: An iterator of SeqRecord objects.
    :returns: An ordered dictionary with the lengths of the SeqRecord objects.
    :rtype: OrderedDict

    """
    lengths = OrderedDict((record.id, len(record)) for record in input_iter)
    return lengths


def prob_to_phred(error_prob, max_q=93, qround=True):
    """Convert error probability into phred score.

    :param error_prob: Base error probability.
    :param max_q: Maximum quality value.
    :param qround: Round calculated score.
    :returns: Phred score.
    :rtype: int
    """
    if error_prob == 0:
        return max_q
    q = -10 * np.log10(error_prob)
    if qround:
        return int(round(min(max_q, q)))
    else:
        return min(max_q, q)


def phred_to_prob(phred):
    """Convert phred score into error probability.

    :param phred: Phred quality score.
    :returns: Error probability.
    :rtype: float
    """
    return np.power(10, -phred / 10.0)


def mean_qscore(scores, qround=True):
    """ Returns the phred score corresponding to the mean of the probabilities associated with the phred scores provided.

    :param scores: Iterable of phred scores.
    :param qround: Round after calculating mean score.
    :returns: Phred score corresponding to the average error rate, as estimated from the input phred scores.
    """
    if len(scores) == 0:
        return 0.0
    sum_prob = 0.0
    for val in scores:
        sum_prob += phred_to_prob(val)
    mean_prob = sum_prob / float(len(scores))
    mean_phred = prob_to_phred(mean_prob, qround=qround)
    return mean_phred


def quality_array_to_string(quality_list):
    """Convert list of phred quality values to string.

    :param quality_list: List of phred quality scores.
    :returns: Quality string.
    :rtype: str
    """
    return pysam.qualities_to_qualitystring(quality_list)


def quality_string_to_array(quality_string):
    """Convert quality string into a list of phred scores.

    :param quality_string: Quality string.
    :returns: Array of scores.
    :rtype: array
    """
    return pysam.qualitystring_to_array(quality_string)


def read_alignment(input_file, format='fasta'):
    """
    Load multiple alignment from file.

    :param input_file: Input file name.
    :returns: The alignment read from the input file.
    :rtype: MultipleSeqAlignment

    """
    msa = AlignIO.read(input_file, format)
    return msa


def alignment_stats(ref, query, gap_character='-'):
    """
    Calculate statistics from two aligned sequences.

    :param ref: Reference sequence.
    :param query: Query sequence.
    :param gap_character: Gap symbol.
    :returns: AlnStats namedtuple.
    :rtype: namedtuple
    """
    if len(ref) != len(query):
        raise Exception('Aligned sequences differ in length!')
    AlnStats = namedtuple(
        'AlnStats', 'length substitutions deletions insertions accuracy')
    substitutions = 0
    deletions = 0
    insertions = 0
    for ref_sym, query_sym in izip(ref, query):
        if ref_sym != query_sym:
            if query_sym == gap_character:
                deletions += 1
            elif ref_sym == gap_character:
                insertions += 1
            else:
                substitutions += 1
    accuracy = 1 - ((substitutions + deletions + insertions) / float(len(ref)))
    return AlnStats(len(ref), substitutions, deletions, insertions, accuracy)
