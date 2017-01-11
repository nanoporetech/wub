# -*- coding: utf-8 -*-

import numpy as np
import functools
from collections import namedtuple

from wub.util import seq as seq_util

uniform_probs = [0.25, 0.25, 0.25, 0.25]

strand_directions = ['+', '-']

MutatedSeq = namedtuple(
    'MutatedSeq', 'seq real_qual real_subst real_del real_ins cigar')

cigar_operations = {'match': 'M', 'substitution': 'M', 'insertion': 'I', 'deletion': 'D'}


def sample_direction(forward_prob):
    return np.random.choice(strand_directions, p=[forward_prob, 1 - forward_prob])


def random_base(probs=uniform_probs):
    """Generate a random DNA base.

    :param probs: Probabilities of sampling a base, in the ACGT order.
    :returns: A sampled base.
    :rtype: str
    """
    return np.random.choice(seq_util.bases, p=probs)


def random_base_except(excluded, probs=uniform_probs):
    """Generate a random base according to the specified probabilities with the exclusion of the specified base.

    :param excluded: Exclude this base from sampling.
    :param probs: Base sampling probabilities in the ACGT order.
    :returns: A sampled base.
    :rtype: str
    """
    if len(probs) != len(seq_util.bases):
        raise ValueError('Probability vector has wrong length!')
    # Filter out excluded base:
    bp_dict = dict((x, y)
                   for x, y in zip(seq_util.bases, probs) if x != excluded)
    filtered_bases = bp_dict.keys()
    norm_probs = np.array(bp_dict.values(), dtype=float)
    # Re-normalise probabilities:
    norm_probs = norm_probs / np.sum(norm_probs)
    return np.random.choice(filtered_bases, p=norm_probs)


def simulate_sequence(length, probs=uniform_probs):
    """Simulate sequence of specified length and base composition.

    :param length: Length of simulated sequence.
    :param probs: Base composition vector in the ACGT order.
    :returns: Simulated sequence.
    :rtype: str
    """
    return ''.join(np.random.choice(seq_util.bases, size=length, p=probs))


def sample_error_type(error_weights):
    """Sample error type from error weights dictionary.

    :param error_weights: A dcitionary with (type, probability) pairs.
    :returns: Error type
    :rtype: str
    """
    return np.random.choice(error_weights.keys(), p=error_weights.values())


def cigar_list_to_string(cigar_list):
    """Sample error type from error weights dictionary.

    :param error_weights: A dcitionary with (type, probability) pairs.
    :returns: Error type
    :rtype: str
    """

    tmp = map(lambda x: str(x[0]) + str(x[1]), cigar_list)
    return ''.join(tmp)


def compress_raw_cigar_list(raw_cigar):
    """Sample error type from error weights dictionary.

    :param error_weights: A dcitionary with (type, probability) pairs.
    :returns: Error type
    :rtype: str
    """

    raw_cigar[0] = [raw_cigar[0]]

    def cigar_op_compose(a, b):
        x = a.pop()
        if x[1] == b[1]:
            a.append((x[0] + b[0], x[1]))
        else:
            a.extend([x, b])
        return a

    cigar = functools.reduce(cigar_op_compose, raw_cigar)
    return cigar


def simulate_sequencing_errors(sequence, error_rate, error_weights):
    """Simulate substitutions, deletions and insertions.

    :param sequence: Input sequence.
    :param error_rate: Total error rate.
    :param error_weights: A dictionary with error types as keys and probabilities as values.
    The possible error types are: substitution, deletion, insertion.
    :returns: A named tuple with elements: mutated sequence, realised quality, number of realised substitutions,
    number of realised deletions, number of realised insertions, cigar string.
    :rtype: namedtuple
    """
    if len(sequence) == 0:
        raise Exception('Cannot simulate sequencing errors on empty sequence!')

    new_bases = []

    realised_substitutions = 0
    realised_deletions = 0
    realised_insertions = 0
    raw_cigar_list = []

    for position, base in enumerate(sequence):
        if np.random.uniform() < error_rate:
            error_type = sample_error_type(error_weights)

            if error_type == 'substitution':
                new_base = random_base_except(base)
                realised_substitutions += 1
                raw_cigar_list.append((1, cigar_operations[error_type]))

            elif error_type == 'deletion':
                new_base = ''
                realised_deletions += 1
                raw_cigar_list.append((1, cigar_operations[error_type]))

            elif error_type == 'insertion':
                new_base = base + random_base()
                realised_insertions += 1
                raw_cigar_list.append((1, cigar_operations['match']))
                raw_cigar_list.append((1, cigar_operations[error_type]))

            else:
                raise Exception("Unhandled error type: {}".format(error_type))
        else:
            raw_cigar_list.append((1, cigar_operations['match']))
            new_base = base
        new_bases.append(new_base)

    new_sequence = ''.join(new_bases)
    cigar = cigar_list_to_string(compress_raw_cigar_list(raw_cigar_list))

    realised_events = realised_substitutions + \
        realised_deletions + realised_insertions
    realised_quality = seq_util.prob_to_phred(
        round(float(realised_events) / float(len(sequence)), 3))
    mutated_record = MutatedSeq(
        new_sequence, realised_quality, realised_substitutions, realised_deletions, realised_insertions, cigar)
    return mutated_record


def add_errors(seq, nr_errors, error_type):
    """Introduce a specified number of errors in the target sequence at random positions.

    :param seq: Input DNA sequence.
    :param nr_errors: Number of mismatches to introduce.
    :returns: Mutated sequence.
    :rtype: str
    """
    seq = list(seq)
    positions = np.random.choice(np.arange(len(seq)), size=nr_errors, replace=False)
    if error_type == 'substitution':
        for pos in positions:
            seq[pos] = random_base_except(seq[pos])
    elif error_type == 'deletion':
        for pos in positions:
            seq[pos] = ''
    elif error_type == 'insertion':
        for pos in positions:
            seq[pos] = seq[pos] + random_base()
    else:
        raise Exception('Invalid error type')
    return ''.join(seq)
