# -*- coding: utf-8 -*-

import numpy as np

from wub.util import seq as seq_util

uniform_probs = [0.25, 0.25, 0.25, 0.25]


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
    bp_dict = dict((x, y) for x, y in zip(seq_util.bases, probs) if x != excluded)
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
