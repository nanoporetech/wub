# -*- coding: utf-8 -*-

import numpy as np
from collections import OrderedDict

from wub.util import seq as seq_util
from wub.simulate import seq as sim_seq
from wub.simulate import dist


def simulate_genome(number_chromosomes, mean_length, gamma_shape, low_truncation, high_truncation, base_frequencies):
    """Generator function for simulating chromosomes in a genome. Chromosome lengths are sampled from a truncated gamma distribution.

    :param number_chromosomes: Number of simulated chromosomes.
    :param mean_length: Mean length of simulated chromosomes.
    :param gamma_shape: Shape parameter of the chromosome length distribution.
    :param low_truncation: Minimum chromosome length.
    :param high_truncation: Maximum chromosome length.
    :param base_frequencies: Comma-separated list of relative base frequencies in the ACGT order.
    :returns: A generator of SeqRecord objects.
    :rtype: generator

    """
    chrom_info = OrderedDict(
        ('chr' + str(i), int(dist.sample_truncated_gamma(mean_length, gamma_shape, low_truncation, high_truncation))) for i in xrange(number_chromosomes))
    sim_iter = (seq_util.new_dna_record(sim_seq.simulate_sequence(length, base_frequencies), name) for name, length in chrom_info.iteritems())
    return sim_iter
