# -*- coding: utf-8 -*-

import six
import sys
import numpy as np
from collections import OrderedDict, namedtuple

from wub.util import seq as seq_util
from wub.simulate import seq as sim_seq
from wub.simulate import dist

Fragment = namedtuple('Fragment', 'chrom uid start end seq')


def simulate_genome(number_chromosomes, mean_length, gamma_shape, low_truncation, high_truncation, base_frequencies):
    """Generator function for simulating chromosomes in a genome.
    Chromosome lengths are sampled from a truncated gamma distribution.

    :param number_chromosomes: Number of simulated chromosomes.
    :param mean_length: Mean length of simulated chromosomes.
    :param gamma_shape: Shape parameter of the chromosome length distribution.
    :param low_truncation: Minimum chromosome length.
    :param high_truncation: Maximum chromosome length.
    :param base_frequencies: Array of base frequencies in the ACGT order.
    :returns: A generator of SeqRecord objects.
    :rtype: generator

    """
    chrom_info = OrderedDict(
        ('chr' + str(i),
         int(dist.sample_truncated_gamma(mean_length, gamma_shape, low_truncation, high_truncation)))
        for i in six.xrange(number_chromosomes))
    sim_iter = (seq_util.new_dna_record(sim_seq.simulate_sequence(length, base_frequencies), name)
                for name, length in six.iteritems(chrom_info))
    return sim_iter


def sample_chromosome(chromosomes):
    """Sample a random chromosome.

    :param chromosomes: A collection of SeqRecord object.
    :returns: A randomly sampled element from the input collection.
    :rtype: SeqRecord
    """
    indexes = range(len(chromosomes))
    pick = np.random.choice(indexes)
    return chromosomes[pick]


def simulate_fragment(chromosome, mean_length, gamma_shape, low_truncation, high_truncation, fragment_number):
    """Simulate a fragment from a chromosome.

    :param chromosome: Chromosome to simulate fragment from, SeqRecord object.
    :param mean_length: Mean length of simulated fragment.
    :param gamma_shape: Shape parameter of length distribution.
    :param low_truncation: Minimum read length.
    :param high_truncation: Maximum read length.
    :param fragment_number: The unique identifier of fragment in simulation (number of fragment).
    :returns: A named tuple with chromosome id, fragment number, start, end and sequence.
    :rtype: namedtuple
    """
    fragment_length = int(dist.sample_truncated_gamma(
        mean_length, gamma_shape, low_truncation, high_truncation))
    upper_boundary = len(chromosome) - fragment_length
    # Special case when upper boundary is less than the read length. Maybe
    # should handle this by rejection?
    if upper_boundary < fragment_length:
        start = 0
        end = len(chromosome)
    else:
        start = np.random.randint(0, upper_boundary)
        end = start + fragment_length
    fragment_sequence = chromosome.seq[start:end]
    return Fragment(chromosome.id, fragment_number, start, end, fragment_sequence)


def simulate_fragments(chromosomes, mean_length, gamma_shape, low_truncation, high_truncation, number_fragments):
    """Simulate a fragments from a set of chromosomes. Chromosomes are picked randomly for each fragment.

    :param chromosomes: Chromosomes to simulate fragment from, a list of SeqRecord objects.
    :param mean_length: Mean length of simulated fragments.
    :param gamma_shape: Shape parameter of length distribution.
    :param low_truncation: Minimum read length.
    :param high_truncation: Maximum read length.
    :param number_fragments: Number of fragments to simulate.
    :returns: An iterator named tuples with chromosome id, fragment number, start, end and sequence.
    :rtype: generator
    """
    fragment_uid = 0
    while True:
        if fragment_uid >= number_fragments:
            break
        chromosome = sample_chromosome(chromosomes)
        fragment = simulate_fragment(
            chromosome, mean_length, gamma_shape, low_truncation, high_truncation, fragment_uid)
        if (fragment.end - fragment.start) > 0:
            fragment_uid += 1
            yield fragment
        else:
            sys.stderr.write(
                "Skipped zero length fragment! Consider increase minimum read length!\n")
