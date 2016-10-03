#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

from Bio import SeqIO
import numpy as np

from wub.simulate import seq as sim_seq
from wub.simulate import genome as sim_genome
from wub.util import parse as parse_util
from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Sample fragments from the input genome and simulate sequencing errors. Read lengths are drawn from the specified truncated gamma distribution. Chromosomes are sampled randomly for each read.

    The format of the read names is the following:
    r<unique_id>_<chromosome>_<frag_start>_<frag_end>/q<realised_quality>/s<realised_substiutions>/d<realised_deletions>/i<realised_insertions>
    """)
parser.add_argument(
    '-n', metavar='nr_reads', type=int, help="Number of simulated reads (1).", default=1)
parser.add_argument('-m', metavar='mean_length', type=int,
                    help="Mean read length (5000).", default=5000)
parser.add_argument(
    '-a', metavar='gamma_shape', type=float, help="Read length distribution: gamma shape parameter (1).", default=1.0)
parser.add_argument(
    '-l', metavar='low_trunc', type=int, help="Read length distribution: lower truncation point (None).", default=None)
parser.add_argument(
    '-u', metavar='high_trunc', type=int, help="Read length distribution: upper truncation point (None).", default=None)
parser.add_argument(
    '-e', metavar='error_rate', type=float, help="Total rate of substitutions insertions and deletions (0.1).", default=0.1)
parser.add_argument('-w', metavar='error_weights', type=str,
                    help="Relative frequency of substitutions,insertions,deletions (1,1,4).", default="1,1,4")
parser.add_argument(
    '-q', metavar='mock_quality', type=int, help="Mock base quality for fastq output (40).", default=40)
parser.add_argument('input_fasta', nargs='?', help='Input genome in fasta format (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_fastq', nargs='?', help='Output fastq (default: stdout)',
                    type=argparse.FileType('w'), default=sys.stdout)

# Possible extension to sequencing error simulation:
# - Save true alignment in SAM format.
# - Biased error probabilities and biased insert composition.
# - Simulate insert length from distribution.
# - Save simulation details to CSV file.
# - Simulate deletions longer than one? Does it make sense? Might be challenging.


def simulate_sequencing(chromosomes, mean_length, gamma_shape, low_truncation, high_truncation, mock_quality, number_reads):
    """Simulate sequenceing.
    :param chromosomes: Input chromosomes (list of SeqRecord obejcts).
    :param mean_length: Mean read length.
    :param gamma_shape: Shape paramter of the read length distribution.
    :param low_truncation: Minimum read length.
    :param high_truncation: Maximum read length.
    :param mock_quality: Mock base quality for fastq output.
    :param number_reads: Number of reads to simulate.
    """
    for fragment in sim_genome.simulate_fragments(chromosomes, mean_length, gamma_shape, low_truncation, high_truncation, number_reads):
        # Simulate sequenceing errors:
        mutated_record = sim_seq.simulate_sequencing_errors(
            fragment.seq, args.e, error_weights)
        # Construct read name:
        read_name = "r{}_{}_{}_{}".format(
            fragment.uid, fragment.chrom, fragment.start, fragment.end)
        read_name = "{}/q{}/s{}/d{}/i{}".format(read_name, mutated_record.real_qual,
                                                mutated_record.real_subst, mutated_record.real_del, mutated_record.real_ins)

        yield seq_util.new_dna_record(mutated_record.seq, read_name, [mock_quality] * len(mutated_record.seq))

if __name__ == '__main__':
    args = parser.parse_args()

    # Read in chromosomes of the input genome:
    chromosomes = list(seq_util.read_seq_records(args.input_fasta))

    # Process error weights:
    error_weights = np.array(parse_util.separated_list_to_floats(args.w))
    # Normalise relative base frequencies to probabilities:
    error_weights = parse_util.normalise_array(error_weights)
    error_weights = dict(
        zip(['substitution', 'insertion', 'deletion'], error_weights))

    simulation_iterator = simulate_sequencing(
        chromosomes, args.m, args.a, args.l, args.u, args.q, args.n)

    seq_util.write_seq_records(
        simulation_iterator, args.output_fastq, format='fastq')
