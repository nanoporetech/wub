#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

import numpy as np
import pysam
from collections import OrderedDict

from wub.simulate import seq as sim_seq
from wub.simulate import genome as sim_genome
from wub.util import parse as parse_util
from wub.util import seq as seq_util
from wub.bam import sam_writer

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Sample fragments from the input genome and simulate sequencing errors.
    Read lengths are drawn from the specified truncated gamma distribution. Chromosomes are sampled randomly
    for each read.

    The format of the read names is the following:
    r<unique_id>_<chromosome>_<frag_start>_<frag_end>_<strand>/q<realised_quality>/s<realised_substiutions>/d<realised_deletions>/i<realised_insertions>
    """)
parser.add_argument(
    '-n', metavar='nr_reads', type=int, help="Number of simulated reads (1).", default=1)
parser.add_argument('-m', metavar='mean_length', type=int,
                    help="Mean read length (5000).", default=5000)
parser.add_argument(
    '-a', metavar='gamma_shape', type=float, help="Read length distribution: gamma shape parameter (1).", default=1.0)
parser.add_argument(
    '-l', metavar='low_trunc', type=int, help="Read length distribution: lower truncation point (100).", default=100)
parser.add_argument(
    '-u', metavar='high_trunc', type=int, help="Read length distribution: upper truncation point (None).", default=None)
parser.add_argument(
    '-e', metavar='error_rate', type=float,
    help="Total rate of substitutions insertions and deletions (0.1).", default=0.1)
parser.add_argument('-w', metavar='error_weights', type=str,
                    help="Relative frequency of substitutions,insertions,deletions (1,1,4).", default="1,1,4")
parser.add_argument(
    '-b', metavar='strand_bias', type=float, help="Strand bias: the ratio of forward and reverse reads (0.5).", default=0.5)
parser.add_argument(
    '-q', metavar='mock_quality', type=int, help="Mock base quality for fastq output (40).", default=40)
parser.add_argument(
    '-s', metavar='true_sam', type=str, help="Save true alignments in this SAM file (None).", default=None)
parser.add_argument('-z', metavar='random_seed', type=int,
                    help="Random seed (None).", default=None)
parser.add_argument('input_fasta', nargs='?', help='Input genome in fasta format (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_fastq', nargs='?', help='Output fastq (default: stdout)',
                    type=argparse.FileType('w'), default=sys.stdout)

# Possible extension to sequencing error simulation:
# - Biased error probabilities and biased insert composition.
# - Suggestion by cwright: simulate non-uniformity of errors.
# - Simulate insert length from distribution.
# - Save simulation details to CSV file.
# - Simulate deletions longer than one? Does it make sense? Might be challenging.


def build_sam_header(chromosomes):
    """Build SAM header from chromosomes.

    :param chromosomes: List of chromosomes as SeqRecord obejcts.
    :returns: Header information.
    :rtype: OrderedDict
    """
    header = OrderedDict()
    header['HD'] = [OrderedDict([('VN', '1.5'), ('SO', 'unsorted')])]
    chrom_info = []
    for chrom in chromosomes:
        chrom_info.append(OrderedDict([('SN', chrom.id), ('LN', len(chrom))]))
    header['SQ'] = chrom_info
    return header


def simulate_sequencing(chromosomes, mean_length, gamma_shape, low_truncation,
                        high_truncation, error_rate, error_weights, strand_bias, mock_quality, number_reads, sam_writer=False):
    """Simulate sequencing.

    :param chromosomes: Input chromosomes (list of SeqRecord obejcts).
    :param mean_length: Mean read length.
    :param gamma_shape: Shape parameter of the read length distribution.
    :param low_truncation: Minimum read length.
    :param high_truncation: Maximum read length.
    :param error_rate: Total error rate.
    :param error_weights: "Relative frequency of substitutions,insertions,deletions.
    :param strand_bias: Strand bias: the ratio of forward and reverse reads.
    :param mock_quality: Mock base quality for fastq output.
    :param number_reads: Number of reads to simulate.
    :param sam_writer: SAM writer object.
    """
    for fragment in sim_genome.simulate_fragments(chromosomes, mean_length, gamma_shape,
                                                  low_truncation, high_truncation, number_reads):
        frag_seq = fragment.seq
        # Sample read direction:
        direction = sim_seq.sample_direction(strand_bias)

        # Simulate sequencing errors:
        mutated_record = sim_seq.simulate_sequencing_errors(
            frag_seq, error_rate, error_weights)

        # Special case:
        if len(mutated_record.seq) == 0:
            sys.stderr.write(
                "The length of fragment {} is 0 after simulating sequencing errors!".format(fragment.uid))

        # Construct read name:
        read_name = "r{}_{}_{}_{}_{}".format(
            fragment.uid, fragment.chrom, fragment.start, fragment.end, direction)
        read_name = "{}/q{}/s{}/d{}/i{}".format(read_name, mutated_record.real_qual,
                                                mutated_record.real_subst, mutated_record.real_del,
                                                mutated_record.real_ins)

        mock_qualities = [mock_quality] * len(mutated_record.seq)

        sam = None
        if sam_writer:
            # Construct SAM flag:
            flag = 0 | 0x2
            if direction == '-':
                flag = flag | 0x10

            # Construct SAM record:
            sam = sam_writer.new_sam_record(qname=read_name, flag=flag, rname=fragment.chrom,
                                            pos=fragment.start + 1, mapq=93, cigar=mutated_record.cigar,
                                            rnext='*', pnext=0, tlen=0, seq=mutated_record.seq, qual=pysam.qualities_to_qualitystring(mock_qualities))

        read_seq = mutated_record.seq
        if direction == '-':
            read_seq = seq_util.reverse_complement(mutated_record.seq)

        yield seq_util.new_dna_record(read_seq, read_name, mock_qualities), sam

if __name__ == '__main__':
    args = parser.parse_args()

    # Set random seed:
    if args.z is not None:
        np.random.seed(args.z)

    # Read in chromosomes of the input genome:
    chromosomes = list(seq_util.read_seq_records(args.input_fasta))

    # Process error weights:
    error_weights = np.array(parse_util.separated_list_to_floats(args.w))
    # Normalise error weights to probabilities:
    error_weights = parse_util.normalise_array(error_weights)
    error_weights = dict(
        zip(['substitution', 'insertion', 'deletion'], error_weights))

    sw = None
    if args.s is not None:
        sw = sam_writer.SamWriter(args.s, build_sam_header(chromosomes))

    simulation_iterator = simulate_sequencing(
        chromosomes, args.m, args.a, args.l, args.u, args.e, error_weights, args.b, args.q, args.n, sw)

    for simmed, sam in simulation_iterator:
        seq_util.write_seq_records(
            simmed, args.output_fastq, format='fastq')
        if sw is not None:
            sw.write(sam)

    if sw is not None:
        sw.close()
