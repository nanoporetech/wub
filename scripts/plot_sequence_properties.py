#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

from wub.util import seq as seq_util
from wub.vis import report

import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import seaborn as sns
warnings.resetwarnings()

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Plot histograms of lengths and quality values.')
parser.add_argument(
    '-f', metavar='format', type=str, help="Input format (fastq).", default='fastq')
parser.add_argument(
    '-b', metavar='bins', type=int, help="Number of bins on histograms (50).", default=50)
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report pdf (plot_sequence_properties.pdf).", default='plot_sequence_properties.pdf')
parser.add_argument(
    '-j', help="Produce joint plot of lengths and mean quality values (False).", default=False, action="store_true")
parser.add_argument('input_fastx', nargs='?', help='Input (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)


if __name__ == '__main__':
    args = parser.parse_args()

    in_format = args.f
    input_iterator = seq_util.read_seq_records(
        args.input_fastx, format=in_format)

    # Could be more efficient with dictionaries if we did not have to
    # deal with the joint plot.
    lengths = []
    mean_qualities = []

    for record in input_iterator:
        lengths.append(len(record))
        if in_format == 'fastq':
            mean_quality = seq_util.mean_qscore(record.letter_annotations["phred_quality"])
            mean_qualities.append(mean_quality)

    plotter = report.Report(args.r)

    plotter.plot_histograms(
        {'lengths': lengths}, title="Distribution of sequence lengths", xlab="Length", ylab="Count", legend=False)

    if in_format == 'fastq':
        plotter.plot_histograms(
            {'qualities': mean_qualities}, title="Distribution of mean base qualities", xlab="Mean base quality", ylab="Count", legend=False)
        if args.j:
            plotter.plot_arrays({'scatter': (lengths, mean_qualities)}, title="Sequence length vs. mean base quality",
                                xlab="Sequence length", ylab="Mean base quality", legend=False)

    plotter.close()
