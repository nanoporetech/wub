#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import numpy as np
from os import path
from wub.vis import report
from wub.util import seq as seq_util

import warnings
warnings.filterwarnings('ignore')

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Plot histograms of length distributions from multiple sequence files.""")
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF.", default="multi_length_hist.pdf")
parser.add_argument(
    '-f', metavar='in_format', type=str, help="Input format (fastq).", default="fastq")
parser.add_argument(
    '-b', metavar='nr_bins', type=int, help="Number of bins (50).", default=50)
parser.add_argument(
    '-l', metavar='min_len', type=int, help="Minimum read length (None).", default=None)
parser.add_argument(
    '-u', metavar='max_len', type=int, help="Maximum read length (None).", default=None)
parser.add_argument(
    '-L', action="store_true", help="Log transform lengths.", default=False)
parser.add_argument(
    'in_files', metavar='input_counts', nargs='*', type=str, help="Input sequence files.")


def _get_lengths(in_file, in_format, min_length, max_length, do_log):
    """ Iterate over input and accumulate sequence lengths. """
    input_iterator = seq_util.read_seq_records(in_file, format=in_format)
    lengths = []
    for record in input_iterator:
        length = len(record)
        # Filter for minimum read length:
        if (min_length is not None) and (length < min_length):
            continue
        # Filter for maximum read length:
        if (max_length is not None) and (length > max_length):
            continue
        if do_log:
            length = np.log(length)
        lengths.append(length)
    input_iterator.close()
    return lengths


if __name__ == '__main__':
    args = parser.parse_args()
    plotter = report.Report(args.r)

    if len(args.in_files) == 0:
        sys.stderr.write("No input files given!\n")
        sys.exit(1)

    data_map = {}
    for in_file in args.in_files:
        name = path.basename(in_file).rsplit('.', 1)[0]
        data_map[name] = _get_lengths(in_file, args.f, args.l, args.u, args.L)

    if args.L:
        xlab = 'log(read length)'
    else:
        xlab = 'read length'

    plotter.plot_histograms(data_map, title='Read length distributions', xlab=xlab, ylab='Count', bins=args.b, alpha=0.7, legend_loc='best', legend=True, vlines=None)

    plotter.close()
