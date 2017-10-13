#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import numpy as np

from wub.util import seq as seq_util
from wub.vis import report

import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import seaborn as sns
warnings.resetwarnings()
_ = sns

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Plot the mean quality values across non-overlapping windows in the input sequences.')
parser.add_argument(
    '-w', metavar='win_size', type=int, help="Window size (50).", default=50)
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report pdf (plot_qualities.pdf).", default='plot_qualities.pdf')
parser.add_argument('input_fastx', nargs='?', help='Input (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)


def _smooth_qualitites(quals, winsize):
    """ Smooth out qualities by taking average of non-overlapping windows. """
    smooth_quals = []
    for i in range(0, len(quals) - winsize, winsize):
        smooth_quals.append(np.mean(quals[i:i + winsize]))
    smooth_quals = np.array(smooth_quals, dtype=float)
    return smooth_quals


if __name__ == '__main__':
    args = parser.parse_args()

    input_iterator = seq_util.read_seq_records(
        args.input_fastx, format="fastq")

    plotter = report.Report(args.r)

    for record in input_iterator:
        quals = np.array(record.letter_annotations["phred_quality"])
        smooth_quals = _smooth_qualitites(quals, args.w)
        pos = np.arange(len(smooth_quals))
        data_map = {'Mean qualities': (pos, smooth_quals)}
        plotter.plot_arrays(data_map, marker='-', title=record.id, xlab="Window", ylab="Mean quality")

    plotter.close()
