#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import os
import sys
import pandas as pd
from collections import OrderedDict
from scipy.stats import gaussian_kde
from scipy.optimize import minimize_scalar

from wub.util import misc
from wub.vis import report
from wub.bam import stats

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Produce accuracy statistics of the input BAM file. Calculates global accuracy and identity and various per-read statistics.
    The input BAM file must be sorted by coordinates and indexed.
    """)
parser.add_argument(
    '-c', metavar='region', type=str, help="BAM region (None).", required=False, default=None)
parser.add_argument(
    '-g', metavar='global_tsv', type=str, default=None, help="Tab separated file to save global statistics (None).", required=False)
parser.add_argument(
    '-l', metavar='read_tsv', type=str, default=None, help="Tab separated file to save per-read statistics (None).", required=False)
parser.add_argument(
    '-t', metavar='bam_tag', type=str, default=None, help="Dataset tag (BAM basename).", required=False)
parser.add_argument(
    '-q', metavar='aqual', type=int, default=0, help="Minimum alignment quality (0).")
parser.add_argument(
    '-e', action="store_true", default=False, help="Include hard and soft clipps in alignment length when calculating accuracy (False).")
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF (bam_accuracy.pdf).", default="bam_accuracy.pdf")
parser.add_argument(
    '-p', metavar='results_pickle', type=str, help="Save pickled results in this file (None).", default=None)
parser.add_argument(
    '-Q', action="store_true", help="Be quiet and do not print progress bar (False).", default=False)
parser.add_argument(
    'bam', metavar='bam', type=str, help="Input BAM file.")


def estimate_mode(acc):
    """ Estimate the mode of a set of float values between 0 and 1.

    :param acc: Data.
    :returns: The mode of the sample
    :rtype: float
    """
    # Taken from sloika.
    if len(acc) > 1:
        da = gaussian_kde(acc)
        optimization_result = minimize_scalar(lambda x: -da(x), bounds=(0, 1), method='brent')
        if optimization_result.success:
            try:
                mode = optimization_result.x[0]
            except IndexError:
                mode = optimization_result.x
        else:
            sys.stderr.write("Mode computation failed")
            mode = 0
    else:
        mode = acc[0]
    return mode


def base_stats_qc(st, report):
    """ Plot base statistics.

    :param st: Statistics dict.
    :param report: Plotter object.
    :returns: None
    """

    bs = st.copy()
    del bs['accuracy']
    del bs['identity']
    plotter.plot_bars_simple(
        bs, title="Basewise statistics", xlab="Type", ylab="Count")
    plotter.plot_bars_simple(OrderedDict([('Identity ({})'.format(st['identity']), st['identity']), ('Accuracy ({})'.format(
        st['accuracy']), st['accuracy'])]), title="Precision statistics: length weighted means", xlab="Type", ylab="Count")


def read_precision_qc(st, report):
    """ Plot read precision statistics.

    :param st: Statistics dict.
    :param report: Plotter object.
    :returns: Mode of accuracy and identity.
    :rtype: dict
    """
    accuracy_mode = estimate_mode(st['accuracy'])
    report.plot_histograms(OrderedDict([('Accuracy', st[
        'accuracy'])]), title="Distribution of per-read accuracies", xlab="Accuracy", ylab="Count", legend=True,
        vlines={'Mode:{0:.4f}'.format(accuracy_mode): accuracy_mode})

    identity_mode = estimate_mode(st['identity'])
    report.plot_histograms(OrderedDict([('Identity', st[
        'identity'])]), title="Distribution of per-read identitities", xlab="Identity", ylab="Count", legend=True,
        vlines={'Mode:{0:.4f}'.format(identity_mode): identity_mode})

    modes = {'accuracy_mode': accuracy_mode, 'identity_mode': identity_mode}
    return modes


if __name__ == '__main__':
    args = parser.parse_args()
    tag = args.t if args.t is not None else os.path.basename(args.bam)

    plotter = report.Report(args.r)

    read_stats = stats.read_stats(
        args.bam, region=args.c, min_aqual=args.q, with_clipps=args.e, verbose=not args.Q)
    read_stats['tag'] = tag
    base_stats = read_stats['base_stats']
    precision_stats = read_stats['read_stats']

    base_stats_qc(base_stats, plotter)
    modes = read_precision_qc(precision_stats, plotter)

    plotter.close()

    global_stats = OrderedDict([
        ('Accuracy', [read_stats['base_stats']['accuracy']]),
        ('AccuracyMode', modes['accuracy_mode']),
        ('Identity', [read_stats['base_stats']['identity']]),
        ('IdentityMode', modes['identity_mode']),
        ('Mapped', [read_stats['mapped']]),
        ('Unmapped', [read_stats['unmapped']]),
        ('Tag', [read_stats['tag']]), ])
    global_stats = pd.DataFrame(global_stats)

    if args.g is not None:
        global_stats.to_csv(args.g, sep="\t", index=False)

    if args.l is not None:
        read_df = pd.DataFrame(precision_stats)
        read_df.to_csv(args.l, sep="\t", index=False)

    if args.p is not None:
        misc.pickle_dump(read_stats, args.p)
