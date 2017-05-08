#!/usr/bin/env python
# -*- coding: utf-8 -*-

import six
import argparse
import numpy as np
from collections import OrderedDict
import pandas as pd
from os import path
from wub.vis import report
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import seaborn as sns
warnings.resetwarnings()
_ = sns

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Correlate counts produced by multiple runs of bam_count_reads.py.""")
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF (bam_multi_qc.pdf).", default="correlate_counts.pdf")
parser.add_argument(
    '-x', action="store_true", help="Produce pairwise jointplots.", default=False)
parser.add_argument(
    'counts', metavar='input_counts', nargs='*', type=str, help="Input counts as tab separated files.")


def load_counts(counts):
    """Load statistics from pickle files.

    :param pickles: List of pickle files.
    :returns: OrderedDict of stats per dataset.
    :rtype: OrderedDict
    """
    stats = OrderedDict()
    for count_file in counts:
        stats[name] = pd.read_csv(count_file, sep="\t")
    return stats


if __name__ == '__main__':
    args = parser.parse_args()
    plotter = report.Report(args.r)

    stats = load_counts(args.counts)

    plotter.close()
