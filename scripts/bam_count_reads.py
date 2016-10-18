#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse

import pandas as pd
from collections import OrderedDict

from wub.bam import read_counter
from wub.util import misc

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Count reads mapping to each reference in a BAM file.""")
parser.add_argument(
    '-a', metavar='min_aqual', type=int, help="Minimum mapping quality (0).", default=0)
parser.add_argument(
    '-p', metavar='results_pickle', type=str, help="Save pickled results in this file (None).", default=None)
parser.add_argument(
    '-t', metavar='tsv_file', type=str, help="Save results in tsv format in this file (bam_count_reads.tsv).", default="bam_count_reads.tsv")
parser.add_argument(
    'bam', metavar='bam_file', type=str, help="Input BAM file.")

if __name__ == '__main__':
    args = parser.parse_args()

    counts = read_counter.count_reads(args.bam, args.a)
    data = OrderedDict([('Reference', counts.keys()), ('Count', counts.values())])
    data_frame = pd.DataFrame(data)

    if args.t is not None:
        data_frame.to_csv(args.t, sep='\t')

    if args.p is not None:
        misc.pickle_dump(data, args.p)
