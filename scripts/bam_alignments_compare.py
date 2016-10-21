#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import pandas as pd
from collections import OrderedDict

from wub.util import misc
from wub.bam import compare as bam_compare

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Compare alignments stored in two BAM files.
    The two BAM files must have the same set of reads in the same order (name sorted).
    """)
parser.add_argument(
    '-a', metavar='min_aqual', type=int, help="Minimum mapping quality (0).", default=0)
parser.add_argument(
    '-p', metavar='results_pickle', type=str, help="Save pickled results in this file (None).", default=None)
parser.add_argument(
    '-t', metavar='tsv_file', type=str, help="Save results in tsv format in this file (bam_count_reads.tsv).", default="bam_count_reads.tsv")
parser.add_argument(
    '-s', action="store_true", help="Sort results before saving in tsv format (False).", default=False)
parser.add_argument(
    'bam_one', metavar='bam_one', type=str, help="First input BAM file.")
parser.add_argument(
    'bam_two', metavar='bam_two', type=str, help="Second input BAM file.")

if __name__ == '__main__':
    args = parser.parse_args()

    stats = bam_compare.bam_compare(args.bam_one, args.bam_two, in_format='BAM')

    if args.t is not None:
        misc.pickle_dump(dict(stats), args.t)
