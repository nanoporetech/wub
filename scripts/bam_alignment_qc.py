#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import os
import pandas as pd
from collections import OrderedDict

from wub.util import misc
from wub.vis import report
from wub.bam import stats
from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Produce alignment based QC plots of the input BAM file.
    The input BAM file must be sorted by coordinates and indexed.
    """)
parser.add_argument(
    '-f', metavar='reference', type=str, help="Reference fasta.", required=True)
parser.add_argument(
    '-c', metavar='region', type=str, help="BAM region (None).", required=False, default=None)
parser.add_argument(
    '-t', metavar='bam_tag', type=str, default=None, help="Dataset tag (BAM basename).", required=False)
parser.add_argument(
    '-q', metavar='aqual', type=int, default=5, help="Minimum alignment quality (5).")
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF (bam_alignments_compare.pdf).", default="bam_alignments_compare.pdf")
parser.add_argument(
    '-p', metavar='results_pickle', type=str, help="Save pickled results in this file (bam_alignment_qc.pk).", default="bam_alignment_qc.pk")
parser.add_argument(
    'bam', metavar='bam', type=str, help="Input BAM file.")

if __name__ == '__main__':
    args = parser.parse_args()

    # print stats.read_stats(args.bam, args.c)
    # print stats.pileup_stats(args.bam, args.c)
    references = seq_util.read_seq_records_dict(args.f)
    print stats.error_and_read_stats(args.bam, references, region=args.c)
