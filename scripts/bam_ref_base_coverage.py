#!/usr/bin/env python
# -*- coding: utf-8 -*-

import six
import argparse
import tqdm
import os
import pandas as pd
from Bio import SeqIO
from wub.bam import stats as bam_stats

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Calculate percent covered reference lengths.""")
parser.add_argument(
    '-f', metavar='reference', type=str, help="Reference fasta.", required=True)
parser.add_argument(
    '-c', metavar='region', type=str, help="BAM region (None).", required=False, default=None)
parser.add_argument(
    '-t', metavar='tsv', type=str, default="bam_ref_base_coverage.tsv", help="Output tab separated file (bam_ref_base_coverage.tsv).", required=False)
parser.add_argument(
    '-m', metavar='min_cov', type=int, default=1, help="Minimum base coverage for a position to be counted (1).")
parser.add_argument(
    '-Q', action="store_true", help="Be quiet and do not show progress bars.", default=False)
parser.add_argument(
    'bam', metavar='bam', type=str, help="Input BAM file.")


if __name__ == '__main__':
    args = parser.parse_args()
    verbose = not args.Q
    tag = args.t
    if tag is None:
        tag = os.path.basename(args.bam)

    # Load reference lengths:
    references = SeqIO.index(args.f, format='fasta')
    chrom_lengths = {name: len(so) for name, so in six.iteritems(references)}
    # Parse fragments:
    st = bam_stats.pileup_stats(args.bam, region=args.c, verbose=verbose, with_quals=False)['coverage']

    res = {}
    for chrom, chrom_length in six.iteritems(chrom_lengths):
        # No coverage:
        if chrom not in st:
            res[chrom] = 0.0
        else:
            nr_hits = 0
            # Iterate over covered positions and count valid hits:
            for pos, cov in six.iteritems(st[chrom]):
                if cov >= args.m:
                    nr_hits += 1
            # Calculate percent covered reference length:
            res[chrom] = float(nr_hits * 100) / chrom_length

    # Convert results to sorted data frame:
    df = pd.DataFrame({'Chrom': list(res.keys()), 'Percent_cov': list(res.values())})
    df.sort_values(['Percent_cov'], ascending=[0], inplace=True)
    df.to_csv(args.t, sep="\t", index=False)
