#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

import pandas as pd
from collections import OrderedDict

from wub.bam import read_counter
from wub.util import misc
from wub.util import seq as seq_util

import tqdm

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Count reads mapping to each reference in a BAM file.""")
parser.add_argument(
    '-a', metavar='min_aqual', type=int, help="Minimum mapping quality (0).", default=0)
parser.add_argument(
    '-f', metavar='in_format', type=str, help="Input format (BAM).", default='BAM')
parser.add_argument(
    '-z', metavar='ref_fasta', type=str, help="Reference fasta. GC content and length columns are added if present (None).", default=None)
parser.add_argument(
    '-p', metavar='results_pickle', type=str, help="Save pickled results in this file (None).", default=None)
parser.add_argument(
    '-t', metavar='tsv_file', type=str, help="Save results in tsv format in this file (bam_count_reads.tsv).", default="bam_count_reads.tsv")
parser.add_argument(
    '-s', action="store_true", help="Sort results before saving in tsv format (False).", default=False)
parser.add_argument(
    '-Q', action="store_true", help="Be quiet and do not print progress bar (False).", default=False)
parser.add_argument(
    'bam', metavar='bam_file', type=str, help="Input BAM file.")

if __name__ == '__main__':
    args = parser.parse_args()

    counts = read_counter.count_reads(
        args.bam, in_format=args.f, min_aln_qual=args.a, verbose=not args.Q)

    if args.s:
        counts = OrderedDict(
            sorted((item for item in counts.iteritems()), key=lambda x: x[1], reverse=True))

    data = OrderedDict()

    if args.z is not None:
        lengths, gc_contents = {}, {}
        ref_iter = seq_util.read_seq_records(args.z)
        if not args.Q:
            sys.stderr.write("Calculating sequence features:\n")
            ref_iter = tqdm.tqdm(ref_iter)

        for ref in ref_iter:
            # Augment counts dictionary with missing reference entries:
            if ref.id not in counts:
                counts[ref.id] = 0
            lengths[ref.id] = len(ref)
            gc_contents[ref.id] = seq_util.gc_content(str(ref.seq))
        data['Length'] = [lengths[tr] for tr in counts.iterkeys()]
        data['GC'] = [gc_contents[tr] for tr in counts.iterkeys()]

    data['Reference'] = counts.keys()
    data['Count'] = counts.values()
    data_frame = pd.DataFrame(data)

    if args.t is not None:
        data_frame.to_csv(args.t, sep='\t', index=False)

    if args.p is not None:
        misc.pickle_dump(data, args.p)
