#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import sys
import pandas as pd
from collections import OrderedDict
import tqdm

from wub.bam import common as bam_common

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Produce a tab separated file of alignment lengths and other information.
    Rows are sorted by number of aligned reference bases unless the -x option is specified.
    """)
parser.add_argument(
    '-t', metavar='tsv_file', type=str, help="Tab separated file to save alignment lengths (bam_alignment_length.tsv).", required=False, default="bam_alignment_length.tsv")
parser.add_argument(
    '-q', metavar='aqual', type=int, default=0, help="Minimum alignment quality (0).")
parser.add_argument(
    '-x', action="store_true", help="Sort by number of read bases instead of number of aligned reference bases.", default=False)
parser.add_argument(
    '-Q', action="store_true", help="Be quiet and do not print progress bar (False).", default=False)
parser.add_argument(
    'bam', metavar='bam', type=str, help="Input BAM file.")


if __name__ == '__main__':
    args = parser.parse_args()
    verbose = not args.Q

    bam_reader = bam_common.pysam_open(args.bam, in_format='BAM')

    if verbose:
        sys.stdout.write(
            "Gathering read and alignment lengths from file: {}\n".format(args.bam))
        try:
            total_reads = bam_reader.mapped + bam_reader.unmapped
        except:
            total_reads = None
        bam_reader = tqdm.tqdm(bam_reader, total=total_reads)

    read_names = []
    ref_names = []
    ref_lengths = []
    read_lengths = []
    aln_lengths = []
    mapping_quals = []

    # Gather alignment information:
    for record in bam_reader:
        if (not record.is_unmapped) and (record.mapping_quality > args.q):
            read_names.append(record.query_name)
            ref_names.append(record.reference_name)
            read_lengths.append(len(record.query_sequence))
            aln_lengths.append(record.query_alignment_length)
            ref_lengths.append(record.reference_length)
            mapping_quals.append(record.mapping_quality)

    # Construct data frame:
    data = OrderedDict([('read_name', read_names),
                        ('aligned_ref_bases', ref_lengths),
                        ('aligned_read_bases', aln_lengths),
                        ('read_length', read_lengths),
                        ('reference', ref_names),
                        ('mapping_quality', mapping_quals)
                        ])

    df = pd.DataFrame(data)
    del data, read_names, ref_names, mapping_quals
    del read_lengths, aln_lengths, ref_lengths

    # Sort data frame and save tsv:
    sort_by = 'aligned_ref_bases'
    if args.x:
        sort_by = 'aligned_read_bases'

    df.sort_values([sort_by], ascending=[0], inplace=True)
    df.to_csv(args.t, sep="\t", index=False)
