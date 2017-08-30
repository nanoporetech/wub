#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import pandas as pd
import pysam
import tqdm
from collections import OrderedDict

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Produce a tab separated file with read identifiers and the corresponding references, sorted by reference.""")
parser.add_argument(
    '-t', metavar='read_tsv', type=str, default="bam_ref_tab.tsv", help="Tab separated file to save reference table.", required=False)
parser.add_argument(
    '-Q', action="store_true", help="Be quiet and do not print progress bar (False).", default=False)
parser.add_argument(
    'bam', metavar='bam', type=str, help="Input BAM file.")


def process_reads(alignment_file, in_format='BAM', verbose=False):
    """Process reads and extract the corresponding reference.

    :param alignment_file: BAM file.
    :param verbose: Verbosity flag.
    :returns: pandas dataframe with reads and references.
    :rtype: dict
    """
    reads, refs = [], []
    if in_format == 'BAM':
        mode = "rb"
    elif in_format == 'SAM':
        mode = "r"
    else:
        raise Exception("Invalid format: {}".format(in_format))

    aln_iter = pysam.AlignmentFile(alignment_file, mode)

    if verbose and in_format == "BAM":
        try:
            total_reads = aln_iter.mapped + aln_iter.unmapped
        except:
            total_reads = None
        sys.stdout.write(
            "Gathering read statistics from file: {}\n".format(alignment_file))
        if in_format == "BAM":
            aln_iter = tqdm.tqdm(aln_iter, total=total_reads)

    for segment in aln_iter:
        if segment.is_unmapped:
            continue
        refs.append(segment.reference_name)
        reads.append(segment.query_name)

    aln_iter.close()

    data = OrderedDict([('Read', reads), ('Reference', refs)])
    df = pd.DataFrame(data)

    return df


if __name__ == '__main__':
    args = parser.parse_args()
    verbose = not args.Q

    df = process_reads(args.bam, verbose=verbose)
    df.sort_values(['Reference'], ascending=[0], inplace=True)
    df.to_csv(args.t, sep="\t", index=False)
