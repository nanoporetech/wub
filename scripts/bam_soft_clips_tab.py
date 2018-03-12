#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

import pandas as pd
import pysam
import tqdm
from collections import OrderedDict

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Produce a tab separated file with read identifiers and number of soft clipped bases at each end (relative to the original sequence in the fastq).""")
parser.add_argument(
    '-t', metavar='tsv', type=str, default="bam_soft_clips_tab.tsv", help="Output tab separated file.", required=False)
parser.add_argument(
    '-Q', action="store_true", help="Be quiet and do not print progress bar (False).", default=False)
parser.add_argument(
    'bam', metavar='bam', type=str, help="Input BAM file.")


def _get_clips(cigar, is_reverse):
    """ Get clips at the start and end relative to the original sequence. """
    clip_start, clip_end = 0, 0

    # Consider the first CIGAR tuple:
    if cigar[0][0] == 4:
        clip_start = cigar[0][1]

    # Consider the last CIGAR tuple:
    if cigar[-1][0] == 4:
        clip_end = cigar[-1][1]

    # Reverse orientation if necessary:
    if is_reverse:
        clip_start, clip_end = clip_end, clip_start
    return clip_start, clip_end


def process_reads(alignment_file, in_format='BAM', verbose=False):
    """Process reads and extract the corresponding information.

    :param alignment_file: BAM file.
    :param verbose: Verbosity flag.
    :returns: pandas dataframe with reads and soft clip lengths.
    :rtype: pandas.DataFrame
    """
    reads, strand, clip_start, clip_end = [], [], [], []
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
        reads.append(segment.query_name)
        strand.append('-' if segment.is_reverse else '+')
        cs, ce = _get_clips(segment.cigartuples, segment.is_reverse)
        clip_start.append(cs)
        clip_end.append(ce)

    aln_iter.close()

    data = OrderedDict([('Read', reads), ('Strand', strand), ('ClipStart', clip_start), ('ClipEnd', clip_end)])
    df = pd.DataFrame(data)

    return df


if __name__ == '__main__':
    args = parser.parse_args()
    verbose = not args.Q

    df = process_reads(args.bam, verbose=verbose)
    df.to_csv(args.t, sep="\t", index=False)
