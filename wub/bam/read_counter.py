# -*- coding: utf-8 -*-
"""Count reads per reference in BAM/SAM file."""
import sys

import pysam
from collections import defaultdict
import tqdm


def count_reads(alignment_file, in_format='BAM', min_aln_qual=0, verbose=False):
    """Count reads mapping to references in a BAM file.

    :param alignment_file: BAM file.
    :param min_aln_qual: Minimum mapping quality.
    :param verbose: Minimum mapping quality.
    :returns: Dictionary with read counts per reference.
    :rtype: dict
    """
    counts = defaultdict(int)
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
        if segment.mapping_quality >= min_aln_qual:
            counts[segment.reference_name] += 1

    aln_iter.close()

    return dict(counts)


def count_reads_realtime(alignment_file='-', in_format='SAM', min_aln_qual=0, yield_freq=1, verbose=False):
    """Online counting of reads mapping to references in a SAM/BAM stream from stdin.

    :param alignment_file: BAM file (stdin).
    :param min_aln_qual: Minimum mapping quality.
    :param yield_freq: Yield frequency.
    :param verbose: Minimum mapping quality.
    :returns: Generator of dictionary with read counts per reference.
    :rtype: generator
    """
    counts = defaultdict(int)
    if in_format == 'BAM':
        mode = "rb"
    elif in_format == 'SAM':
        mode = "r"
    else:
        raise Exception("Invalid format: {}".format(in_format))

    aln_iter = pysam.AlignmentFile(alignment_file, mode)

    if verbose:
        sys.stdout.write(
            "Online counting of read statistics from file: {}\n".format(alignment_file))
        aln_iter = iter(tqdm.tqdm(aln_iter))

    nr_mapped = 0
    while True:
        try:
            segment = aln_iter.next()
        except StopIteration:
            # Final yield:
            yield counts
            return

        if segment.is_unmapped:
            continue
        if segment.mapping_quality >= min_aln_qual:
            counts[segment.reference_name] += 1
        nr_mapped += 1

        if nr_mapped % yield_freq == 0:
            yield counts

    aln_iter.close()
