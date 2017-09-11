#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import argparse
import tqdm

import os
import sys
import pandas as pd

import pysam
from wub.vis import report
from wub.util import seq as seq_util

import warnings
warnings.simplefilter("ignore")
import seaborn as sns

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""
    Produce a plot of GC content of aligned read and reference portion versus their mean quality values.
    """)
parser.add_argument('-f', metavar='reference', type=str, help="Reference fasta.", required=True)
parser.add_argument(
    '-q', metavar='aqual', type=int, default=0, help="Minimum alignment quality (0).")
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF (bam_gc_vs_qual.pdf).", default="bam_gc_vs_qual.pdf")
parser.add_argument(
    '-Q', action="store_true", help="Be quiet and do not show progress bars.", default=False)
parser.add_argument(
    'bam', metavar='bam', type=str, help="Input BAM file.")


def _process_reads(alignment_file, refs, in_format='BAM', min_aln_qual=0, verbose=False):
    """
    Gather information about the GC content and mean quality value of aligned portions of the reads.
    """
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
            "Gathering GC content vs. quality information from file: {}\n".format(alignment_file))
        if in_format == "BAM":
            aln_iter = tqdm.tqdm(aln_iter, total=total_reads)

    rgcs, gcs, quals = [], [], []
    for segment in aln_iter:
        if segment.is_unmapped:
            continue
        if segment.mapping_quality >= min_aln_qual:
            # Calculate GC content of aligned read portion:
            aln_seq = segment.query_alignment_sequence
            gcs.append(seq_util.gc_content(aln_seq))

            # Calculate GC content of aligned reference:
            ref_seq = refs[segment.reference_name].seq[segment.reference_start:segment.reference_end]
            rgcs.append(seq_util.gc_content(ref_seq))

            # Calculate mean quality score of aligned read portion:
            aln_quals = segment.query_alignment_qualities
            quals.append(seq_util.mean_qscore(aln_quals, qround=False))

    aln_iter.close()

    df = pd.DataFrame({'GC_content': gcs, 'MeanQuality': quals, 'GC_content_ref': rgcs})

    return df


if __name__ == '__main__':
    args = parser.parse_args()
    verbose = not args.Q
    tag = os.path.basename(args.bam)

    references = seq_util.read_seq_records_dict(args.f)
    data = _process_reads(args.bam, references, min_aln_qual=args.q, verbose=verbose)

    # Plot GC content of aligned read portion vs. mean quality.
    plotter = report.Report(args.r)
    sns.jointplot("GC_content", "MeanQuality", kind="hex", data=data)
    plotter.plt.tight_layout()
    plotter.pages.savefig()
    plotter.plt.clf()

    # Plot GC content of aligned reference portion vs. mean quality.
    sns.jointplot("GC_content_ref", "MeanQuality", kind="hex", data=data)
    plotter.plt.tight_layout()
    plotter.pages.savefig()
    plotter.plt.clf()

    plotter.close()
