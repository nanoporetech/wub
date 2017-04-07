#!/usr/bin/env python
# -*- coding: utf-8 -*-

import six
import argparse
import os
import pandas as pd
from collections import OrderedDict
from wub.util import misc
from wub.bam import compare as bam_compare
from wub.vis import report

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Compare alignments stored in two BAM files.
    The two BAM files must have the same set of reads in the same order (name sorted).
    """)
parser.add_argument(
    '-w', metavar='coarse_tolerance', type=int, help="Tolerance when performing coarse comparison of alignments (50).", default=50)
parser.add_argument(
    '-g', action="store_true", help="Do strict comparison of alignment flags.", default=False)
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF (bam_alignments_compare.pdf).", default="bam_alignments_compare.pdf")
parser.add_argument(
    '-p', metavar='results_pickle', type=str, help="Save pickled results in this file (bam_alignments_compare.pk).", default="bam_alignments_compare.pk")
parser.add_argument(
    '-t', metavar='tsv_file', type=str, help="Save results in tsv format in this file (None).", default=None)
parser.add_argument(
    '-f', metavar='format', type=str, help="Input format (BAM).", default='BAM')
parser.add_argument(
    '-Q', action="store_true", help="Be quiet and do not print progress bar (False).", default=False)
parser.add_argument(
    'bam_one', metavar='bam_one', type=str, help="First input BAM file.")
parser.add_argument(
    'bam_two', metavar='bam_two', type=str, help="Second input BAM file.")

if __name__ == '__main__':
    args = parser.parse_args()

    stats = bam_compare.bam_compare(
        args.bam_one, args.bam_two, coarse_tolerance=args.w, strict_flags=args.g, in_format=args.f, verbose=not args.Q)

    plotter = report.Report(args.r)

    query_stats = OrderedDict((field, stats[field]) for field in (
        'TotalQueries', 'CoarseMatches', 'DirectionMismatch', 'SeqMismatch'))
    plotter.plot_bars_simple(
        query_stats, title="Per query statistics", xlab="Field", ylab="Count", auto_limit=False)

    query_stats = OrderedDict((field, stats[field]) for field in (
        'TotalQueries', 'StrictFlagMismatch', 'RefMismatch'))
    plotter.plot_bars_simple(
        query_stats, title="Per query statistics (continued)", xlab="Field", ylab="Count", auto_limit=False)

    aligned_bases = OrderedDict(
        (os.path.basename(bam), stats[bam]['AlignedBases']) for bam in stats['BamFiles'])
    plotter.plot_bars_simple(
        aligned_bases, title="Aligned bases", xlab="BAM", ylab="Bases", auto_limit=False)

    aligned_queries = OrderedDict(
        (os.path.basename(bam), stats[bam]['AlignedQueries']) for bam in stats['BamFiles'])
    plotter.plot_bars_simple(
        aligned_queries, title="Aligned queries", xlab="BAM", ylab="Bases", auto_limit=False)

    unaligned_queries = OrderedDict(
        (os.path.basename(bam), stats[bam]['UnalignedQueries']) for bam in stats['BamFiles'])
    plotter.plot_bars_simple(
        unaligned_queries, title="Unaligned queries", xlab="BAM", ylab="Bases", auto_limit=False)

    unaligned_queries = OrderedDict(
        (os.path.basename(bam), stats[bam]['HardClippedBases']) for bam in stats['BamFiles'])
    plotter.plot_bars_simple(
        unaligned_queries, title="Hard clipped bases", xlab="BAM", ylab="Bases", auto_limit=False)

    base_stats = OrderedDict((field, stats[field])
                             for field in ('CommonAlignedBases', 'CommonMatchingBases'))
    plotter.plot_bars_simple(
        base_stats, title="Common aligned base statistics", xlab="Field", ylab="Count", auto_limit=False)

    sim_stats = OrderedDict((field, stats[field]) for field in ['AlignedSimilarity'])
    plotter.plot_bars_simple(
        sim_stats, title="Proportion of bases with matching alignment ({})".format(sim_stats.values()[0]), xlab="Field", ylab="Count", auto_limit=False)

    plotter.plot_histograms({'PerQuerySim': stats[
                            'PerQueryBaseSim']}, title="Distribution of percent bases with matched alignment",
                            xlab="Percent bases with matched alignment", ylab="Count", legend=False)

    plotter.plot_histograms({'PerQuerySimClipped': stats[
                            'PerQueryBaseSimClipped']}, title="Distribution of percent bases with matched alignment (with clipping)",
                            xlab="Percent bases with matched alignment", ylab="Count", legend=False)

    plotter.close()

    if args.p is not None:
        misc.pickle_dump(dict(stats), args.p)

    if args.t is not None:
        data_map = stats.copy()
        del data_map['PerQueryBaseSim']
        del data_map['PerQueryBaseSimClipped']
        for bam in data_map['BamFiles']:
            del data_map[bam]
        del data_map['BamFiles']
        data_map = OrderedDict((key, [value]) for key, value in six.iteritems(data_map))
        data_frame = pd.DataFrame(data_map)
        data_frame.to_csv(args.t, sep="\t", index=False)
