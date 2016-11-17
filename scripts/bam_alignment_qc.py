#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import os
import pandas as pd
import numpy as np
from collections import OrderedDict
import itertools

from wub.util import misc
from wub.vis import report
from wub.bam import stats
from wub.util import seq as seq_util
from wub.util import misc

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
    '-n', metavar='context_sizes', type=str, help="Left and right context sizes (1,1).", required=False, default="1,1")
parser.add_argument(
    '-x', action="store_true", help="Do not plot per-reference information.", default=False)
parser.add_argument(
    '-t', metavar='bam_tag', type=str, default=None, help="Dataset tag (BAM basename).", required=False)
parser.add_argument(
    '-q', metavar='aqual', type=int, default=5, help="Minimum alignment quality (5).")
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF (bam_alignment_qc.pdf).", default="bam_alignment_qc.pdf")
parser.add_argument(
    '-p', metavar='results_pickle', type=str, help="Save pickled results in this file (bam_alignment_qc.pk).", default="bam_alignment_qc.pk")
parser.add_argument(
    'bam', metavar='bam', type=str, help="Input BAM file.")


def find_max_pos(d):
    """ Find maximum position in reference. """
    return max(d.keys())


def find_max_qual(d):
    """ Find maximum quality value in parsed dataset. """
    return max(itertools.chain.from_iterable(d.itervalues()))


def mean_qual_per_pos(d):
    """ Calculate mean quality scores per position """
    mq = {}
    for pos, quals in d.iteritems():
        mq[pos] = seq_util.mean_qscore(quals)
    return mq


def stats_to_matrix(rst):
    """ Convert dict with quality values per position to data matrix. """
    max_pos = find_max_pos(rst)
    max_qual = find_max_qual(rst)
    mat = np.zeros((max_qual + 1, max_pos + 1), dtype=float)
    for pos, quals in rst.iteritems():
        for q in quals:
            mat[q][pos] += 1
    return mat


def normalise_data(d):
    nd = {}
    for k, v in d.iteritems():
        total = float(sum(v.values()))
        nd[k] = {}
        for ik, iv in v.iteritems():
            nd[k][ik] = iv / total
    return nd

def ref_qual_qc(st, report):
    """ Plot per reference statistics. """
    ref_qc = {}
    for ref, stats in st['qualities'].iteritems():
        mat = stats_to_matrix(stats)
        report.plot_heatmap(mat, title="Quality values across {}".format(
            ref), xlab="Position", ylab="Quality bin")
        mq = mean_qual_per_pos(stats)
        report.plot_dicts(OrderedDict([('Dummy', mq)]), title="Mean quality values across {}".format(
            ref), xlab="Position", ylab="Mean quality", legend=False)
        report.plot_dicts(OrderedDict([('Dummy', st['coverage'][ref])]), title="Base coverage across {}".format(
            ref), xlab="Position", ylab="Coverage", legend=False)


def read_qual_qc(st, report):
    """ Plot histograms of various read statistics. """
    plotter.plot_bars_simple(OrderedDict([('Unmapped', st['unmapped']), ('Mapped', st[
                             'mapped']), ('Mqual < {}'.format(args.q), len(st['mqfail_alignment_lengths']))]), title="Read statistics", xlab="Type", ylab="Count")
    report.plot_histograms(OrderedDict([('Unmapped', st['unaligned_quals']), ('Mapped', st[
                           'aligned_quals'])]), title="Distribution of mean quality values in fractions", xlab="Mean base quality", ylab="Count")
    report.plot_histograms(OrderedDict([('Unmapped', st[
                           'unaligned_lengths'])]), title="Distribution of read lengths in the unmapped fraction", xlab="Read length", ylab="Count", legend=False)
    report.plot_histograms(OrderedDict([('Mapped', st[
                           'alignment_lengths'])]), title="Distribution of read lengths in the mapped fraction", xlab="Read length", ylab="Count", legend=False)
    report.plot_histograms(OrderedDict([('Mapped', st[
                           'mqfail_alignment_lengths'])]), title="Distribution of read lengths in the mapping with quality less than {}".format(args.q), xlab="Read length", ylab="Count", legend=False)
    report.plot_histograms(OrderedDict([('Mapped', st[
                           'alignment_lengths'])]), title="Distribution of alignment lengths", xlab="Alignment length", ylab="Count", legend=False)
    report.plot_histograms(OrderedDict([('MappingQuals', st[
                           'mapping_quals'])]), title="Distribution of mapping qualities", xlab="Mapping quality", ylab="Count", legend=False)
    report.plot_arrays(OrderedDict([('Dummy', (st['alignment_lengths'], st[
                       'aligned_quals']))]), title="Alignment lengths vs. mean base qualities", xlab="Alignment length", ylab="Mean base quality", legend=False)


def indel_dist_qc(st, report):
    """ Plot histograms of indel distribution statistics. """
    report.plot_dicts(OrderedDict([('Dummy', st['deletion_lengths'])]), title="Distribution of deletion lengths",
                      xlab="Deletion length", ylab="Count", legend=False, hist_style=True)
    report.plot_dicts(OrderedDict([('Dummy', st['insertion_lengths'])]), title="Distribution of insertion lengths",
                      xlab="Insertion length", ylab="Count", legend=False, hist_style=True)
    plotter.plot_bars_simple(
        st['insertion_composition'], title="Base composition of insertions", xlab="Base", ylab="Count")


def enumerate_contexts(csizes):
    template = [seq_util.bases] * (sum(csizes) + 1)
    contexts = sorted([''.join(c) for c in itertools.product(*template)], key=lambda x:
                      (x[csizes[0]:len(x) - csizes[1]], x[0:csizes[0]], x[len(x) - csizes[1]:len(x)]))
    return contexts


def error_stat_qc(st, report, csizes, ommit_diagonal=False):
    """Plot error counts in context."""
    contexts = enumerate_contexts(csizes)
    bases_plus = list(seq_util.bases)
    bases_plus.extend(['-', '*'])
    nd = normalise_data(st) 

    z = np.zeros((len(bases_plus), len(contexts)), dtype=float)
    for conti, cont in enumerate(contexts):
        for bi, b in enumerate(bases_plus):
            z[bi][conti] = nd[cont][b]
            central_base = cont[csizes[0]:len(cont) - csizes[1]]
            if central_base == b:
                if ommit_diagonal:
                    z[bi][conti] = 0.0

    report.plot_pcolor(z, xticks=contexts, yticks=bases_plus,
                       colormap=report.plt.cm.Blues, invert_yaxis=True)


if __name__ == '__main__':
    args = parser.parse_args()
    tag = args.t
    if tag is None:
        tag = os.path.basename(args.bam)
    context_sizes = args.n.split(",")
    context_sizes = (int(context_sizes[0]), int(context_sizes[1]))

    plotter = report.Report(args.r)

    references = seq_util.read_seq_records_dict(args.f)
    err_read_stats = stats.error_and_read_stats(
        args.bam, references, region=args.c, context_sizes=context_sizes)
    read_stats = err_read_stats['read_stats']
    error_stats = err_read_stats['events']
    indel_stats = err_read_stats['indel_dists']

    read_qual_qc(read_stats, plotter)
    error_stat_qc(error_stats, plotter, context_sizes, ommit_diagonal=True)
    indel_dist_qc(indel_stats, plotter)

    pileup_stats = None
    if not args.x:
        pileup_stats = stats.pileup_stats(args.bam, args.c)
        ref_qual_qc(pileup_stats, plotter)

    plotter.close()

    # Dump results of parsing into output pickle:
    rd = {'tag': tag, 'read_stats': read_stats, 'error_stats': error_stats,
          'indel_stats': indel_stats, 'pileup_stats': pileup_stats}
    misc.pickle_dump(rd, args.p)
