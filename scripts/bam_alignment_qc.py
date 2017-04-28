#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import six
import argparse
import tqdm

import os
import numpy as np
from collections import OrderedDict, defaultdict
import itertools

from wub.util import misc
from wub.vis import report
from wub.bam import stats as bam_stats
from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Produce alignment based QC plots of the input BAM file.
    The input BAM file must be sorted by coordinates and indexed.

    It produces the following global plots:
        * Read statistics: number of mapped, unmapped and low mapping quality reads.
        * Distribution of mean quality values in the mapped and unmapped fractions.
        * Distribution of read lengths in the unmapped fraction.
        * Distribution of read lengths in the mapped fraction.
        * Distribution of read lengths in the mapping with quality less than -q
        * Distribution of alignment lengths.
        * Distribution of mapping qualities.
        * Plot of alignment lengths vs. mean base qualities.
        * Basewise statistics: total alignment length, number of insertions, deleltions, matches and mismatches.
        * Precision statistics: accuracy and identity.
        * Frequency of errors in the context specifed by the left and right context sizes (-n). Definition of context: for substitutions the event is happening from the "central base", in the case of indels the events are located between the central base and the base before. The columns of the heatmap are normalised to sum to one and then the diagonal element are set to zero.
        * Distribution of deletion lengths.
        * Distribution of insertion lengths.
        * Base composition of insertions.

    The following plots are produced for every reference unless disabled via -x:
        * Distribution of quality values across the reference as a heatmap.
        * Mean quality values across the reference.
        * Base coverage across the reference.

    The tool saves the gathered statistics in a pickle file, which can be fed to `bam_multi_qc.py` to compare different samples.
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
    '-q', metavar='aqual', type=int, default=0, help="Minimum alignment quality (0).")
parser.add_argument(
    '-i', metavar='qual_ints', type=int, default=6, help="Number of quality intervals (6).")
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF (bam_alignment_qc.pdf).", default="bam_alignment_qc.pdf")
parser.add_argument(
    '-p', metavar='results_pickle', type=str, help="Save pickled results in this file (bam_alignment_qc.pk).", default="bam_alignment_qc.pk")
parser.add_argument(
    '-Q', action="store_true", help="Be quiet and do not show progress bars.", default=False)
parser.add_argument(
    'bam', metavar='bam', type=str, help="Input BAM file.")


def find_max_pos(d):
    """ Find maximum position in reference.

    :param d: Dictionary wit position-based values.
    :returns: Maximum position.
    :rtype: int
    """
    return max(list(d.keys()))


def find_max_qual(d):
    """ Find maximum quality value in parsed dataset.

    :param d: Dictionary with quality values.
    :returns: Maximum quality value across all positions.
    :rtype: int
    """
    return max(itertools.chain.from_iterable(list(six.itervalues(d))))


def mean_qual_per_pos(d):
    """ Calculate mean quality scores per position.

    :param d: Dictionary of quality values per positions.
    :returns: Dictionary with mean quality score per position.
    :rtype: dict
    """
    mq = {}
    for pos, quals in six.iteritems(d):
        mq[pos] = seq_util.mean_qscore(quals, qround=False)
    return mq


def stats_to_matrix(rst):
    """Convert dictionary with quality values per position to data matrix.

    :param rst: Dictionary of quality values per position.
    :returns: Matrix of quality values per position.
    :rtype: np.ndarray
    """
    max_pos = find_max_pos(rst)
    max_qual = find_max_qual(rst)
    mat = np.zeros((max_qual + 1, max_pos + 1), dtype=float)
    for pos, quals in six.iteritems(rst):
        for q in quals:
            mat[q][pos] += 1
    return mat


def normalise_data(d):
    """ Normalise nested dictionary such as the values in the inner dictionary sum to one.

    :param d: Nested dictionary.
    :returns: Normalised nested dictionary.
    :rtype: dict
    """
    nd = {}
    for k, v in six.iteritems(d):
        total = float(sum(v.values()))
        nd[k] = {}
        for ik, iv in six.iteritems(v):
            nd[k][ik] = iv / total
    return nd


def ref_qual_qc(st, report, verbose):
    """ Plot per reference statistics.

    :param st: Dictionary with statistics.
    :param report: Report object.
    :param verbose: Verbosity level.
    """
    quals_iter = six.iteritems(st['qualities'])
    if verbose:
        print("Generating per-reference plots.")
        quals_iter = tqdm.tqdm(quals_iter, total=len(st['qualities']))

    for ref, stats in quals_iter:
        mat = stats_to_matrix(stats)
        report.plot_heatmap(mat, title="Quality values across {}".format(
            ref), xlab="Position", ylab="Quality bin")
        mq = mean_qual_per_pos(stats)
        report.plot_dicts(OrderedDict([('Dummy', mq)]), title="Mean quality values across {}".format(
            ref), xlab="Position", ylab="Mean quality", legend=False)
        report.plot_dicts(OrderedDict([('Dummy', st['coverage'][ref])]), title="Base coverage across {}".format(
            ref), xlab="Position", ylab="Coverage", legend=False)


def base_stats_qc(st, report):
    """ Plot base statistics.

    :param st: Dictionary with statistics.
    :param report: Report object.
    """

    bs = st.copy()
    del bs['accuracy']
    del bs['identity']
    plotter.plot_bars_simple(
        bs, title="Basewise statistics", xlab="Type", ylab="Count")
    plotter.plot_bars_simple(OrderedDict([('Identity ({})'.format(st['identity']), st['identity']), ('Accuracy ({})'.format(
        st['accuracy']), st['accuracy'])]), title="Precision statistics", xlab="Type", ylab="Count")


def read_qual_qc(st, report, qual_intervals=5):
    """ Plot histograms of various read statistics.

    :param st: Dictionary with statistics.
    :param report: Report object.
    """
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

    # Process alignment qualities and lengths:
    breaks = np.linspace(start=np.floor(np.min(st['aligned_quals'])), stop=np.ceil(
        np.max(st['aligned_quals'])), num=qual_intervals + 1)
    aq_map = OrderedDict()

    decimals = 3
    for i in range(len(breaks) - 1):
        aq_map[(round(breaks[i], decimals), round(breaks[i + 1], decimals))] = []
    intervals = list(aq_map.keys())

    for i, aln_qual in enumerate(st['aligned_quals']):
        index = np.searchsorted(breaks, aln_qual) - 1
        aq_map[intervals[index]].append(st['alignment_lengths'][i])

    report.plot_boxplots(aq_map, title="Mean base qualities vs. alignment lengths",
                         xlab="Mean base quality", ylab="Alignment lengths", xticks_rotation=45)


def indel_dist_qc(st, report):
    """ Plot histograms of indel distribution statistics.

    :param st: Dictionary with statistics.
    :param report: Report object.
    """
    report.plot_dicts(OrderedDict([('Dummy', st['deletion_lengths'])]), title="Distribution of deletion lengths",
                      xlab="Deletion length", ylab="Count", legend=False, hist_style=True)
    report.plot_dicts(OrderedDict([('Dummy', st['insertion_lengths'])]), title="Distribution of insertion lengths",
                      xlab="Insertion length", ylab="Count", legend=False, hist_style=True)
    plotter.plot_bars_simple(
        st['insertion_composition'], title="Base composition of insertions", xlab="Base", ylab="Count")


def enumerate_contexts(csizes):
    """Enumerate all base context with the specified left and right context sizes. The context length is csizes[0] + csizes[1] + 1.

    :param csizes: Tupple of left and right context sizes.
    :returns: List of all context sorted by the central base (and then lexicographically).
    :rtype: list
    """
    template = [seq_util.bases] * (sum(csizes) + 1)
    contexts = sorted([''.join(c) for c in itertools.product(*template)], key=lambda x:
                      (x[csizes[0]:len(x) - csizes[1]], x[0:csizes[0]], x[len(x) - csizes[1]:len(x)]))
    return contexts


def error_stat_qc(st, report, csizes, ommit_diagonal=False):
    """Plot error counts in context.

    :param st: Dictionary with statistics.
    :param report: Report object.
    """
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
                       colormap=report.plt.cm.Blues, invert_yaxis=True, title="", xlab="From context", ylab="To base")


if __name__ == '__main__':
    args = parser.parse_args()
    verbose = not args.Q
    tag = args.t
    if tag is None:
        tag = os.path.basename(args.bam)
    context_sizes = args.n.split(",")
    context_sizes = (int(context_sizes[0]), int(context_sizes[1]))

    plotter = report.Report(args.r)

    references = seq_util.read_seq_records_dict(args.f)

    err_read_stats = bam_stats.error_and_read_stats(
        args.bam, references, region=args.c, context_sizes=context_sizes, min_aqual=args.q, verbose=verbose)
    read_stats = err_read_stats['read_stats']
    error_stats = err_read_stats['events']
    base_stats = err_read_stats['base_stats']
    indel_stats = err_read_stats['indel_dists']

    read_qual_qc(read_stats, plotter, args.i)
    base_stats_qc(base_stats, plotter)
    error_stat_qc(error_stats, plotter, context_sizes, ommit_diagonal=True)
    indel_dist_qc(indel_stats, plotter)

    pileup_stats = None
    if not args.x:
        pileup_stats = bam_stats.pileup_stats(args.bam, region=args.c, verbose=verbose)
        ref_qual_qc(pileup_stats, plotter, verbose)

    plotter.close()

    # Dump results of parsing into output pickle:
    rd = {'tag': tag, 'read_stats': read_stats, 'error_stats': error_stats,
          'indel_stats': indel_stats, 'pileup_stats': pileup_stats, 'base_stats': base_stats}
    misc.pickle_dump(rd, args.p)
