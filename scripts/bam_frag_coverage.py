#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import tqdm

import os
import numpy as np
from collections import OrderedDict
import itertools
from Bio import SeqIO

from wub.util import misc
from wub.vis import report
from wub.bam import stats as bam_stats
from wub.util import seq as seq_util

import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import seaborn as sns
warnings.resetwarnings()
_ = sns

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Produce aggregated and individual plots of fragment coverage.
    """)
parser.add_argument(
    '-f', metavar='reference', type=str, help="Reference fasta.", required=True)
parser.add_argument(
    '-c', metavar='region', type=str, help="BAM region (None).", required=False, default=None)
parser.add_argument(
    '-x', action="store_true", help="Plot per-reference information.", default=False)
parser.add_argument(
    '-t', metavar='bam_tag', type=str, default=None, help="Dataset tag (BAM basename).", required=False)
parser.add_argument(
    '-q', metavar='aqual', type=int, default=0, help="Minimum alignment quality (0).")
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF (bam_alignment_qc.pdf).", default="bam_frag_coverage.pdf")
parser.add_argument(
    '-p', metavar='results_pickle', type=str, help="Save pickled results in this file (bam_alignment_qc.pk).", default="bam_frag_coverage.pk")
parser.add_argument(
    '-Q', action="store_true", help="Be quiet and do not show progress bars.", default=False)
parser.add_argument(
    'bam', metavar='bam', type=str, help="Input BAM file.")


def _set_properties_and_close(plotter, fig, title, xlab, ylab):
    """Utility method to set title, axis labels and close the figure.
    """
    plotter.plt.xlabel(xlab)
    plotter.plt.ylabel(ylab)
    plotter.plt.title(title)
    plotter.plt.legend(loc='best')
    plotter.pages.savefig(fig)
    plotter.plt.close(fig)


def _process_ref_coverage(plotter, cov, strand, scale_pos, scale_cov):
    x = np.arange(len(cov), dtype=float)
    if scale_pos:
        x = x / len(x)
    y = cov.astype(float)
    if scale_cov:
        y = y / np.sum(y)
    if strand == 'rev':
        y = -y
    return x, y


def _plot_frag_coverage(st, chroms, plotter, scale_pos=True, scale_cov=False, title="", bins='auto'):
    fig = plotter.plt.figure()

    if bins == 'auto':
        bins = int(np.mean(chroms.values()))

    X = np.arange(bins, dtype=float)
    x = X / len(X)
    cov_fwd = np.zeros((bins), dtype=float)
    cov_rev = np.zeros((bins), dtype=float)

    for chrom in chroms.keys():
        if chrom in st['frags_fwd']:
            fx, fy = _process_ref_coverage(
                plotter, cov=st['frags_fwd'][chrom], strand='fwd', scale_pos=scale_pos, scale_cov=scale_cov)
            cov_fwd += np.interp(x, fx, fy)
        if chrom in st['frags_rev']:
            rx, ry = _process_ref_coverage(
                plotter, cov=st['frags_rev'][chrom], strand='rev', scale_pos=scale_pos, scale_cov=scale_cov)
            cov_rev += np.interp(x, rx, ry)

    plotter.plt.plot(X, np.log(cov_fwd), '-', label='+')
    plotter.plt.plot(X, -np.log(-cov_rev), '-',label='-')
    _set_properties_and_close(plotter, fig, title=title, xlab="Position", ylab="Fragment coverage")
    return {'global_cov_fwd': cov_fwd, 'global_cov_rev': cov_rev}


if __name__ == '__main__':
    args = parser.parse_args()
    verbose = not args.Q
    tag = args.t
    if tag is None:
        tag = os.path.basename(args.bam)

    plotter = report.Report(args.r)

    references = SeqIO.index(args.f, format='fasta')
    chrom_lengths = {name: len(so) for name, so in references.iteritems()}

    st = bam_stats.frag_coverage(args.bam, chrom_lengths, args.c, args.q, verbose)

    res = {'chrom_covs': {}}

    res['global_cov'] = _plot_frag_coverage(
        st, chrom_lengths, plotter, title="Global fragment coverage for {}".format(tag))

    if args.x:
        for chrom, length in chrom_lengths.iteritems():
            res['chrom_covs'][chrom] = _plot_frag_coverage(
                st, {chrom: length}, plotter, title="Fragment coverage for {}:{}".format(tag, chrom))

    plotter.close()

    # Dump results of parsing into output pickle:
    misc.pickle_dump(res, args.p)
