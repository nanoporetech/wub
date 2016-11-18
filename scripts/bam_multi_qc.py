#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import os
import numpy as np
from collections import OrderedDict, defaultdict


from wub.vis import report
from wub.util import misc

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Compare alignment QC statistics of multiple samples.
    """)
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF (bam_multi_qc.pdf).", default="bam_multi_qc.pdf")
parser.add_argument(
    '-x', action="store_true", help="Do not plot reference statistics.", default=False)
parser.add_argument(
    'pickles', metavar='input_pickles', nargs='*', type=str, help="Input pickles.")


def load_stats(pickles):
    stats = {}
    for pickle_file in pickles:
        data = misc.pickle_load(pickle_file)
        stats[data['tag']] = data
    return stats


def normalise_dict(d):
    total = float(sum(d.itervalues()))
    for k, v in d.iteritems():
        d[k] = v / total
    return d


def mean_dict(d):
    for k, v in d.iteritems():
        d[k] = np.mean(v)
    return d


if __name__ == '__main__':
    args = parser.parse_args()
    plotter = report.Report(args.r)

    stats = load_stats(args.pickles)

    mapping_stat_types = OrderedDict([
        ('mapped', {'title': "Number of mapped reads",
                             'xlab': 'Dataset', 'ylab': 'Count'}),
        ('unmapped', {
            'title': "Number of unmapped reads", 'xlab': 'Dataset', 'ylab': 'Count'}), ])

    for field, props in mapping_stat_types.iteritems():
        data_map = OrderedDict(
            [(tag, d['read_stats'][field]) for tag, d in stats.iteritems()])
        plotter.plot_bars_simple(
            data_map, title=props['title'], xlab=props['xlab'], ylab=props['ylab'])

    read_stat_types = OrderedDict([
        ('unaligned_quals', {'title': "Distribution of mean quality values in the unaligned fraction",
                             'xlab': 'Mean base quality', 'ylab': 'Count'}),
        ('aligned_quals', {'title': "Distribution of mean quality values in the aligned fraction",
                           'xlab': 'Mean base quality', 'ylab': 'Count'}),
        ('unaligned_lengths', {
            'title': "Distribution of read lengths in the unaligned fraction", 'xlab': 'Read length', 'ylab': 'Count'}),
        ('aligned_lengths', {
            'title': "Distribution of read lengths in the aligned fraction", 'xlab': 'Read length', 'ylab': 'Count'}),
        ('alignment_lengths', {
            'title': "Distribution of alignment lengths", 'xlab': 'Alignment length', 'ylab': 'Count'}),
        ('mapping_quals', {
            'title': "Distribution of mapping qualities", 'xlab': 'Mapping quality', 'ylab': 'Count'}), ])

    for field, props in read_stat_types.iteritems():
        data_map = OrderedDict(
            [(tag, d['read_stats'][field]) for tag, d in stats.iteritems()])
        plotter.plot_histograms(
            data_map, title=props['title'], xlab=props['xlab'], ylab=props['ylab'])

    base_stat_types = OrderedDict([
        ('accuracy', {'title': "Alignment accuracy",
                      'xlab': 'Dataset', 'ylab': 'Accuracy'}),
        ('identity', {
            'title': "Alignment identity", 'xlab': 'Dataset', 'ylab': 'Identity'}), ])

    for field, props in base_stat_types.iteritems():
        data_map = OrderedDict(
            [(tag, d['base_stats'][field]) for tag, d in stats.iteritems()])
        plotter.plot_bars_simple(
            data_map, title=props['title'], xlab=props['xlab'], ylab=props['ylab'])

    if args.x is not None:
        ref_stat_types = OrderedDict([
            ('coverage', {'title': "Relative coverage: ",
                          'xlab': 'Position', 'ylab': 'Relative coverage'}),
            ('qualities', {
                'title': "Mean qualities per position: ", 'xlab': 'Position', 'ylab': 'Mean quality'}), ])
        cov_stats = defaultdict(dict)
        import pprint
        pp = pprint.PrettyPrinter()

        for tag, d in stats.iteritems():
            for ref, cov in d['pileup_stats']['coverage'].iteritems():
                cov_stats[ref][tag] = normalise_dict(cov)

        for ref, data_map in cov_stats.iteritems():
            plotter.plot_dicts(data_map, title=ref_stat_types['coverage'][
                               'title'] + ref, xlab=ref_stat_types['coverage']['xlab'], ylab=ref_stat_types['coverage']['ylab'])

        qual_stats = defaultdict(dict)
        for tag, d in stats.iteritems():
            for ref, quals in d['pileup_stats']['qualities'].iteritems():
                cov_stats[ref][tag] = mean_dict(quals)

        for ref, data_map in cov_stats.iteritems():
            plotter.plot_dicts(data_map, title=ref_stat_types['qualities'][
                               'title'] + ref, xlab=ref_stat_types['qualities']['xlab'], ylab=ref_stat_types['qualities']['ylab'])

    plotter.close()
