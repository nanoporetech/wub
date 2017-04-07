#!/usr/bin/env python
# -*- coding: utf-8 -*-

import six
import argparse
import numpy as np
from collections import OrderedDict, defaultdict
import pandas as pd
from wub.vis import report
from wub.util import misc
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import seaborn as sns
warnings.resetwarnings()
_ = sns

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Compare alignment QC statistics of multiple samples.

    It takes a list of pickle files produced by `bam_alignment_qc.py` and produces plots comparing the following properties
    of the input samples:
        * Number of mapped reads.
        * Number of unmapped reads.
        * Distribution of mean quality values in the unaligned fraction.
        * Distribution of mean quality values in the aligned fraction.
        * Distribution of read lengths in the unaligned fraction.
        * Distribution of read lengths in the aligned fraction.
        * Distribution of alignment lengths.
        * Distribution of mapping qualities.
        * Alignment accuracy.
        * Alignment identity.
        * Distribution of deletion lengths.
        * Distribution of insertion lengths.

    Per reference plots (can be disabled by -x):
        * Relative coverage across reference.
        * Mean qualities per position.

    """)
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF (bam_multi_qc.pdf).", default="bam_multi_qc.pdf")
parser.add_argument(
    '-x', action="store_true", help="Do not plot reference statistics.", default=False)
parser.add_argument(
    'pickles', metavar='input_pickles', nargs='*', type=str, help="Input pickles.")


def load_stats(pickles):
    """Load statistics from pickle files.

    :param pickles: List of pickle files.
    :returns: OrderedDict of stats per dataset.
    :rtype: OrderedDict
    """
    stats = OrderedDict()
    for pickle_file in pickles:
        data = misc.pickle_load(pickle_file)
        stats[data['tag']] = data
    return stats


def normalise_dict(d):
    """ Normalise nested dictionary such as the values in the inner dictionary sum to one.

    :param d: Nested dictionary.
    :returns: Normalised nested dictionary.
    :rtype: dict
    """
    total = float(sum(six.itervalues(d)))
    for k, v in six.iteritems(d):
        d[k] = v / total
    return d


def mean_dict(d):
    """Calculate means from values in a dictionary.
    :param d: A dictionary of lists.
    :returns: A dictionary of means.
    :rtype: dict
    """
    for k, v in six.iteritems(d):
        d[k] = np.mean(v)
    return d


def plot_bars_sns(data_map, title, xlab, ylab, plotter):
    """Barplot using seaborn backend.
    :param data_map: A dictionary of labels and values.
    :param title: Plot title.
    :param xlab: X axis label.
    :param ylab: Y axis label.
    :param plotter: A wub.vis.report.Report instance.
    """
    data = pd.DataFrame({'Value': data_map.values(), 'Label': data_map.keys(),
                         'x': np.arange(len(data_map))})
    ax = sns.barplot(x="x", y="Value", hue="Label", data=data, hue_order=data_map.keys())
    ax.set_title(title)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_xticks([])
    plotter.pages.savefig()
    plotter.plt.clf()


if __name__ == '__main__':
    args = parser.parse_args()
    plotter = report.Report(args.r)

    stats = load_stats(args.pickles)

    # Plot mapping statistics:
    mapping_stat_types = OrderedDict([
        ('mapped', {'title': "Number of mapped reads",
                             'xlab': 'Dataset', 'ylab': 'Count'}),
        ('unmapped', {
            'title': "Number of unmapped reads", 'xlab': 'Dataset', 'ylab': 'Count'}), ])

    for field, props in six.iteritems(mapping_stat_types):
        data_map = OrderedDict(
            [(tag, d['read_stats'][field]) for tag, d in six.iteritems(stats)])
        plot_bars_sns(
            data_map, title=props['title'], xlab=props['xlab'], ylab=props['ylab'], plotter=plotter)

    # Plot read statistics:
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

    for field, props in six.iteritems(read_stat_types):
        data_map = OrderedDict(
            [(tag, d['read_stats'][field]) for tag, d in six.iteritems(stats)])
        plotter.plot_histograms(
            data_map, title=props['title'], xlab=props['xlab'], ylab=props['ylab'])

    # Plot basewise statistics:
    base_stat_types = OrderedDict([
        ('accuracy', {'title': "Alignment accuracy",
                      'xlab': 'Dataset', 'ylab': 'Accuracy'}),
        ('identity', {
            'title': "Alignment identity", 'xlab': 'Dataset', 'ylab': 'Identity'}), ])

    for field, props in six.iteritems(base_stat_types):
        data_map = OrderedDict(
            [(tag, d['base_stats'][field]) for tag, d in six.iteritems(stats)])
        plot_bars_sns(
            data_map, title=props['title'], xlab=props['xlab'], ylab=props['ylab'], plotter=plotter)

    # Plot indel statistics:
    indel_stat_types = OrderedDict([
        ('deletion_lengths', {'title': "Distribution of deletion lengths",
                              'xlab': 'Length', 'ylab': 'Count'}),
        ('insertion_lengths', {
            'title': "Distribution of insertion lengths", 'xlab': 'Length', 'ylab': 'Count'})])

    for field, props in six.iteritems(indel_stat_types):
        data_map = OrderedDict(
            [(tag, d['indel_stats'][field]) for tag, d in six.iteritems(stats)])
        plotter.plot_dicts(
            data_map, title=props['title'], xlab=props['xlab'], ylab=props['ylab'], hist_style=True)

    # Plot per-reference statistics:
    if args.x is not True:
        ref_stat_types = OrderedDict([
            ('coverage', {'title': "Relative coverage: ",
                          'xlab': 'Position', 'ylab': 'Relative coverage'}),
            ('qualities', {
                'title': "Mean qualities per position: ", 'xlab': 'Position', 'ylab': 'Mean quality'}), ])
        cov_stats = defaultdict(dict)

        for tag, d in six.iteritems(stats):
            for ref, cov in six.iteritems(d['pileup_stats']['coverage']):
                cov_stats[ref][tag] = normalise_dict(cov)

        for ref, data_map in six.iteritems(cov_stats):
            plotter.plot_dicts(data_map, title=ref_stat_types['coverage'][
                               'title'] + ref, xlab=ref_stat_types['coverage']['xlab'], ylab=ref_stat_types['coverage']['ylab'])

        qual_stats = defaultdict(dict)
        for tag, d in six.iteritems(stats):
            for ref, quals in six.iteritems(d['pileup_stats']['qualities']):
                cov_stats[ref][tag] = mean_dict(quals)

        for ref, data_map in six.iteritems(cov_stats):
            plotter.plot_dicts(data_map, title=ref_stat_types['qualities'][
                               'title'] + ref, xlab=ref_stat_types['qualities']['xlab'], ylab=ref_stat_types['qualities']['ylab'])

    plotter.close()
