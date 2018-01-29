#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
from collections import OrderedDict
import pandas as pd
from wub.vis import report
from wub.util import misc
import matplotlib.pyplot as plt
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import seaborn as sns
warnings.resetwarnings()
_ = sns

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Compare alignment QC statistics of multiple samples.
    """)
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF (plot_gffcmp_stats.pdf).", default="plot_gffcmp_stats.pdf")
parser.add_argument(
    '-p', metavar='pickle_out', type=str, help="Output pickle file.", default="plot_gffcmp_stats.pk")
parser.add_argument(
    'input', metavar='input_txt', type=str, help="Input gffcompare stats file.")


def plot_bars_sns(data_map, title, xlab, ylab, plotter):
    """Barplot using seaborn backend.
    :param data_map: A dictionary of labels and values.
    :param title: Plot title.
    :param xlab: X axis label.
    :param ylab: Y axis label.
    :param plotter: A wub.vis.report.Report instance.
    """
    data = pd.DataFrame({'Value': list(data_map.values()), 'Label': list(data_map.keys()),
                         'x': np.arange(len(data_map))})
    ax = sns.barplot(x="x", y="Value", hue="Label", data=data, hue_order=list(data_map.keys()))
    ax.set_title(title)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_xticks([])
    plotter.pages.savefig()
    plotter.plt.clf()


def _parse_stat_line(sl):
    res = {}
    tmp = sl.split(':')[1]
    tmp = tmp.split('|')
    res['sensitivity'] = float(tmp[0].strip())
    res['precision'] = float(tmp[1].strip())
    return res


def _parse_matching_line(l):
    tmp = l.split(':')[1].strip()
    return int(tmp)


def _parse_mn_line(l):
    res = {}
    tmp = l.split(':')[1].strip()
    tmp = tmp.split('/')
    res['value'] = int(tmp[0])
    tmp = tmp[1].split('(')
    res['value_total'] = int(tmp[0].strip())
    res['percent'] = float(tmp[1].split('%)')[0])
    return res


def _parse_total_line(l):
    res = {}
    tmp = l.split(':')[1].strip()
    tmp = tmp.split('in')
    res['transcripts'] = int(tmp[0].strip())
    tmp = tmp[1].split('loci')
    res['loci'] = int(tmp[0].strip())
    tmp = int(tmp[1].split('(')[1].split(' ')[0])
    res['me_transcripts'] = tmp
    return res


def parse_gffcmp_stats(txt):

    sensitivity = []
    precision = []
    level = []

    matching = OrderedDict()

    missed_level = []
    missed = []
    missed_total = []
    missed_percent = []

    novel_level = []
    novel = []
    novel_total = []
    novel_percent = []

    total_target = []
    total_loci = []
    total_transcripts = []
    total_multiexonic = []

    fh = open(txt, 'r')
    for line in fh:
        line = line.strip()
        if len(line) == 0:
            continue
        # Parse totals:
        if line.startswith('#     Query mRNAs'):
            total_target.append('Query')
            r = _parse_total_line(line)
            total_loci.append(r['loci'])
            total_transcripts.append(r['transcripts'])
            total_multiexonic.append(r['me_transcripts'])

        if line.startswith('# Reference mRNAs '):
            total_target.append('Reference')
            r = _parse_total_line(line)
            total_loci.append(r['loci'])
            total_transcripts.append(r['transcripts'])
            total_multiexonic.append(r['me_transcripts'])

        # Parse basic statistics:
        if line.startswith('Base level'):
            st = _parse_stat_line(line)
            level.append('Base')
            sensitivity.append(st['sensitivity'])
            precision.append(st['precision'])
        if line.startswith('Exon level'):
            st = _parse_stat_line(line)
            level.append('Exon')
            sensitivity.append(st['sensitivity'])
            precision.append(st['precision'])
        if line.startswith('Intron level'):
            st = _parse_stat_line(line)
            level.append('Intron')
            sensitivity.append(st['sensitivity'])
            precision.append(st['precision'])
        if line.startswith('Intron chain level'):
            st = _parse_stat_line(line)
            level.append('Intron chain')
            sensitivity.append(st['sensitivity'])
            precision.append(st['precision'])
        if line.startswith('Transcript level'):
            st = _parse_stat_line(line)
            level.append('Transcript')
            sensitivity.append(st['sensitivity'])
            precision.append(st['precision'])
        if line.startswith('Locus level'):
            st = _parse_stat_line(line)
            level.append('Locus')
            sensitivity.append(st['sensitivity'])
            precision.append(st['precision'])

        # Parse match statistics:
        if line.startswith('Matching intron chains'):
            m = _parse_matching_line(line)
            matching['Intron chains'] = [m]
        if line.startswith('Matching transcripts'):
            m = _parse_matching_line(line)
            matching['Transcripts'] = [m]
        if line.startswith('Matching loci'):
            m = _parse_matching_line(line)
            matching['Loci'] = [m]

        # Parse missing statistics:
        if line.startswith('Missed exons'):
            missed_level.append('Exons')
            r = _parse_mn_line(line)
            missed.append(r['value'])
            missed_total.append(r['value_total'])
            missed_percent.append(r['percent'])
        if line.startswith('Missed introns'):
            missed_level.append('Introns')
            r = _parse_mn_line(line)
            missed.append(r['value'])
            missed_total.append(r['value_total'])
            missed_percent.append(r['percent'])
        if line.startswith('Missed loci'):
            missed_level.append('Loci')
            r = _parse_mn_line(line)
            missed.append(r['value'])
            missed_total.append(r['value_total'])
            missed_percent.append(r['percent'])

        # Parse novel statistics:
        if line.startswith('Novel exons'):
            novel_level.append('Exons')
            r = _parse_mn_line(line)
            novel.append(r['value'])
            novel_total.append(r['value_total'])
            novel_percent.append(r['percent'])
        if line.startswith('Novel introns'):
            novel_level.append('Introns')
            r = _parse_mn_line(line)
            novel.append(r['value'])
            novel_total.append(r['value_total'])
            novel_percent.append(r['percent'])
        if line.startswith('Novel loci'):
            novel_level.append('Loci')
            r = _parse_mn_line(line)
            novel.append(r['value'])
            novel_total.append(r['value_total'])
            novel_percent.append(r['percent'])

    fh.close()

    df_stats = pd.DataFrame(OrderedDict(
        [('Sensitivity', sensitivity), ('Precision', precision)]), index=level)
    df_match = pd.DataFrame(matching, index=['Matching'])
    df_miss = pd.DataFrame(OrderedDict(
        [('Total', missed_total), ('Missed', missed), ('Percent missed', missed_percent)]), index=missed_level)
    df_novel = pd.DataFrame(OrderedDict(
        [('Total', novel_total), ('Novel', novel), ('Percent novel', novel_percent)]), index=novel_level)

    df_total = pd.DataFrame(OrderedDict(
        [('Loci', total_loci), ('Transcripts', total_transcripts), ('Multiexonic', total_multiexonic)]), index=total_target)

    return df_stats, df_match, df_miss, df_novel, df_total


if __name__ == '__main__':
    args = parser.parse_args()

    stats, match, miss, novel, total = parse_gffcmp_stats(args.input)

    plotter = report.Report(args.r)

    # Plot overview panel:

    plt.figure(1)
    plt.subplot(2, 2, 1)
    total.plot(ax=plt.gca(), kind='barh', sharex=False, title='Totals')
    plt.tight_layout()

    plt.subplot(2, 2, 2)
    stats.plot(ax=plt.gca(), kind='barh', legend=True, sharex=False, title='Performance').legend(loc='best')
    plt.tight_layout()

    plt.subplot(2, 2, 3)
    miss.copy().drop('Percent missed', axis=1).plot(ax=plt.gca(), kind='barh', legend=True, sharex=False, title='Missed')
    plt.tight_layout()

    plt.subplot(2, 2, 4)
    novel.copy().drop('Percent novel', axis=1).plot(ax=plt.gca(), kind='barh', legend=True, sharex=False, title='Novel')
    plt.tight_layout()
    plotter.pages.savefig()

    # Plot individual panels:

    total.plot(kind='barh', subplots=True, legend=False, sharex=False)
    plt.tight_layout()
    plotter.pages.savefig()

    stats.plot(kind='barh', subplots=True, legend=False, sharex=False)
    plt.tight_layout()
    plotter.pages.savefig()

    match.plot(kind='barh', subplots=True, legend=False)
    plt.tight_layout()
    plotter.pages.savefig()

    miss.plot(kind='barh', subplots=True, legend=False, sharex=False)
    plt.tight_layout()
    plotter.pages.savefig()

    novel.plot(kind='barh', subplots=True, legend=False, sharex=False)
    plt.tight_layout()
    plotter.pages.savefig()

    plotter.close()

    if args.p is not None:
        p = {'total': total, 'stats': stats, 'match': match, 'miss': miss, 'novel': novel}
        misc.pickle_dump(p, args.p)
