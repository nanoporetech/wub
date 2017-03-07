#!/usr/bin/env python
__author__ = 'prughani'

import argparse
import pandas as pd
import os
from wub.read_stats import contig_stats as cstats
from wub.util import misc



def main():

    savepath = args.savepath
    fastx = args.fastx
    tag = args.tag


    if savepath is None:
        savepath = os.getcwd()
    else:
        savepath = misc.mkdir(savepath)


    if tag is None:
        tag = misc.get_fname(fastx)

    if misc._getextension(fastx) == 'fastq':
        fq = True
    else:
        fq = False


    rawdata = cstats.GC_per_read(cstats.readfast(fastx), fq=fq)


    print os.path.join(savepath, '{}_summary.stats'.format(tag))


    if args.raw:
        rawdata.to_csv(os.path.join(savepath, '{}_raw.stats'.format(tag)))

    summary = cstats.get_stats(df=rawdata)
    summary.to_csv(os.path.join(savepath, '{}_summary.stats'.format(tag)))
    print summary.round(2).to_string()

    if args.report:
        from wub.vis import report
        Plotter  = report.Report(os.path.join(savepath, '{}.pdf'.format(tag)))

        rawdata = rawdata.sort_values('Seqlen', ascending=True)

        rawdata['cumsum'] = rawdata["Seqlen"].cumsum()
        rawdata['norm'] = 100.0 * rawdata['cumsum']/rawdata['cumsum'].max()

        Plotter.plot_line(data=rawdata, x='Seqlen', y='norm', title='Normalized cumulative plot', xlab='length (bp)', ylab="normalized (%)",)

        # df1.sort_values('Seqlen', ascending=False)
        # df1["cumsum1"] = df1['Seqlen'].cumsum()
        # Plotter.plot_line(data=rawdata, x='Cumsum1', y=df1.reset_index().index, title='Ordered cumulative sum plot', xlab="contigs ordered largest to smallest", ylab='cumulative sum')

        Plotter.plot_scatter(data=rawdata, x='GC content (%)', y='Seqlen', title='GC content vs length plot',
                             xlab="GC content (%)", ylab="length (bp)", alpha=0.5, ylim=0, xlim=0)
        if 'mean_q' in rawdata:

            Plotter.plot_scatter(data=rawdata, x='mean_q', y='Seqlen', title='Mean Q score vs length',
                                 xlab='Mean Q', ylab='length', alpha=0.5, xlim=rawdata['mean_q'].min()-0.5 , ylim=rawdata['Seqlen'].min()-0.5 )

        Plotter.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Calculates the GC content and N50')

    parser.add_argument('--fastx', '-i',
                        metavar='FILE',
                        required=True,
                        help='input file fastq or fasta')

    parser.add_argument('--raw', '-a',
                        action='store_true',
                        required=False,
                        help='save raw the gc content per read/contig. default[False]')

    parser.add_argument('--savepath', '-s',
                        metavar='DIR',
                        required=False,
                        default=None,
                        help='output dir. default[cwd]')

    parser.add_argument('--report', '-r',
                        # metavar="TRUE",
                        action='store_true',
                        required=False,
                        default=None,
                        help = "Report PDF default[False]")

    parser.add_argument('--tag', '-n',
                        metavar='STR',
                        required=False,
                        default=None,
                        help='output name or tag. default[input name]')

    args = parser.parse_args()
    print args

    main()