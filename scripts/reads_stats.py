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
    tag =  args.tag

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

    if args.raw:
        rawdata.to_csv(os.path.join(savepath, '{}_raw.stats'.format(tag)))

    summary = cstats.get_stats(df=rawdata)
    summary.to_csv(os.path.join(savepath, '{}_summary.stats'.format(tag)))
    print summary.to_string()



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Calculates the GC content and N50')

    parser.add_argument('--fastx', '-i',
                        metavar='FILE',
                        required=True,
                        help='input file fastq or fasta')

    parser.add_argument('--raw', '-r',
                        action='store_true',
                        required=False,
                        help='save raw the gc content per read/contig. default[False]')

    parser.add_argument('--savepath', '-s',
                        metavar='DIR',
                        required=False,
                        default=None,
                        help='output dir. default[cwd]')

    parser.add_argument('--tag', '-n',
                        metavar='STR',
                        required=False,
                        default=None,
                        help='output name or tag. default[input name]')

    args = parser.parse_args()
    print args

    main()