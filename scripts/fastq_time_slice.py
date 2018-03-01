#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

import pandas as pd
from wub.util import seq as seq_util
import datetime

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Filter a fastq file by starting time.""")
parser.add_argument(
    '-t', metavar='time_tsv', type=str, help="Tab separeted file produced by fastq_time_tab.py.", required=True)
parser.add_argument(
    '-s', metavar='start_perc', type=float, help="Start of slice as percent of total time.", required=False, default=0.0)
parser.add_argument(
    '-e', metavar='end_perc', type=float, help="End of slice as percent of total time.", required=False, default=100.0)
parser.add_argument('input_fastq', nargs='?', help='Input fastq (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_fastq', nargs='?', help='Output fastq (default: stdout)',
                    type=argparse.FileType('w'), default=sys.stdout)


def _time_slice(input_iter, start_perc, end_perc, time_df):
    """ Filter for fastq records falling in the specified time range. """
    first = time_df.index.min()
    last = time_df.index.max()

    s = first + ((last - first) * start_perc) / 100.0
    e = first + ((last - first) * end_perc) / 100.0

    for rec in input_iter:
        desc = rec.description.split()
        tmp_start = desc[4].split("=")[1]
        start_time = datetime.datetime.strptime(tmp_start, "%Y-%m-%dT%H:%M:%SZ")
        if start_time >= s and start_time <= e:
            yield rec


if __name__ == '__main__':
    args = parser.parse_args()
    time_df = pd.read_csv(args.t, sep="\t", parse_dates=True, index_col="StartTime")

    input_iterator = seq_util.read_seq_records(args.input_fastq, format='fastq')

    output_iterator = _time_slice(input_iterator, args.s, args.e, time_df)

    seq_util.write_seq_records(output_iterator, args.output_fastq, format='fastq')
