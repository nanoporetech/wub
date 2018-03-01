#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import pandas as pd
from collections import OrderedDict
from wub.util import seq as seq_util
import datetime

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Produce a tab separated file with read start times, read and channel numbers sorted by start time.""")
parser.add_argument(
    '-t', metavar='read_tsv', type=str, default="fastq_time_tab.tsv", help="Tab separated file to save read time table.", required=False)
parser.add_argument(
    'fastq', metavar='fastq', type=str, help="Input fastq file.")


if __name__ == '__main__':
    args = parser.parse_args()

    input_iterator = seq_util.read_seq_records(args.fastq, format='fastq')

    read, read_nr, channel, start, length = [], [], [], [], []

    for rec in input_iterator:
        read.append(rec.id)
        length.append(len(rec.seq))
        desc = rec.description.split()

        # Parse out read number:
        tmp_read_nr = int(desc[2].split("=")[1])
        read_nr.append(tmp_read_nr)

        # Parse out channel:
        tmp_channel = int(desc[3].split("=")[1])
        channel.append(tmp_channel)

        # Parse out start time:
        tmp_start = desc[4].split("=")[1]
        tmp_start = datetime.datetime.strptime(tmp_start, "%Y-%m-%dT%H:%M:%SZ")
        start.append(tmp_start)

    df = pd.DataFrame(OrderedDict([('Read', read), ('Channel', channel),
                                   ('ReadNumber', read_nr), ('StartTime', start), ("ReadLength", length)]))

    df.sort_values(by="StartTime", inplace=True)
    df = df.set_index("StartTime")

    df.to_csv(args.t, sep="\t")
