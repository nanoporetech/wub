#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

from wub.util import seq as seq_util
import pandas as pd
from collections import OrderedDict

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Generate a table of read names and mean quality values.')
parser.add_argument(
    '-t', metavar='tsv', type=str, help="Output tab separated file.", default='fastq_qual_tab.tsv')
parser.add_argument('input_fastq', nargs='?', help='Input fastq (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)


if __name__ == '__main__':
    args = parser.parse_args()

    input_iterator = seq_util.read_seq_records(
        args.input_fastq, format='fastq')

    read = []
    mean_qualities = []

    for record in input_iterator:
        read.append(record.id)
        mean_quality = seq_util.mean_qscore(record.letter_annotations["phred_quality"], qround=False)
        mean_qualities.append(mean_quality)

    df = pd.DataFrame(OrderedDict([('Read', read), ('MeanQual', mean_qualities)]))
    df = df.set_index("Read")
    df.to_csv(args.t, sep="\t")
