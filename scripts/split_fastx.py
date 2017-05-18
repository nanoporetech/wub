#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
from os import path

from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Split sequence records in file to one record per file.')
parser.add_argument(
    '-i', metavar='in_format', type=str, help="Input format (fastq).", default='fastq')
parser.add_argument(
    '-o', metavar='out_format', type=str, help="Output format (fastq).", default='fastq')
parser.add_argument('input_fastx', nargs='?', help='Input file (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_dir', nargs='?', help='Output directory (default: .)', default='.')


if __name__ == '__main__':
    args = parser.parse_args()

    input_iterator = seq_util.read_seq_records(
        args.input_fastx, format=args.i)

    for record in input_iterator:
        bn = path.basename(args.input_fastx.name)
        fh = open(path.join(args.output_dir, "{}_{}".format(record.id, bn)), 'w')
        seq_util.write_seq_records([record], fh, format=args.o)
        fh.flush()
        fh.close()
