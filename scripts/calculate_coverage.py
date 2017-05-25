#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

from wub.util import seq as seq_util
from wub.util import misc


# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Calculate total number of bases and genome coverage if genome size is given.')
parser.add_argument(
    '-f', metavar='format', type=str, help="Input format (fastq).", default='fastq')
parser.add_argument(
    '-s', metavar='genome_size', type=int, help="Genome size (None).", default=None)
parser.add_argument(
            '-p', metavar='results_pickle', type=str, help="Save pickled results in this file.", default=None)
parser.add_argument('input_fastx', nargs='?', help='Input (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)


if __name__ == '__main__':
    args = parser.parse_args()

    in_format = args.f
    input_iterator = seq_util.read_seq_records(
        args.input_fastx, format=in_format)

    total_bases = 0
    for record in input_iterator:
        total_bases += len(record)
    results = {'total_bases': total_bases}
    print("Total bases\t{}".format(total_bases))

    if args.s is not None:
        results['genome_size'] = args.s
        results['coverage'] = float(total_bases) / args.s
        print("Genome size\t{}".format(results['genome_size']))
        print("Coverage\t{}".format(results['coverage']))

    if args.p is not None:
        misc.pickle_dump(results, args.p)
