#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import pysam
from wub.bam import filter as bam_filter
from wub.bam import common as bam_common

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Filter SAM/BAM records by score or other criteria.
    WARNING: the input records must be sorted by name or the filtering will not work
    as expected.
    """)
parser.add_argument(
    '-f', metavar='format', type=str, help="Input/output format (SAM).", default='SAM')
parser.add_argument(
    '-s', metavar='strategy', type=str, help="Filtering strategy: top_per_query, query_coverage (top_per_query).",
    default="top_per_query", choices=['top_per_query', 'query_coverage'])
parser.add_argument(
    '-q', metavar='query_cover', type=float, help="Minimum query coverage fraction (0.8).", default=0.8)
parser.add_argument(
    'infile', metavar='input_file', type=str, help="Input file.")
parser.add_argument(
    'outfile', metavar='output_file', type=str, help="Output SAM file.")

if __name__ == '__main__':
    args = parser.parse_args()

    input_iter = bam_common.pysam_open(args.infile, args.f)

    if args.s == 'top_per_query':
        output_iter = bam_filter.filter_top_per_query(input_iter.fetch(until_eof=True))
    elif args.s == 'query_coverage':
        output_iter = bam_filter.filter_query_coverage(input_iter.fetch(until_eof=True), args.q)
    else:
        raise Exception('Filtering strategy not implemented!')

    writer = pysam.AlignmentFile(args.outfile, "wh", template=input_iter, header=input_iter.header)
    for record in output_iter:
        writer.write(record)

    writer.close()
