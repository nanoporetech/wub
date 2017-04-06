#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
from __future__ import print_function
from wub.util import cmd as cmd_util
from wub.wrappers import dnadiff
from wub.util import misc

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Compare a set of reference sequences (genome) to another set (target assembly) using mummer's dnadiff.
    It prints the alignment results to stdout. All parsed results can be saved in a pickle file.
    """)
parser.add_argument(
    '-p', metavar='results_pickle', type=str, help="Save pickled results in this file (None).", default=None)
parser.add_argument(
    '-r', metavar='raw_file', type=str, help="Save dnadiff report in this file (None).", default=None)
parser.add_argument(
    '-d', metavar='work_dir', type=str, help="Use this working directory instead of a temporary directory (None).", default=None)
parser.add_argument(
    '-k', action="store_true", help="Keep dnadiff result files (False).", default=False)
parser.add_argument(
    '-v', action="store_true", help="Print out dnadiff output (False).", default=False)
parser.add_argument(
    'ref', metavar='reference_fasta', type=str, help="Reference fasta.")
parser.add_argument(
    'target', metavar='target_fasta', type=str, help="Target fasta.")

if __name__ == '__main__':
    args = parser.parse_args()

    cmd_util.ensure_executable('dnadiff')
    cmd_util.ensure_executable('delta-filter')
    cmd_util.ensure_executable('show-diff')
    cmd_util.ensure_executable('show-snps')
    cmd_util.ensure_executable('show-coords')
    cmd_util.ensure_executable('nucmer')

    results, raw_report, log = dnadiff.dnadiff(args.ref, args.target, args.d, not args.k)

    if args.v:
        sys.stdout.write(log)

    if args.r is not None:
        with open(args.r, 'w') as out_handle:
            out_handle.write(raw_report)

    if args.p is not None:
        misc.pickle_dump(results, args.p)

    for section, properties in results['Alignments'].iteritems():
        print(section, ":\t\tref\tquery")
        for name, prop in properties.iteritems():
            print("\t{}\t{}\t{}".format(name, prop.ref, prop.query))
