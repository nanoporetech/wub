#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse

from wub.util import cmd as cmd_util
from wub.wrappers import dnadiff

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Compare a set of reference sequences (genome) to another set (target assembly) using mummer's dnadiff.
    """)
parser.add_argument(
    '-d', metavar='work_dir', type=str, help="Use this working directory instead of a temporary directory (None).", default=None)
parser.add_argument(
    '-k', action="store_true", help="Keep dnadiff result files (False).", default=False)
parser.add_argument(
    '-t', metavar='details_tsv', type=str, help="Save details of dnadiff comparison in this tab-separated file (None).", default=None)
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report with alignment details plot (None).", default=None)
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

    print dnadiff.dnadiff(args.ref, args.target, args.d, not args.k)
