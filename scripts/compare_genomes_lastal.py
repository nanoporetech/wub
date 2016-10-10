#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse

from wub.mappers import lastal
from wub.util import parse
from wub.vis import report
import pandas as pd

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Compare a set of reference sequences (genome) to another set (target assembly) using lastal alignment.
    Accuracy is the total number of matched bases divided by total alignment length. Coverage is total reference covered
    by alignment divided by total length of reference.

    Caveats:
     - The lastal alignments are filtered by default (use -f to disable) so only the best scoring alignment is kept per query. Hence some shorter valid
     alignments might be discarded causing an underestimation of coverage.
     - The estimated accuracy is dependent on the scoring of gaps and mismatches. By default gap open and gap extend penalties are set to equal.
    """)
parser.add_argument(
    '-l', metavar='lastal_args', type=str, help="Parameters passed to lastal in the <arg>:value,... format (a:1,b:1).", default="a:1,b:1")
parser.add_argument(
    '-t', metavar='details_tsv', type=str, help="Save details of lastal alignment in this tab-separated file (None).", default=None)
parser.add_argument(
    '-f', help="Do *not* filter for best alignment per query.", default=False, action="store_true")
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report with alignment details plot (None).", default=None)
parser.add_argument('ref', metavar='reference_fasta', type=str, help="Reference fasta.")
parser.add_argument('target', metavar='target_fasta', type=str, help="Target fasta.")

if __name__ == '__main__':
    args = parser.parse_args()

    filter_alignments = not args.f
    lastal_args = parse.args_string__to_dict(args.l)
    stats = lastal.compare_genomes_lastal(args.ref, args.target, lastal_options=lastal_args, filter_alns=filter_alignments, cleanup=True)

    global_accuracy = (stats['aln_length'].sum() - stats['substitutions'].sum() -
                       stats['deletions'].sum() - stats['insertions'].sum()) / float(stats['aln_length'].sum())
    global_coverage = stats['ref_aln_len'].sum() / float(stats['ref_len'].sum())

    sys.stdout.write("Accuracy\tCoverage\n")
    sys.stdout.write("{}\t{}\n".format(global_accuracy, global_coverage))

    if args.t is not None:
        stats.to_csv(args.t, sep='\t')

    if args.r is not None:
        plotter = report.Report(args.r)
        data = {'': (stats['coverage'], stats['accuracy'])}
        plotter.plot_arrays(
            data, title="Alignment properties", xlab='Coverage', ylab='Accuracy', legend=False)
        plotter.close()
