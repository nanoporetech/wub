#!/usr/bin/env python
# -*- coding: utf-8 -*-

import six
import argparse
import sys
import pandas as pd
from collections import OrderedDict, defaultdict
from wub.bam import read_counter
from wub.util import misc
from wub.util import seq as seq_util

import tqdm

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Count reads mapping to each reference in a BAM file.""")
parser.add_argument(
    '-a', metavar='min_aqual', type=int, help="Minimum mapping quality (0).", default=0)
parser.add_argument(
    '-f', metavar='in_format', type=str, help="Input format (BAM).", default='BAM')
parser.add_argument(
    '-z', metavar='ref_fasta', type=str, help="Reference fasta. GC content and length columns are added if present (None).", default=None)
parser.add_argument(
    '-k', metavar="words", type=str, help="Include word frequencies of specifed length in output (1,2).", default="1,2")
parser.add_argument(
    '-p', metavar='results_pickle', type=str, help="Save pickled results in this file (None).", default=None)
parser.add_argument(
    '-t', metavar='tsv_file', type=str, help="Save results in tsv format in this file (bam_count_reads.tsv).", default="bam_count_reads.tsv")
parser.add_argument(
    '-Q', action="store_true", help="Be quiet and do not print progress bar (False).", default=False)
parser.add_argument(
    '-R', action="store_true", help="Count reads from SAM stream in stdin. Only read count fields are written. Header required! (False).", default=False)
parser.add_argument(
    '-F', metavar='yield_freq', type=int, help="Yield counts after every -Fth mapped record when doing online counting (100).", default=100)
parser.add_argument('bam', nargs='?', help='Input file (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)


def _offline_counter(args):
    """ Offline counting from SAM/BAM file. """
    # Offline counting from SAM/BAM file:
    counts = read_counter.count_reads(
        args.bam.name, in_format=args.f, min_aln_qual=args.a, verbose=not args.Q)
    counts = OrderedDict(six.iteritems(counts))

    calc_words = [int(k) for k in args.k.split(",")]

    data = OrderedDict()

    # Calculate sequence properties:
    if args.z is not None:
        lengths, gc_contents, word_freqs = {}, {}, defaultdict(
            lambda: defaultdict(dict))
        ref_iter = seq_util.read_seq_records(args.z)
        if not args.Q:
            sys.stderr.write("Calculating sequence features:\n")
            ref_iter = tqdm.tqdm(ref_iter)

        for ref in ref_iter:
            # Augment counts dictionary with missing reference entries:
            if ref.id not in counts:
                counts[ref.id] = 0
            lengths[ref.id] = len(ref)
            gc_contents[ref.id] = seq_util.gc_content(str(ref.seq))
            if args.k is not None:
                for word_size in calc_words:
                    bf = seq_util.word_composition(ref.seq, word_size)
                    for word, count in six.iteritems(bf):
                        word_freqs[word_size][ref.id][
                            word] = float(count) / len(ref)

        data['Length'] = [lengths[tr] for tr in six.iterkeys(counts)]
        data['GC_content'] = [gc_contents[tr] for tr in six.iterkeys(counts)]

    data['Reference'] = counts.keys()
    data['Count'] = counts.values()

    # Calculate word frequencies:
    if args.k is not None and args.z:
        for ks in calc_words:
            for word in word_freqs[ks].values().next().keys():
                tmp = []
                for ref in counts.keys():
                    tmp.append(word_freqs[ks][ref][word])
                data[word] = tmp

    data_frame = pd.DataFrame(data)
    data_frame = data_frame.sort(['Count', 'Reference'], ascending=False)

    if args.t is not None:
        data_frame.to_csv(args.t, sep='\t', index=False)

    if args.p is not None:
        misc.pickle_dump(data, args.p)


def _online_counter(args):
    """ Online counting from SAM stream. """
    # Open counts stream:
    counts_iter = read_counter.count_reads_realtime(
        alignment_file='-', in_format=args.f, min_aln_qual=args.a, verbose=not args.Q, yield_freq=args.F)

    for counts in counts_iter:
        data_frame = pd.DataFrame(
            OrderedDict([('Reference', counts.keys()), ('Count', counts.values())]))
        data_frame = data_frame.sort(['Count', 'Reference'], ascending=False)

        if args.t is not None:
            data_frame.to_csv(args.t, sep='\t', index=False)
        if args.p is not None:
            misc.pickle_dump(counts, args.p)


if __name__ == '__main__':
    args = parser.parse_args()

    if not args.R:
        # Offline counting from SAM/BAM file:
        if args.bam == sys.stdin:
            raise Exception("Input file not specified!")
        _offline_counter(args)
    else:
        # Online counting from SAM on stdin.
        _online_counter(args)
