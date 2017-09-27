#!/usr/bin/env python
# -*- coding: utf-8 -*-

import six
import argparse
import sys
from os import path

from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Split sequence records in file to one record per file or batches of records.')
parser.add_argument(
    '-i', metavar='in_format', type=str, help="Input format (fastq).", default='fastq')
parser.add_argument(
    '-o', metavar='out_format', type=str, help="Output format (fastq).", default='fastq')
parser.add_argument(
    '-b', metavar='batch_size', type=int, help="Batch size (None).", default=None)
parser.add_argument('input_fastx', nargs='?', help='Input file (default: stdin).',
                    type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('output_dir', nargs='?', help='Output directory (default: .)', default='.')


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.
    Taken from the biopython wiki: http://biopython.org/wiki/Split_large_file

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.

    :param iterator: Input iterator.
    :param batch_size: Batch size.
    :returns: Generator of lists.
    :rtype: generator
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = six.next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


if __name__ == '__main__':
    args = parser.parse_args()

    input_iterator = seq_util.read_seq_records(
        args.input_fastx, format=args.i)

    if args.b is None:
        # Splitting one record per file:
        for record in input_iterator:
            bn = path.basename(args.input_fastx.name)
            ext = bn.rsplit('.', 1)[-1]
            fh = open(path.join(args.output_dir, "{}.{}".format(record.id, ext)), 'w')
            seq_util.write_seq_records([record], fh, format=args.o)
            fh.flush()
            fh.close()
    else:
        # Split into batches:
        input_iterator = batch_iterator(input_iterator, args.b)
        i = 0
        for records in input_iterator:
            bn = path.basename(args.input_fastx.name)
            fh = open(path.join(args.output_dir, "batch_{}_{}".format(i, bn)), 'w')
            seq_util.write_seq_records(records, fh, format=args.o)
            fh.flush()
            fh.close()
            i += 1
