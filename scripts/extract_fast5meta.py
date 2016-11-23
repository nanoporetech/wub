#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import tqdm
import concurrent.futures

import os
from os.path import join
import pandas as pd
import numpy as np
from collections import defaultdict
import itertools

from wub.util import fast5


def all_fast5(fast5_dir):
    all_files = []
    for root, dirs, files in os.walk(fast5_dir):
        pool = [join(root, name) for name in files if name.endswith('.fast5')]
        all_files.extend(pool)
    return all_files


def get_fast5_meta(fname):
    _, channel, read_number, meta = fast5.load_read_data(fname)
    meta['channel'] = channel
    meta['read_number'] = read_number
    meta['file'] = fname
    return meta

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Extract metadata from FAST5 files and save as tab separated file.
    """)
parser.add_argument(
    '-Q', action="store_true", help="Be quiet and do not show progress bars.", default=False)
parser.add_argument(
    '-t', metavar='max_threads', type=int, default=None, help="Maximum number of worker threads (None).")
parser.add_argument(
    'fast5dir', metavar='fast5dir', type=str, help="Directory with FAST5 files.")
parser.add_argument(
    'tsv', metavar='tsv', type=str, help="Output tab separated file.")


if __name__ == '__main__':
    args = parser.parse_args()
    verbose = not args.Q

    if verbose:
        print "Locating FAST5 files in directory: {}".format(args.fast5dir)

    files = all_fast5(args.fast5dir)

    stats = defaultdict(list)

    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        # Launch worker threads:
        future_to_file = {executor.submit(get_fast5_meta, f5): fast5 for f5 in files}
        future_iter = concurrent.futures.as_completed(future_to_file)

        if verbose:
            print "Parsing FAST5 files:"
            future_iter = tqdm.tqdm(future_iter, total=len(files))

        # Iterate over results:
        for future in future_iter:
            f5 = future_to_file[future]
            try:
                data = future.result()
                for name, value in data.iteritems():
                    stats[name].append(value)
            except Exception as exc:
                print('%r generated an exception: %s' % (f5, exc))
            else:
                pass

    df = pd.DataFrame(stats)
    df.to_csv(args.tsv, sep="\t", index=False)
