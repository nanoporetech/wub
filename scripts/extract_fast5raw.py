#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
from collections import OrderedDict

from wub.util import fast5
from wub.vis import report

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Extract raw data from FAST5 file and save as tab separated file. Plot signal if requested.
    """)
parser.add_argument(
    '-t', metavar='tsv', type=str, default="fast5raw.tsv", help="Output tab separated file (fast5raw.tsv).")
parser.add_argument(
    '-r', metavar='report_pdf', type=str, default=None, help="Output tab separated file (None).")
parser.add_argument(
    '-z', metavar='zoom', type=str, default=None, help="Zoom on region 'start:end' (None).")
parser.add_argument(
    'fast5', metavar='fast5', type=str, help="FAST5 file.")


if __name__ == '__main__':
    args = parser.parse_args()

    _, channel, read_number, meta = fast5.load_read_data(args.fast5, raw=True)
    raw_data = meta['raw']

    start, end = 0, len(raw_data)
    if args.z is not None:
        start, end = [int(x) for x in args.z.split(":")]
        raw_data = raw_data[start:end]
    indexes = np.arange(start, end)

    raw_frame = pd.DataFrame(OrderedDict([('Index', indexes), ('RawData', raw_data)]))
    del raw_data
    del indexes
    raw_frame.to_csv(args.t, sep="\t", index=False)

    if args.r is not None:
        plotter = report.Report(args.r)
        plotter.plot_arrays(data_map={'Data': (raw_frame['Index'], raw_frame[
                            'RawData'])}, marker='-', title="Scaled raw signal", xlab="Index", ylab="Value", legend=False)
        plotter.close()
