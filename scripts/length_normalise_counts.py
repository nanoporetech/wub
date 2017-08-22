#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
from collections import OrderedDict
from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Calculate RPKM values from raw counts and a transcriptome reference.""")
parser.add_argument(
    '-f', metavar='in_trs', type=str, help="Input transcriptome.", required=True)
parser.add_argument('input_counts', nargs=1, help='Input count file.',
                    type=str, default=None)
parser.add_argument('output_count', nargs=1, help='Output RPKM file.',
                    type=str, default=None)


def _load_transcript_lengths(fasta):
    """ Load transcript lengths. """
    res = {}
    for record in seq_util.read_seq_records(fasta):
        res[record.id] = len(record.seq)
    return res


if __name__ == '__main__':
    args = parser.parse_args()

    # Load transcript lengths:
    trs_lens = _load_transcript_lengths(args.f)

    # Load input counts:
    in_df = pd.read_csv(args.input_counts[0], sep="\t")

    # Calculate scaling factor:
    million_factor = np.sum(in_df["Count"]) / float(10**6)

    # Normalise counts:
    refs, rpkms = [], []
    for row in in_df.itertuples():
        refs.append(row.Reference)
        rpkms.append(row.Count / (million_factor *
                                  (trs_lens[row.Reference] / 1000.0)))

    out_data = OrderedDict([('Reference', refs), ('Count', rpkms)])
    out_df = pd.DataFrame(out_data)
    out_df.to_csv(args.output_count[0], sep="\t", index=False)
