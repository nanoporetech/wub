#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import os
import sys
from collections import defaultdict, OrderedDict
from wub.vis import report
import seaborn as sns
from scipy.stats import spearmanr
import numpy as np

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Scatter plot of two set of counts.')
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF.", required=False, default="plot_counts_correlation.pdf")
parser.add_argument(
    '-T', metavar='tags', type=str, help="Data tags: tag1,tag2.", required=False, default=None)
parser.add_argument(
    '-t', metavar='merged_data', type=str, help="Merged data TSV.", required=False, default=None)
parser.add_argument(
    '-o', metavar='Correlation_tsv', type=str, help="Correlation TSV.", required=False, default=None)
parser.add_argument(
    'counts_one', metavar='counts_one', type=str, help="Input tab separated file.")
parser.add_argument(
    'counts_two', metavar='counts_two', type=str, help="Input tab separated file.")


def _create_tagged_column(df, tag):
    df[tag] = df["Count"]
    df = df.drop("Count", axis=1)
    return df


if __name__ == '__main__':
    args = parser.parse_args()

    data_one = pd.read_csv(args.counts_one, sep="\t")
    data_two = pd.read_csv(args.counts_two, sep="\t")

    # Set data tags:
    tags = args.T
    if tags is not None:
        tags = args.T.split(",")
    else:
        t1 = os.path.basename(args.counts_one).rsplit(".", 1)[0]
        t2 = os.path.basename(args.counts_two).rsplit(".", 1)[0]
        tags = [t1, t2]
    
    # Set column names:
    data_one = _create_tagged_column(data_one, tags[0])
    data_two = _create_tagged_column(data_two, tags[1])

    data_merged = pd.merge(data_one, data_two, on=["Reference"], how="outer")
    data_merged = data_merged.fillna(0.0)

    plotter = report.Report(args.r)

    g = sns.jointplot(tags[0], tags[1], data=data_merged, stat_func=spearmanr, kind="reg")
    plt.tight_layout()
    plotter.pages.savefig()

    plotter.close()

    if args.t is not None:
        data_merged.to_csv(args.t, sep="\t", index=False)

    if args.o is not None:
        rho, pval = spearmanr(data_merged[tags[0]], data_merged[tags[1]])
        res = pd.DataFrame(OrderedDict([("rho", [rho]), ("pval", [pval])]))
        res.to_csv(args.o, sep="\t", index=False)
