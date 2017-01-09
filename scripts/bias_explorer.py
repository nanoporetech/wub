#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import numpy as np

import pandas as pd
from wub.vis import report
from wub.util import seq as seq_util
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
import itertools


import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    # import seaborn as sns
warnings.resetwarnings()

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Explore the effects of length and GC content on transcript coverage using multivariate regression.

    This script takes the output of bam_count_reads.py (with the -z option) and uses a NegativeBinomial GLM to fit transcript
    length and GC content (quadratic fit) on the counts. If the true transcript proportions (Target) are inculded using the -t option
    then they are included in the model. I the -b switch is not specified then the dinucleotide frequnecies are alos included as predictors.

    Then a quadratic fit is performed with the GC content and Length against the counts individually.
    Various plots are produced in the report PDF.
    """)
parser.add_argument(
    '-t', metavar='target_tsv', type=str, help="Tab separated file containing the target counts (\"true concentrations\") (None).", default=None)
parser.add_argument(
    '-b', action='store_true', help="Do not include dinucleotide frequencies as predictors (False).", default=False)
parser.add_argument(
    '-q', action='store_true', help="Include square of dinucleotide frequencies as predictors (False).", default=False)
parser.add_argument(
    '-l', action='store_true', help="Label points and make plots ugly (False).", default=False)
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF (bias_explorer.pdf).", default='bias_explorer.pdf')
parser.add_argument(
    '-m', metavar='merged_tsv', type=str, help="Save merged data frame in this file (None).", default=None)
parser.add_argument(
    'counts', metavar='counts', type=str, help="Tab separated file with counts and features. Produced by bam_count_read.py with the -z option.")


def label_points(target, md, predicted, fontsize=6):
    for ref, x, y in zip(md["Reference"], md[target], predicted):
        plotter.plt.text(x, y, ref, fontsize=fontsize)


def global_model(md, with_target, label=False):
    """Fit all predictors on counts.
    :param md: Input data frame.
    :param with_target: Include Target if true.
    """
    md = md.sort(['Count'])
    if with_target:
        formula = 'Count ~ Target + Length + Length2 + GC_content + GC_content2'
    else:
        formula = 'Count ~ Length + Length2 + GC_content + GC_content2'

    # Add dinucleotide frequencies:
    if not args.b:
        for kmer in itertools.product(*([seq_util.bases] * 2)):
            kmer = ''.join(kmer)
            if args.q:
                formula += " + {} + {}2".format(kmer, kmer)
            else:
                formula += " + {}".format(kmer)

    print "\nFitting: ", formula, "\n"
    res = smf.glm(
        formula=formula, data=md, family=sm.families.NegativeBinomial()).fit()
    print res.summary()
    print "Null deviance: ", res.null_deviance, "Null deviance/Deviance: ", res.null_deviance / res.deviance

    if args.m is not None:
        md.to_csv(args.m, sep="\t", index=False)

    slope, intercept, r_value, p_value, std_err = stats.linregress(
        md["Count"], res.fittedvalues)
    plotter.plt.plot(md["Count"], res.fittedvalues, 'o')
    y_values = [slope * i + intercept for i in md["Count"]]
    plotter.plt.plot(
        md["Count"], y_values, 'g-', label="r={:.3f}, p={:.3f}".format(r_value, p_value))
    plotter.plt.legend(loc='best')
    plotter.plt.title("Actual vs. predicted read counts")
    plotter.plt.xlabel("Count")
    plotter.plt.ylabel("Predicted count")
    if label:
        label_points("Count", md, res.fittedvalues)
    plotter.pages.savefig()
    plotter.plt.close()


def gc_model(md, label=False):
    """Quadratic fit of GC content on counts.
    :param md: Input data frame.
    """
    md = md.sort(['GC_content'])
    formula = 'Count ~ GC_content + GC_content2'
    print
    print "\nFitting: ", formula, "\n"
    results = smf.glm(
        formula=formula, data=md, family=sm.families.NegativeBinomial()).fit()
    print results.summary()
    print "Null deviance: ", results.null_deviance, "Null deviance/Deviance: ", results.null_deviance / results.deviance

    plotter.plt.plot(md["GC_content"], md["Count"], 'o', label='data')
    plotter.plt.plot(md["GC_content"], results.fittedvalues, '-', label='Predicted')
    if label:
        label_points("GC_content", md, md["Count"])
    plotter.plt.title("GC content vs. read counts")
    plotter.plt.xlabel("GC content")
    plotter.plt.ylabel("Count")
    plotter.plt.legend(loc='best')
    plotter.pages.savefig()
    plotter.plt.close()


def base_model(md, base, label=False):
    """Quadratic fit of GC content on counts.
    :param md: Input data frame.
    """
    md = md.sort([base])
    formula = 'Count ~ {} + {}2'.format(base, base)
    print
    print "\nFitting: ", formula, "\n"
    results = smf.glm(
        formula=formula, data=md, family=sm.families.NegativeBinomial()).fit()
    print results.summary()
    print "Null deviance: ", results.null_deviance, "Null deviance/Deviance: ", results.null_deviance / results.deviance

    plotter.plt.plot(md[base], md["Count"], 'o', label='data')
    plotter.plt.plot(md[base], results.fittedvalues, '-', label='Predicted')
    if label:
        label_points(base, md, md["Count"])
    plotter.plt.title("{} frequency vs. read counts".format(base))
    plotter.plt.xlabel("{} frequency".format(base))
    plotter.plt.ylabel("Count")
    plotter.plt.legend(loc='best')
    plotter.pages.savefig()
    plotter.plt.close()


def length_model(md, label):
    """Quadratic fit of transcript length on counts.
    :param md: Input data frame.
    """
    md = md.sort(['Length'])
    md = md.sort(['Length'])
    formula = 'Count ~ Length + Length2'
    print "\nFitting: ", formula, "\n"
    results = smf.glm(
        formula=formula, data=md, family=sm.families.NegativeBinomial()).fit()
    print results.summary()
    print "Null deviance: ", results.null_deviance, "Null deviance/Deviance: ", results.null_deviance / results.deviance

    plotter.plt.plot(md["Length"], md["Count"], 'o', label='data')
    plotter.plt.plot(
        md["Length"], results.fittedvalues, '-', label='Predicted')
    if label:
        label_points("Length", md, md["Count"])
    plotter.plt.legend(loc='best')
    plotter.plt.title("Length vs. read counts")
    plotter.plt.xlabel("Length")
    plotter.plt.ylabel("Count")
    plotter.pages.savefig()
    plotter.plt.close()


if __name__ == '__main__':
    args = parser.parse_args()

    plotter = report.Report(args.r)

    # Load data:
    data = pd.read_csv(args.counts, sep="\t")

    # Augment data frame with true proportions if present:
    with_target = False
    if args.t is not None:
        target = pd.read_csv(args.t, sep="\t")
        target['Target'] = target['Count']
        del target['Count']
        md = target.merge(data, how='outer', on='Reference').dropna()
        with_target = True
    else:
        md = data

    # Augment data frame with square dinucleotide frequencies:
    if not args.b:
        for kmer in itertools.product(*([seq_util.bases] * 2)):
            kmer = ''.join(kmer)
            md[kmer + "2"] = md[kmer] ** 2

    # Matrix scatter plot of counts and features:
    fields = ['Count', 'Length', 'GC_content']
    pd.tools.plotting.scatter_matrix(data[fields], diagonal="kde")
    plotter.pages.savefig()
    plotter.plt.close()

    # Add squared features for GC and Length:
    md["GC_content2"] = md["GC_content"] ** 2
    md["Length2"] = md["Length"] ** 2

    global_model(md, with_target, args.l)
    gc_model(md, args.l)
    length_model(md, args.l)
    if not args.b:
        for base in itertools.product(*([seq_util.bases] * 2)):
            base = ''.join(base)
            base_model(md, base, args.l)

    plotter.close()
