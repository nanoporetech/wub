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


import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    # import seaborn as sns
warnings.resetwarnings()

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Explore the effects of length and GC content on transcript coverage using multivariate regression.

    This script takes the output of bam_count_reads.py (with the -z option) and uses a NegativeBinomial GLM to fit transcript
    length and GC contenti (quadratic fit) on the counts. If the true transcript proportions (Target) are inculded using the -t option#
    then they are included in the model.

    Then a quadratic fir is performed with the GC content and Length against the counts individually.
    Various plots are produced in the report PDF.
    """)
parser.add_argument(
    '-t', metavar='target_tsv', type=str, help="Tab separated file containing the target counts (\"true concentrations\") (None).", default=None)
parser.add_argument(
    '-b', action='store_true', help="Do not include individual base frequencies as predictors (False).", default=False)
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF (bias_explorer.pdf).", default='bias_explorer.pdf')
parser.add_argument(
    'counts', metavar='counts', type=str, help="Tab separated file with counts and features. Produced by bam_count_read.py with the -z option.")


def global_model(md, with_target):
    """Fit all predictors on counts.
    :param md: Input data frame.
    :param with_target: Include Target if true.
    """
    md = md.sort(['Count'])
    if with_target:
        formula = 'Count ~ Target + Length + GC + GC2'
    else:
        formula = 'Count ~ Length + GC + GC2'

    if not args.b:
        for base in seq_util.bases:
            formula += " + {} + {}2".format(base, base)

    print "\nFitting: ", formula, "\n"
    res = smf.glm(formula=formula, data=md, family=sm.families.NegativeBinomial()).fit()
    print res.summary()
    print "Null deviance: ", res.null_deviance, "Null deviance/Deviance: ", res.null_deviance / res.deviance

    plotter.plt.plot(md["Count"], res.fittedvalues, 'o')
    plotter.plt.title("Actual vs. predicted read counts")
    plotter.plt.xlabel("Count")
    plotter.plt.ylabel("Predicted count")
    plotter.pages.savefig()
    plotter.plt.close()
    print stats.spearmanr(md["Count"], res.fittedvalues)
    print stats.pearsonr(md["Count"], res.fittedvalues)


def gc_model(md):
    """Quadratic fit of GC content on counts.
    :param md: Input data frame.
    """
    md = md.sort(['GC'])
    formula = 'Count ~ GC + GC2'
    print
    print "\nFitting: ", formula, "\n"
    results = smf.glm(
        formula=formula, data=md, family=sm.families.NegativeBinomial()).fit()
    print results.summary()
    print "Null deviance: ", results.null_deviance, "Null deviance/Deviance: ", results.null_deviance / results.deviance

    plotter.plt.plot(md["GC"], md["Count"], 'o', label='data')
    plotter.plt.plot(md["GC"], results.fittedvalues, '-', label='Predicted')
    plotter.plt.title("GC content vs. read counts")
    plotter.plt.xlabel("GC content")
    plotter.plt.ylabel("Count")
    plotter.plt.legend(loc='best')
    plotter.pages.savefig()
    plotter.plt.close()

def base_model(md, base):
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
    plotter.plt.title("{} frequency vs. read counts".format(base))
    plotter.plt.xlabel("{} frequency".format(base))
    plotter.plt.ylabel("Count")
    plotter.plt.legend(loc='best')
    plotter.pages.savefig()
    plotter.plt.close()



def length_model(md):
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
    plotter.plt.legend(loc='best')
    plotter.plt.title("Length content vs. read counts")
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

    # Augment data frame with square base frequencies:
    if not args.b:
        for base in seq_util.bases:
            md[base + "2"] = md[base] ** 2

    # Matrix scatter plot of counts and features:
    fields = ['Count', 'Length', 'GC']
    pd.tools.plotting.scatter_matrix(data[fields], diagonal="kde")
    plotter.pages.savefig()
    plotter.plt.close()

    md["GC2"] = md["GC"] ** 2
    md["Length2"] = md["Length"] ** 2

    global_model(md, with_target)
    gc_model(md)
    length_model(md)
    if not args.b:
        for base in seq_util.bases:
            base_model(md, base)

    plotter.close()
