#!/usr/bin/env python
__author__ = 'prughani'

import argparse

from Bio import SeqIO
import sys, os, math
import pandas as pd


# test = '/nfs/us-home/DATA/share/prughani/assemblies/arabidopsis/ref/cat_arabidopsis_thaliana.fasta'
# test = '/home/OXFORDNANOLABS/prughani/md8/MD_8_Reads_7_28_16_basecalled_gt8_filtered.fastq'

def _getextention(fast):
    '''
    finds and check for the correct extension

    :param fast: fastq or fasta file
    :return: "fastq" or "fasta"
    '''

    extension = os.path.splitext(os.path.basename(fast))[1]
    if extension in ('.fa', '.fasta'):
        extension = "fasta"
    elif extension in ('.fq', '.fastq'):
        extension = "fastq"
    else:
        raise Exception('Incorrect file format')
        exit()
        # print >> sys.stderr, "Incorrect file format"
    return extension

def readfast(fast):
    '''
    reads a fasta or fastq file

    :param fast: fastq or fasta
    :return: list of records with attr
    '''

    extension = _getextention(fast)
    for rec in SeqIO.parse(open(fast), extension):

        yield rec

def nanoQ(rec):
    '''
    recalculates the mean Q score from the phred quality

    :param rec: fastq records
    :return: mean Qscore rounded to 2dp
    '''
    probs = 0
    for q in rec.letter_annotations["phred_quality"]:
        e = float(math.pow(10.0, -1 * (float(q) / 10.0)))
        probs += e
    av_prob = float(probs / len(rec.seq))
    av_q = float(-10.0 * (math.log10(av_prob)))
    qmean = (round(av_q, 2))
    return qmean



def N50v2(df, col, percent=50):
    '''
    Calculate the N50 by default however, by changing percent to 75 N75 can be calculated

    :param df: dataframe with seqlen column
    :param col: column with sequence length
    :param percent: percentage to be calculated
    :return:
    :rtype: str
    '''

    df1 = df.copy()
    df1 = df1.sort_values(col, ascending=False).reset_index(drop=True)
    df1['cumsum'] = df1[col].cumsum()
    n50 = df1['cumsum'].max() * percent / 100

    ## need to get the mean if n % 2 == 0

    #     if n % 2 == 0:
    #         val = df1.where(df1['cumsum']>=n)[col].dropna()[0:2]
    #         print val.mean()

    #     else:
    #         val = df1.where(df1['cumsum']>=n)[col].dropna().head(1).tolist()[0]
    return df1.where(df1['cumsum'] >= n50)[col].dropna().head(1).tolist()[0]





def GC_per_read(seq_rec, fq=False):



    d = []

    bases=["A", "T", "C", "G", 'N']
    # total_lengths = 0
    for rec in seq_rec:
        tmp = {"SeqID": rec.id, "Seqlen": len(rec.seq), "A": 0, "T": 0, "G": 0, "C": 0, "N": 0}
        for base in bases:
            tmp[base] += rec.seq.count(base)

        if fq:
            tmp['mean_q'] = nanoQ(rec)

        d.append(tmp)

    raw  = pd.DataFrame(d).set_index('SeqID')
    raw['GC content (%)'] = raw.apply(lambda x :float((x['G']) + x['C']) / x['Seqlen'] * 100.0, axis=1)
    for base in bases:
        raw[base+' (%)'] = (raw[base] / raw["Seqlen"]) * 100.0
    raw["other base"] = raw['Seqlen'] - raw[bases].sum(axis=1)
    return raw


def get_stats(df):

    stats = pd.Series({})
    df = df.copy()
    bases = ["A", "T", "C", "G", 'N']

    total_len = int(df["Seqlen"].sum())
    total_bases = df[bases].sum().sum()


    stats['N75'] = N50v2(df, 'Seqlen', 75)
    stats['N50'] = N50v2(df, 'Seqlen', 50)
    stats['N25'] = N50v2(df, 'Seqlen', 25)


    stats['Max contig'] = df['Seqlen'].max()
    stats['Min contig'] = df['Seqlen'].min()    -
    stats['Avg length'] = df['Seqlen'].mean()
    stats['length SD'] = df['Seqlen'].std()

    stats['total bases (Mb)'] = total_len / 100000.0
    stats['other bases (Mb)'] = (total_len - total_bases)/100000.0
    stats['No. contigs'] = df['Seqlen'].count()

    stats["greater then 10 Kb"] = df[df['Seqlen']>=10000.0].Seqlen.count()
    stats["greater then 100 Kb"] = df[df['Seqlen'] >= 100000.0].Seqlen.count()
    stats["greater then 500 Kb"] = df[df['Seqlen'] >= 500000.0].Seqlen.count()
    stats["greater then 1 Mb"] = df[df['Seqlen'] >= 1000000.0].Seqlen.count()

    if 'mean_q' in df.columns:
        stats['Max Qscore'] = df['mean_q'].max()
        stats['Min Qscore'] = df['mean_q'].min()
        stats['Avg Qscore'] = df['mean_q'].mean()
        stats['Qscore SD'] = df['mean_q'].std()


    stats["GC content"] = float(df[['G',"C"]].sum().sum())/ total_len * 100.0
    for base in bases:

        stats[base + ' (%)'] =  float(df[base].sum()) / total_len * 100.0


    return stats.round(2)

