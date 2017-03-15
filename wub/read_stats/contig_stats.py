from Bio import SeqIO
import pandas as pd
import numpy as np
from wub.util.misc import _getextension
from wub.util.seq import mean_qscore

# test = '/nfs/us-home/DATA/share/prughani/assemblies/arabidopsis/ref/cat_arabidopsis_thaliana.fasta'
# test = '/home/OXFORDNANOLABS/prughani/md8/MD_8_Reads_7_28_16_basecalled_gt8_filtered.fastq'


def readfast(fast):
    '''
    reads a fasta or fastq file

    :param fast: fastq or fasta
    :return: list of records with attr
    '''

    extension = _getextension(fast)
    for rec in SeqIO.parse(open(fast), extension):

        yield rec


def N50(df, col, percent=50):
    '''
    Calculate the N50 by default however, by changing percent to 75, N75 can be calculated

    :param df: dataframe with seqlen column
    :param col: column with sequence length
    :param percent: percentage to be calculated
    :return: N50 Value
    :rtype: int
    '''

    csum = np.array(df[col], dtype=int)
    csum.sort()
    csum = csum[::-1]
    csum = csum.cumsum()
    n50 = csum.max() * percent / 100
    # result = df.where(df['cumsum'] >= n50)[col].dropna().head(1).tolist()[0]
    for i, cs in enumerate(csum):
        if cs >= n50:
            break
    return df[col][len(df[col]) - i]


def L50(df, col, percent=50):
    '''
    Calculate the L50 by default however, by changing percent to 75, N75 can be calculated

    :param df: dataframe with seqlen column
    :param col: column with sequence length
    :param percent: percentage to be calculated
    :return: N50 Value
    :rtype: int
    '''

    df1 = _cumsum(df, col).copy()
    return df1[df1 >= N50(df, col, percent)][col].count()


def GC_per_read(seq_rec, fq=False):
    ''' Calculates the number of bases per sequence, GC content and mean Q score if fastq is given

    :param seq_rec: sequence records with attr from biopython
    :param fq: boolean
    :return: dataframe
    :rtype: dataframe
    '''

    d = []

    bases = ["A", "T", "C", "G", 'N']
    # total_lengths = 0
    for rec in seq_rec:
        tmp = {"SeqID": rec.id, "Seqlen": len(rec.seq), "A": 0, "T": 0, "G": 0, "C": 0, "N": 0}
        for base in bases:
            tmp[base] += rec.seq.count(base)

        if fq:
            tmp['mean_q'] = round(mean_qscore(rec.letter_annotations[
                                  "phred_quality"], qround=False), 2)

        d.append(tmp)

    raw = pd.DataFrame(d).set_index('SeqID')
    raw['GC content (%)'] = raw.apply(lambda x: float(
        (x['G']) + x['C']) / x['Seqlen'] * 100.0, axis=1)

    for base in bases:
        raw[base + ' (%)'] = (raw[base] / raw["Seqlen"]) * 100.0
    raw["other base"] = raw['Seqlen'] - raw[bases].sum(axis=1)
    return raw


def get_stats(df):
    '''
    calcualtes the summary stats

    :param df: dataframe from GC_per_read
    :return: summary Series
    :rtype: Series
    '''

    stats = pd.Series({})
    df = df.copy()

    Mbase = 1000000.0

    bases = ["A", "T", "C", "G", 'N']

    total_len = int(df["Seqlen"].sum())
    total_bases = df[bases].sum().sum()

    stats['N75'] = N50(df, 'Seqlen', 75)
    stats['N50'] = N50(df, 'Seqlen', 50)
    stats['N25'] = N50(df, 'Seqlen', 25)

    stats['L75'] = L50(df, "Seqlen", 75)
    stats['L50'] = L50(df, "Seqlen", 50)
    stats['L25'] = L50(df, "Seqlen", 25)

    stats['Max contig'] = df['Seqlen'].max()
    stats['Min contig'] = df['Seqlen'].min()
    stats['Avg length'] = df['Seqlen'].mean()
    stats['Length SD'] = df['Seqlen'].std()
    stats['Total length (Mb)'] = total_len / Mbase

    stats['Total bases (Mb)'] = total_len / Mbase
    stats['Other bases (Mb)'] = (total_len - total_bases) / Mbase
    stats['No. contigs'] = df['Seqlen'].count()

    stats["Greater then 10 Kb"] = df[df['Seqlen'] >= 10000.0].Seqlen.count()
    stats["Greater then 100 Kb"] = df[df['Seqlen'] >= 100000.0].Seqlen.count()
    stats["Greater then 500 Kb"] = df[df['Seqlen'] >= 500000.0].Seqlen.count()
    stats["Greater then 1 Mb"] = df[df['Seqlen'] >= 1000000.0].Seqlen.count()

    stats['Yield > 10kb (Mb)'] = df[df['Seqlen'] >= 10000.0]['Seqlen'].sum() / Mbase
    stats['Yield > 50kb (Mb)'] = df[df['Seqlen'] >= 50000.0]['Seqlen'].sum() / Mbase

    if 'mean_q' in df.columns:
        stats['Max Qscore'] = df['mean_q'].max()
        stats['Min Qscore'] = df['mean_q'].min()
        stats['Avg Qscore'] = df['mean_q'].mean()
        stats['Qscore SD'] = df['mean_q'].std()
        stats['Yield >Q6 (Mb)'] = df[df['mean_q'] >= 6.0]['Seqlen'].sum() / Mbase
        stats['Yield >Q9 (Mb)'] = df[df['mean_q'] >= 9.0]['Seqlen'].sum() / Mbase

    stats["GC content"] = float(df[['G', "C"]].sum().sum()) / total_len * 100.0
    for base in bases:

        stats[base + ' (%)'] = float(df[base].sum()) / total_len * 100.0

    return stats.round(2)
