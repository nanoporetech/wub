# -*- coding: utf-8 -*-

from collections import defaultdict

from wub.bam import common as bam_common
from wub.util import seq as seq_util


def read_stats(bam, region=None):
    """ Parse reads in BAM file and record various statistics. """
    res = defaultdict(list)
    bam_reader = bam_common.pysam_open(bam, in_format='BAM')
    if region is None:
        bam_iter = bam_reader.fetch(until_eof=True)
    else:
        bam_iter = bam_reader.fetch(region=region, until_eof=True)
    for r in bam_iter:
        if r.is_unmapped:
            res['unaligned_quals'].append(seq_util.mean_qscore(r.query_qualities))
            res['unaligned_lengths'].append(r.infer_query_length())
        else:
            res['aligned_quals'].append(seq_util.mean_qscore(r.query_qualities))
            res['alignment_lengths'].append(r.query_alignment_length)
            res['aligned_lengths'].append(r.infer_query_length())
    bam_reader.close()
    return res


def pileup_stats(bam, region=None):
    """ Parse pileup columns and extract quality values. """
    st = defaultdict(lambda: defaultdict(list))
    cst = defaultdict(lambda: defaultdict(int))
    samfile = bam_common.pysam_open(bam, in_format='BAM')
    for pileupcolumn in samfile.pileup(region=region):
        # print pileupcolumn.reference_name, pileupcolumn.reference_pos, pileupcolumn.nsegments
        cst[pileupcolumn.reference_name][pileupcolumn.reference_pos] = pileupcolumn.nsegments
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # print pileupcolumn.reference_name, pileupcolumn.reference_pos,
                # pileupread.alignment.query_qualities[pileupread.query_position]
                st[pileupcolumn.reference_name][pileupcolumn.reference_pos].append(
                    pileupread.alignment.query_qualities[pileupread.query_position])
    samfile.close()
    return {'qualities': st, 'coverage': cst}


def _register_deletion(e, ref, pos):
    start = pos - 1
    end = pos + 2
    if start < 0 or end > len(ref):
        return
    context = str(ref[start:end])
    e[context]['-'] += 1


def _register_match(e, query, ref, qpos, rpos):
    start = rpos - 1
    end = rpos + 2
    if start < 0 or end > len(ref):
        return
    # if query[qpos] == ref[rpos]:
    #    return
    context = str(ref[start:end])
    e[context][query[qpos]] += 1


def _register_insert(e, ref, rpos):
    start = rpos - 1
    end = rpos + 2
    if start < 0 or end > len(ref):
        return
    context = str(ref[start:end])
    e[context]['*'] += 1


def error_stats(bam, refs, region=None, min_aqual=0):
    events = defaultdict(lambda: defaultdict(int))
    bam_reader = bam_common.pysam_open(bam, in_format='BAM')
    nr_reads = 0
    for r in bam_reader.fetch(region=region, until_eof=True):
        nr_reads += 1
        if r.is_unmapped:
            continue
        if r.mapping_quality < min_aqual:
            continue
        match_pos = 0
        insert = ''
        ref = refs[r.reference_name]
        for t in r.get_aligned_pairs():
            if t[0] is None:
                # deletion
                _register_deletion(events, ref.seq, t[1])
                if insert != '':
                    _register_insert(events, ref.seq, match_pos)
                    insert = ''
                match_pos = t[1]
            elif t[1] is None:
                # insertion
                insert += r.query_sequence[t[0]]
            else:
                # match or mismatch
                _register_match(events, r.query_sequence, ref.seq, t[0], t[1])
                if insert != '':
                    _register_insert(events, ref.seq, match_pos)
                    insert = ''
                match_pos = t[1]
    return events
