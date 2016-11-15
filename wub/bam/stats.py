# -*- coding: utf-8 -*-

from collections import defaultdict

from wub.bam import common as bam_common
from wub.util import seq as seq_util


def _update_read_stats(r, res, min_aqual):
    if r.is_unmapped:
        res['unmapped'] += 1
        res['unaligned_quals'].append(seq_util.mean_qscore(r.query_qualities))
        res['unaligned_lengths'].append(r.infer_query_length())
    elif r.mapping_quality >= min_aqual:
        res['mapped'] += 1
        res['aligned_quals'].append(seq_util.mean_qscore(r.query_qualities))
        res['alignment_lengths'].append(r.query_alignment_length)
        res['aligned_lengths'].append(r.infer_query_length())
    else:
        res['mapped'] += 1
        res['mqfail_aligned_quals'].append(seq_util.mean_qscore(r.query_qualities))
        res['mqfail_alignment_lengths'].append(r.query_alignment_length)
        res['mqfail_aligned_lengths'].append(r.infer_query_length())


def read_stats(bam, min_aqual=0, region=None):
    """ Parse reads in BAM file and record various statistics. """
    res = defaultdict(list)
    res['unmapped'] = 0
    res['mapped'] = 0
    bam_reader = bam_common.pysam_open(bam, in_format='BAM')
    ue = True
    if region is not None:
        ue = False
    bam_iter = bam_reader.fetch(region=region)
    for r in bam_iter:
        _update_read_stats(r, res, min_aqual)

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


def _register_event(e, query, ref, qpos, rpos, etype, context_sizes, fixed_context=None):
    start = rpos - context_sizes[0]
    end = rpos + context_sizes[1] + 1
    context = str(ref[start:end])

    if start < 0 or end > len(ref):
        return
    if etype == 'match':
        base_to = query[qpos]
    elif etype == 'deletion':
        base_to = '-'
        context = fixed_context
    elif etype == 'insertion':
        base_to = '*'
    else:
        raise Exception('Unsupported event type!')

    e[context][base_to] += 1


def _register_insert(insert, rpos, insertion_lengths, insertion_composition):
    insertion_lengths[len(insert)] += 1
    for base, count in seq_util.base_composition(insert).iteritems():
        insertion_composition[base] += count


def _register_deletion(deletion, match_pos, context_sizes, ref, events, deletion_lengths, r, t):
    if deletion[0] > 0:
        right_contex_end = match_pos + context_sizes[1] + 1
        if right_contex_end <= len(ref):
            _register_event(events, query=r.query_sequence, ref=ref.seq, qpos=t[
                0], rpos=t[0], etype='deletion', context_sizes=context_sizes, fixed_context=deletion[1] + str(ref.seq[match_pos:right_contex_end]))
        deletion_lengths[deletion[0]] += 1
        deletion[0] = 0
        deletion[1] = None


def _update_events(r, ref, events, indel_dists, context_sizes):
    match_pos = 0
    insert = ''
    deletion = [0, None]

    for t in r.get_aligned_pairs():
        if t[0] is None:
            # deletion

            if deletion[0] == 0:
                deletion[1] = str(ref.seq[t[1] - context_sizes[0]:t[1]])
            deletion[0] += 1
            match_pos = t[1]  # FIXME: is this the right coordinate?

            if insert != '':
                _register_event(events, query=r.query_sequence, ref=ref.seq, qpos=t[
                                0], rpos=match_pos, etype='insertion', context_sizes=context_sizes)
                _register_insert(
                    insert, match_pos, indel_dists['insertion_lengths'], indel_dists['insertion_composition'])
                insert = ''

        elif t[1] is None:
            # insertion
            insert += r.query_sequence[t[0]]
            _register_deletion(
                deletion, match_pos, context_sizes, ref, events, indel_dists['deletion_lengths'], r, t)

        else:
            # match or mismatch
            _register_event(events, query=r.query_sequence, ref=ref.seq, qpos=t[
                            0], rpos=t[1], etype='match', context_sizes=context_sizes)
            match_pos = t[1]
            if insert != '':
                _register_event(events, query=r.query_sequence, ref=ref.seq, qpos=t[
                                0], rpos=match_pos, etype='insertion', context_sizes=context_sizes)
                _register_insert(
                    insert, match_pos, indel_dists['insertion_lengths'], indel_dists['insertion_composition'])
                insert = ''
            _register_deletion(
                deletion, match_pos, context_sizes, ref, events, indel_dists['deletion_lengths'], r, t)


def error_and_read_stats(bam, refs, context_sizes=(1, 1), region=None, min_aqual=0):
    """WARNING: context overstepping start/end boundaries are not registered."""
    events = defaultdict(lambda: defaultdict(int))
    read_stats = defaultdict(list)
    read_stats['mapped'] = 0
    read_stats['unmapped'] = 0
    indel_dists = {'insertion_lengths': defaultdict(int), 'deletion_lengths': defaultdict(
        int), 'insertion_composition': defaultdict(int)}

    bam_reader = bam_common.pysam_open(bam, in_format='BAM')

    for r in bam_reader.fetch(region=region, until_eof=True):
        _update_read_stats(r, read_stats, min_aqual)
        if r.is_unmapped:
            continue
        if r.mapping_quality < min_aqual:
            continue
        ref = refs[r.reference_name]
        _update_events(r, ref, events, indel_dists, context_sizes)

    return {'events': events, 'read_stats': read_stats, 'indel_dists': indel_dists}
