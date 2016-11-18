# -*- coding: utf-8 -*-

from collections import defaultdict, OrderedDict

from wub.bam import common as bam_common
from wub.util import seq as seq_util


def _update_read_stats(r, res, min_aqual):
    if r.is_unmapped:
        res['unmapped'] += 1
        res['unaligned_quals'].append(seq_util.mean_qscore(r.query_qualities))
        res['unaligned_lengths'].append(r.infer_query_length(False))
    elif r.mapping_quality >= min_aqual:
        res['mapped'] += 1
        res['aligned_quals'].append(seq_util.mean_qscore(r.query_qualities))
        res['alignment_lengths'].append(r.query_alignment_length)
        res['aligned_lengths'].append(r.infer_query_length())
        res['mapping_quals'].append(r.mapping_quality)
    else:
        res['mapped'] += 1
        res['mqfail_aligned_quals'].append(seq_util.mean_qscore(r.query_qualities))
        res['mqfail_alignment_lengths'].append(r.query_alignment_length)
        res['mqfail_aligned_lengths'].append(r.infer_query_length())
        res['mapping_quals'].append(r.mapping_quality)


def read_stats(bam, min_aqual=0, region=None):
    """ Parse reads in BAM file and record various statistics. """
    res = {'unmapped': 0,
           'mapped': 0,
           'unaligned_quals': [],
           'unaligned_lengths': [],
           'aligned_quals': [],
           'alignment_lengths': [],
           'aligned_lengths': [],
           'mqfail_aligned_quals': [],
           'mqfail_alignment_lengths': [],
           'mapping_quals': [],
           }
    base_stats = {'aln_length': 0, 'match': 0, 'mismatch': 0, 'deletion': 0, 'insertion': 0}
    read_stats = OrderedDict([
        ("name", []),
        ("ref", []),
        ("coverage", []),
        ("direction", []),
        ("aln_length", []),
        ("insertion", []),
        ("deletion", []),
        ("mismatch", []),
        ("match", []),
        ("identity", []),
        ("accuracy", [])
    ])

    bam_reader = bam_common.pysam_open(bam, in_format='BAM')
    ue = True
    if region is not None:
        ue = False
    bam_iter = bam_reader.fetch(region=region)
    for r in bam_iter:
        _update_read_stats(r, res, min_aqual)

        bs = stats_from_aligned_read(r)
        if bs is not None:
            for k in base_stats.iterkeys():
                base_stats[k] += bs[k]
            for stat, value in bs.iteritems():
                read_stats[stat].append(value)

    base_stats['identity'] = float(
        base_stats['match']) / (base_stats['match'] + base_stats['mismatch'])
    base_stats['accuracy'] = 1.0 - float(base_stats['insertion'] +
                                         base_stats['deletion'] + base_stats['mismatch']) / base_stats['aln_length']
    res['base_stats'] = base_stats
    res['read_stats'] = read_stats
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
    return {'qualities': dict(st), 'coverage': dict(cst)}


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


def _update_events(r, ref, events, indel_dists, context_sizes, base_stats):
    match_pos = 0
    insert = ''
    deletion = [0, None]

    aligned_pairs = r.get_aligned_pairs()
    # Remove soft clips:
    if r.cigartuples[0][0] == 4:
        aligned_pairs = aligned_pairs[r.cigartuples[0][1]:]
    if r.cigartuples[-1][0] == 4:
        aligned_pairs = aligned_pairs[:-r.cigartuples[-1][1]]

    for t in aligned_pairs:
        if t[0] is None:
            # deletion
            base_stats['deletion'] += 1

            if deletion[0] == 0:
                deletion[1] = str(ref.seq[t[1] - context_sizes[0]:t[1]])
            deletion[0] += 1
            match_pos = t[1]

            if insert != '':
                _register_event(events, query=r.query_sequence, ref=ref.seq, qpos=t[
                                0], rpos=match_pos, etype='insertion', context_sizes=context_sizes)
                _register_insert(
                    insert, match_pos, indel_dists['insertion_lengths'], indel_dists['insertion_composition'])
                insert = ''

        elif t[1] is None:
            # insertion
            base_stats['insertion'] += 1

            insert += r.query_sequence[t[0]]
            _register_deletion(
                deletion, match_pos, context_sizes, ref, events, indel_dists['deletion_lengths'], r, t)

        else:
            # match or mismatch
            if r.query_sequence[t[0]] == ref.seq[t[1]]:
                base_stats['match'] += 1
            else:
                base_stats['mismatch'] += 1

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
    read_stats = {'unmapped': 0,
                  'mapped': 0,
                  'unaligned_quals': [],
                  'unaligned_lengths': [],
                  'aligned_quals': [],
                  'alignment_lengths': [],
                  'aligned_lengths': [],
                  'mqfail_aligned_quals': [],
                  'mqfail_alignment_lengths': [],
                  'mapping_quals': [],
                  }
    indel_dists = {'insertion_lengths': defaultdict(int), 'deletion_lengths': defaultdict(
        int), 'insertion_composition': defaultdict(int)}

    bam_reader = bam_common.pysam_open(bam, in_format='BAM')
    base_stats = {'match': 0, 'mismatch': 0, 'deletion': 0, 'insertion': 0}

    for r in bam_reader.fetch(region=region, until_eof=True):
        _update_read_stats(r, read_stats, min_aqual)
        if r.is_unmapped:
            continue
        if r.mapping_quality < min_aqual:
            continue
        ref = refs[r.reference_name]
        _update_events(r, ref, events, indel_dists, context_sizes, base_stats)

    base_stats['aln_length'] = base_stats['match'] + base_stats['mismatch'] + \
        base_stats['insertion'] + base_stats['deletion']
    base_stats['identity'] = float(
        base_stats['match']) / (base_stats['match'] + base_stats['mismatch'])
    base_stats['accuracy'] = 1.0 - \
        float(base_stats['mismatch'] + base_stats['insertion'] + base_stats['deletion']) / \
        base_stats['aln_length']

    res = {'events': dict(events), 'read_stats': dict(
        read_stats), 'indel_dists': dict(indel_dists), 'base_stats': base_stats}
    return res


def stats_from_aligned_read(read):
    """Create summary information for an aligned read (taken from tang.util.bio).

    :param read: :class:`pysam.AlignedSegment` object
    """
    tags = dict(read.tags)
    try:
        tags.get('NM')
    except:
        raise IOError(
            "Read is missing required 'NM' tag. Try running 'samtools fillmd -S - ref.fa'.")
    name = read.qname
    if read.flag == 4:
        return None
    match = reduce(lambda x, y: x + y[1] if y[0] == 0 else x, read.cigar, 0)
    ins = reduce(lambda x, y: x + y[1] if y[0] == 1 else x, read.cigar, 0)
    delt = reduce(lambda x, y: x + y[1] if y[0] == 2 else x, read.cigar, 0)
    # NM is edit distance: NM = INS + DEL + SUB
    sub = tags['NM'] - ins - delt
    length = match + ins + delt
    iden = 100 * float(match - sub) / match
    acc = 100 - 100 * float(tags['NM']) / length
    coverage = 100 * float(read.query_alignment_length) / read.infer_query_length()
    direction = '-' if read.is_reverse else '+'
    results = OrderedDict([
        ("name", name),
        ("ref", read.reference_name),
        ("coverage", coverage),
        ("direction", direction),
        ("aln_length", length),
        ("insertion", ins),
        ("deletion", delt),
        ("mismatch", sub),
        ("match", match - sub),
        ("identity", iden),
        ("accuracy", acc)
    ])
    return results