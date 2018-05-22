# -*- coding: utf-8 -*-

import six
from six.moves import reduce
import sys
from collections import defaultdict, OrderedDict
import tqdm
import numpy as np

from wub.bam import common as bam_common
from wub.util import seq as seq_util


def _frag_dict_to_array(fd, chrom_lengths):
    """ Convert fragment coverage dictionary into a numpy array. """
    res = {}
    for ref, frags in six.iteritems(fd):
        cov = np.zeros((chrom_lengths[ref],), dtype=int)
        for pos, count in six.iteritems(frags):
            cov[pos] = count
        res[ref] = cov
    return res


def frag_coverage(bam, chrom_lengths, region=None, min_aqual=0, ref_cov=True, verbose=True):
    """ Calculate fragment coverage vectors on the forward and reverse strands.

    :param bam: Input bam file.
    :param chrom_lengths: Dictionary of chromosome names and lengths.
    :param region: Restrict parsing to the specified region.
    :param min_aqual: Minimum mapping quality.
    :param verbose: Display progress bar.
    :returns: Forward and reverse fragment coverage vectors.
    :rtype: dict
    """

    frags_fwd = defaultdict(lambda: defaultdict(int))
    frags_rev = defaultdict(lambda: defaultdict(int))

    aln_ref_cov = (defaultdict(list))

    bam_reader = bam_common.pysam_open(bam, in_format='BAM')
    ue = True
    if region is not None:
        ue = False
    bam_iter = bam_reader.fetch(region=region, until_eof=ue)

    try:
        total_reads = bam_reader.mapped + bam_reader.unmapped
    except:
        total_reads = None
    if verbose and region is None:
        sys.stdout.write(
            "Gathering fragment statistics from file: {}\n".format(bam))
        bam_iter = tqdm.tqdm(bam_iter, total=total_reads)

    for r in bam_iter:
        # Skip unmapped reads:
        if r.is_unmapped:
            continue
        # Skip if mapping quality is too low:
        if r.mapq < min_aqual:
            continue
        pos = r.reference_start
        ref = r.reference_name
        if r.is_reverse:
            frags_rev[r.reference_name][pos] += 1
        else:
            frags_fwd[r.reference_name][pos] += 1

        if ref_cov:
            aln_ref_cov[ref].append(r.reference_length / float(chrom_lengths[ref]))

    frags_fwd = _frag_dict_to_array(frags_fwd, chrom_lengths)
    frags_rev = _frag_dict_to_array(frags_rev, chrom_lengths)

    res = {'frags_fwd': frags_fwd, 'frags_rev': frags_rev, 'ref_cov': aln_ref_cov}
    return res


def _update_read_stats(r, res, min_aqual):
    """ Update read statistics. """
    if r.is_unmapped:
        res['unmapped'] += 1
        res['unaligned_lengths'].append(r.query_length)

        if r.query_qualities is not None:
            res['unaligned_quals'].append(
                seq_util.mean_qscore(r.query_qualities, qround=False))
    elif r.mapping_quality >= min_aqual:
        res['mapped'] += 1
        if r.query_qualities is not None:
            res['aligned_quals'].append(
                seq_util.mean_qscore(r.query_qualities, qround=False))
        res['alignment_lengths'].append(r.query_alignment_length)
        res['aligned_lengths'].append(r.infer_query_length())
        res['mapping_quals'].append(r.mapping_quality)
    else:
        res['mapped'] += 1
        if r.query_qualities is not None:
            res['mqfail_aligned_quals'].append(
                seq_util.mean_qscore(r.query_qualities, qround=False))
        res['mqfail_alignment_lengths'].append(r.query_alignment_length)
        res['mqfail_aligned_lengths'].append(r.infer_query_length())
        res['mapping_quals'].append(r.mapping_quality)


def read_stats(bam, min_aqual=0, region=None, with_clipps=False, verbose=True):
    """ Parse reads in BAM file and record various statistics.

    :param bam: BAM file.
    :param min_aqual: Minimum mapping quality, skip read if mapping quality is lower.
    :param region: smatools region.
    :param with_clipps: Take into account clipps when calculating accuracy.
    :param verbose: Show progress bar.
    :returns: A dictionary with various global and per-read statistics.
    :rtype: dict
    """
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
    base_stats = {'aln_length': 0, 'match': 0, 'mismatch': 0,
                  'deletion': 0, 'insertion': 0, 'clipps': 0}
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
        ("accuracy", []),
        ("clipps", [])
    ])

    bam_reader = bam_common.pysam_open(bam, in_format='BAM')
    ue = True
    if region is not None:
        ue = False
    bam_iter = bam_reader.fetch(region=region, until_eof=ue)

    try:
        total_reads = bam_reader.mapped + bam_reader.unmapped
    except:
        total_reads = None
    if verbose and region is None:
        sys.stdout.write(
            "Gathering read statistics from file: {}\n".format(bam))
        bam_iter = tqdm.tqdm(bam_iter, total=total_reads)

    for r in bam_iter:
        # Update basic read statistics:
        _update_read_stats(r, res, min_aqual)

        # Get detailed statistics from aligned read and
        # updated global stats:
        bs = stats_from_aligned_read(r, with_clipps)
        # bs is None for unaligned reads.
        if bs is not None:
            for k in six.iterkeys(base_stats):
                base_stats[k] += bs[k]
            for stat, value in six.iteritems(bs):
                read_stats[stat].append(value)

    # Calculate global identity and accuracy:
    base_stats['identity'] = float(
        base_stats['match']) / (base_stats['match'] + base_stats['mismatch'])

    clipps = 0
    if with_clipps:
        clipps = base_stats['clipps']

    base_stats['accuracy'] = 1.0 - (float(base_stats['insertion'] +
                                          base_stats['deletion'] + base_stats['mismatch'] + clipps) / base_stats['aln_length'])
    res['base_stats'] = base_stats
    res['read_stats'] = read_stats
    bam_reader.close()
    return res


def pileup_stats(bam, region=None, verbose=True, with_quals=True):
    """ Parse pileup columns and extract quality values.

    :param bam: Input BAM file.
    :param region: samtools region.
    :param verbose: Show progress bar.
    :param with_quals: Return quality values per position.
    :returns: Dictionaries per reference with per-base coverage and quality values.
    :rtype: dict
    """
    st = defaultdict(lambda: defaultdict(list))
    cst = defaultdict(lambda: defaultdict(int))
    samfile = bam_common.pysam_open(bam, in_format='BAM')

    pileup_iter = samfile.pileup(region=region)
    start, end = None, None
    if region is not None:
        tmp = region.split(":")
        _, start, end = tmp[0], int(tmp[1]) - 1, int(tmp[2])
    if verbose:
        sys.stdout.write(
            "Gathering pileup statistics from file: {}\n".format(bam))
        total_bases = sum(samfile.lengths)
        if region is not None:
            tmp = region.split(":")
            total_bases = int(tmp[2]) - int(tmp[1])
        pileup_iter = tqdm.tqdm(pileup_iter, total=total_bases)

    for pileupcolumn in pileup_iter:
        if region is not None and (pileupcolumn.reference_pos < start or pileupcolumn.reference_pos >= end):
            continue
        # print pileupcolumn.reference_name, pileupcolumn.reference_pos,
        # pileupcolumn.nsegments
        cst[pileupcolumn.reference_name][
            pileupcolumn.reference_pos] = pileupcolumn.nsegments
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # print pileupcolumn.reference_name, pileupcolumn.reference_pos,
                # pileupread.alignment.query_qualities[pileupread.query_position]
                if (pileupread.alignment.query_qualities is not None) and with_quals:
                    st[pileupcolumn.reference_name][pileupcolumn.reference_pos].append(
                        pileupread.alignment.query_qualities[pileupread.query_position])
    samfile.close()
    return {'qualities': dict(st), 'coverage': dict(cst)}


def _register_event(e, query, ref, qpos, rpos, etype, context_sizes, fixed_context=None):
    """Register event withing a context."""
    start = rpos - context_sizes[0]
    end = rpos + context_sizes[1] + 1
    context = str(ref[start:end])

    # If context goes outside the reference boundaries,
    # do not register it:
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
    """Register insertion length and base composition."""
    insertion_lengths[len(insert)] += 1
    for base, count in six.iteritems(seq_util.base_composition(insert)):
        insertion_composition[base] += count


def _register_deletion(deletion, match_pos, context_sizes, ref, events, deletion_lengths, r, t):
    """Register deletion and deletion lengths."""
    if deletion[0] > 0:  # We have an active deletion if length is larger than zero.
        right_contex_end = match_pos + context_sizes[1] + 1
        if right_contex_end <= len(ref):
            # Deletions are registered with fixed context as we keep trakc of the left context
            # in the deletion tuple:
            _register_event(events, query=r.query_sequence, ref=ref.seq, qpos=t[
                0], rpos=t[0], etype='deletion', context_sizes=context_sizes, fixed_context=deletion[1] + str(ref.seq[match_pos:right_contex_end]))
        deletion_lengths[deletion[0]] += 1
        # Reset deletion size to zero and context to None: we left deletion
        # event.
        deletion[0] = 0
        deletion[1] = None


def _update_events(r, ref, events, indel_dists, context_sizes, base_stats):
    """Update event details."""
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
            # Increse delelted base count:
            base_stats['deletion'] += 1

            if deletion[0] == 0:
                # We are at the left edge of the deletion, register context:
                deletion[1] = str(ref.seq[t[1] - context_sizes[0]:t[1]])
            deletion[0] += 1
            # Set match position to the last seen reference base:
            match_pos = t[1]

            # Register insert if we are switching from deleltion to insertion:
            if insert != '':
                _register_event(events, query=r.query_sequence, ref=ref.seq, qpos=t[
                                0], rpos=match_pos, etype='insertion', context_sizes=context_sizes)
                _register_insert(
                    insert, match_pos, indel_dists['insertion_lengths'], indel_dists['insertion_composition'])
                insert = ''

        elif t[1] is None:
            # insertion
            # Increase inserted base count:
            base_stats['insertion'] += 1

            # Append base to insert:
            insert += r.query_sequence[t[0]]
            # If switching from deleltion to insertion:
            _register_deletion(
                deletion, match_pos, context_sizes, ref, events, indel_dists['deletion_lengths'], r, t)

        else:
            # match or mismatch
            # Increase mismatch and mathc base counts:
            if r.query_sequence[t[0]] == ref.seq[t[1]]:
                base_stats['match'] += 1
            else:
                base_stats['mismatch'] += 1

            _register_event(events, query=r.query_sequence, ref=ref.seq, qpos=t[
                            0], rpos=t[1], etype='match', context_sizes=context_sizes)
            match_pos = t[1]
            # We are switching from and insertion:
            if insert != '':
                _register_event(events, query=r.query_sequence, ref=ref.seq, qpos=t[
                                0], rpos=match_pos, etype='insertion', context_sizes=context_sizes)
                _register_insert(
                    insert, match_pos, indel_dists['insertion_lengths'], indel_dists['insertion_composition'])
                insert = ''
            # We are switching from a deletion:
            _register_deletion(
                deletion, match_pos, context_sizes, ref, events, indel_dists['deletion_lengths'], r, t)


def error_and_read_stats(bam, refs, context_sizes=(1, 1), region=None, min_aqual=0, verbose=True):
    """Gather read statistics and context-dependend error statistics from BAM file.
    WARNING: context overstepping reference start/end boundaries are not registered.

    Definition of context: for substitutions the event is happening from the "central base", in the case of indels the events are located
    between the central base and the base before.

    :param bam: Input BAM file.
    :param refs: Dictionary of references.
    :param context_sizes: The size of the left and right contexts.
    :param region: samtools regions.
    :param min_qual: Minimum mappign quality.
    :param verbose: Show progress bar.
    :returns: Dictionary with read and error statistics.
    :rtype: dict
    """
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

    read_iter = bam_reader.fetch(region=region, until_eof=True)
    if verbose:
        sys.stdout.write(
            "Gathering read and error statistics from file: {}\n".format(bam))
        try:
            total_reads = bam_reader.mapped + bam_reader.unmapped
        except:
            total_reads = None
        read_iter = tqdm.tqdm(read_iter, total=total_reads)

    for r in read_iter:
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


def stats_from_aligned_read(read, with_clipps=False):
    """Create summary information for an aligned read (modified from tang.util.bio).

    :param read: :class:`pysam.AlignedSegment` object
    :param with_clipps:
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

    # Count clips:
    clipps = reduce(
        lambda x, y: x + y[1] if (y[0] == 4 or y[0] == 5) else x, read.cigar, 0)
    if with_clipps:
        length += clipps

    iden = float(match - sub) / match
    if with_clipps:
        acc = 1.0 - (float(tags['NM'] + clipps) / length)
    else:
        acc = 1.0 - (float(tags['NM']) / length)
    coverage = float(read.query_alignment_length) / read.infer_query_length()
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
        ("accuracy", acc),
        ("clipps", clipps),
    ])
    return results
