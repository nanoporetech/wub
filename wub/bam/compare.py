# -*- coding: utf-8 -*-
"""Compares alignments in two BAM files."""

from itertools import izip, chain
from collections import OrderedDict
from wub.bam import common as bam_common


def is_coarse_match(aln_diff, tolerance):
    """Determine if start and end postions of two alignments are within
    the specified tolerance levels.

    :param aln_diff: Alignment diff structure as returned by compare_alignments.
    :returns: True or False
    :rtype: bool
    """
    if abs(aln_diff['start_pos'][1] - aln_diff['start_pos'][0]) > tolerance:
        return False
    if abs(aln_diff['end_pos'][1] - aln_diff['end_pos'][0]) > tolerance:
        return False
    return True


def bam_compare(aln_one, aln_two, coarse_tolerance=50, strict_flags=False, in_format='BAM'):
    """Count reads mapping to references in a BAM file.

    :param alignment_file: BAM file.
    :param min_aln_qual: Minimum mapping quality.
    :returns: Dictionary with read counts per reference.
    :rtype: dict
    """

    aln_iter_one = bam_common.pysam_open(aln_one, in_format)
    aln_iter_two = bam_common.pysam_open(aln_two, in_format)

    # Comparison summary structure:
    stats = OrderedDict([
        ('BamFiles', [aln_one, aln_two]),
        ('TotalQueries', 0),
        ('DirectionMismatch', 0),
        ('RefMismatch', 0),
        ('StrictFlagMismatch', 0),
        ('SeqMismatch', 0),
        ('CoarseMatches', 0),
        ('CommonAlignedBases', 0),
        ('CommonMatchingBases', 0),
        ('PerQueryBaseSim', []),
        ('PerQueryBaseSimClipped', []),
        (aln_one, {'HardClippedBases': 0, 'SoftClippedBases': 0, 'AlignedBases': 0,
                   'UnalignedQueries': 0, 'AlignedQueries': 0}),
        (aln_two, {'HardClippedBases': 0, 'SoftClippedBases': 0', AlignedBases': 0,
                   'UnalignedQueries': 0, 'AlignedQueries': 0}),
        ('AlignedSimilarity', 0.0),
        ])

    for segments in izip(aln_iter_one.fetch(until_eof=True), aln_iter_two.fetch(until_eof=True)):
        aln_diff = compare_alignments(segments[0], segments[1], strict_flags)
        stats['TotalQueries'] += 1

        # Register hard and soft clipped bases:
        stats[aln_one]['HardClippedBases'] += aln_diff['hard_clipped'][0]
        stats[aln_two]['HardClippedBases'] += aln_diff['hard_clipped'][1]

        stats[aln_one]['SoftClippedBases'] += aln_diff['soft_clipped'][0]
        stats[aln_two]['SoftClippedBases'] += aln_diff['soft_clipped'][1]

        # Both reads are aligned:
        if aln_diff['mapped'] == (True, True):
            stats[aln_one]['AlignedQueries'] += 1
            stats[aln_two]['AlignedQueries'] += 1

            # Reference mismatch:
            if aln_diff['ref_match'] is False:
                stats['RefMismatch'] = + 1
                continue

            # Orientation mismatch:
            if aln_diff['dir_match'] is False:
                stats['DirectionMismatch'] += 1
                continue

            # Flag mismatch:
            if aln_diff['flag_match'] is False:
                stats['StrictFlagMismatch'] += 1
                continue

            # Sequence mismatch:
            if aln_diff['seq_match'] is False:
                stats['SeqMismatch'] += 1

            stats['CommonAlignedBases'] += aln_diff['bases']
            stats['CommonMatchingBases'] += aln_diff['cons_score']
            stats['PerQueryBaseSim'].append(aln_diff['cons_score'] / float(aln_diff['bases']))
            stats['PerQueryBaseSimClipped'].append(float(
                aln_diff['cons_score']) / min(segments[0].infer_query_length(), segments[1].infer_query_length()))

            if is_coarse_match(aln_diff, coarse_tolerance):
                stats['CoarseMatches'] += 1

            stats[aln_one]['AlignedBases'] += aln_diff['bases']
            stats[aln_two]['AlignedBases'] += aln_diff['bases']

        # Read from first BAM is aligned:
        elif aln_diff['mapped'] == (True, False):
            stats[aln_one]['AlignedQueries'] += 1
            stats[aln_one]['AlignedBases'] += aln_diff['bases_one']
            stats[aln_two]['UnalignedQueries'] += 1
        # Read from second BAM is aligned:
        elif aln_diff['mapped'] == (False, True):
            stats[aln_two]['AlignedQueries'] += 1
            stats[aln_two]['AlignedBases'] += aln_diff['bases_two']
            stats[aln_one]['UnalignedQueries'] += 1
        # Both unaligned:
        elif aln_diff['mapped'] == (False, False):
            stats[aln_one]['UnalignedQueries'] += 1
            stats[aln_two]['UnalignedQueries'] += 1

    if stats['CommonAlignedBases'] > 0:
        stats['AlignedSimilarity'] = stats['CommonMatchingBases'] / \
            float(stats['CommonAlignedBases'])
    else:
        stats['AlignedSimilarity'] = 0.0
    return stats


def aligned_pairs_to_matches(aligned_pairs, offset):
    """Convert aligned pairs into a sequence of reference positions.

    :param aligned_pairs: Iterator of aligned pairs.
    :param offset: Offset at the beggining of the sequences.
    :returns: Iterator of reference positions aligned to the sequences positions.
    :rtype: generator
    """
    ref_pos_iter = (pair[1] for pair in aligned_pairs if pair[0] is not None)
    return chain([False] * offset, ref_pos_iter)


def calc_consistency_score(segment_one, segment_two, offset_one, offset_two):
    """Calculate the number of bases aligned to the same reference bases in two
    alignments.
    :param segment_one: Pysam aligned segments.
    :param segment_two: Pysam aligned segments.
    :param offset_one: Hard clipping offset for the first alignment.
    :param offset_two: Hard clipping offset for the second alignment.
    :retruns: Number of matching base alignments.
    :rtype: int
     """
    matches_one = aligned_pairs_to_matches(segment_one.get_aligned_pairs(), offset_one)
    matches_two = aligned_pairs_to_matches(segment_two.get_aligned_pairs(), offset_two)

    score = 0
    for matches in izip(matches_one, matches_two):
        if matches[0] == matches[1]:
            score += 1

    return score


def get_hard_clip_offset(aln):
    """Get hard clipping offset from alignment.

    :param aln: Pysam aligned segment.
    :returns: Hard clipping offset.
    :rtype: int
    """
    op = aln.cigartuples[0]
    if op[0] == 5:
        return op[1]
    return 0


def count_clipped(aln, target_op):
    """Count hard clipped bases in aligned segment.

    :param aln: Pysam aligned segement.
    :param target_op: CIGAR operation.
    :returns: Number of hard clipped bases in segment.
    :rtype: int
    """
    if aln.is_unmapped:
        return 0
    return sum(op[1] for op in (aln.cigartuples[0], aln.cigartuples[-1]) if op[0] == target_op)


def compare_alignments(segment_one, segment_two, strict_flags=False):
    """Count reads mapping to references in a BAM file.

    :param alignment_file: BAM file.
    :param min_aln_qual: Minimum mapping quality.
    :returns: Dictionary with read counts per reference.
    :rtype: dict
    """

    # FIXME: comparison logic does not support multiple alignments, split
    # alignments

    # Check if query names match:
    if segment_one.query_name != segment_two.query_name:
        raise Exception("Mismatched query names: {} {} Note that these should match exactly!".format(
            segment_one.query_name, segment_two.query_name))

    aln_diff = OrderedDict([
        ('mapped', None),
        ('ref_match', None),
        ('dir_match', None),
        ('flag_match', None),
        ('seq_match', None),
        ('bases', None),
        ('start_pos', None),
        ('end_pos', None),
        ('hard_clipped', None),
        ('soft_clipped', None),
        ('cons_score', None),
        ])

    aln_diff['mapped'] = (not segment_one.is_unmapped, not segment_two.is_unmapped)

    # Count hard clipped bases:
    aln_diff['hard_clipped'] = (count_clipped(segment_one, 5), count_clipped(segment_two, 5))

    # Count soft clipped bases:
    aln_diff['soft_clipped'] = (count_clipped(segment_one, 4), count_clipped(segment_two, 4))

    # One or both unmapped:
    if aln_diff['mapped'] == (False, False):
        return aln_diff
    elif aln_diff['mapped'] == (True, False):
        aln_diff['bases_one'] = segment_one.infer_query_length()
        return aln_diff
    elif aln_diff['mapped'] == (False, True):
        aln_diff['bases_two'] = segment_two.infer_query_length()
        return aln_diff

    # Mismatch in reference name:
    if segment_one.reference_name != segment_two.reference_name:
        aln_diff['ref_match'] = False
        return aln_diff
    else:
        aln_diff['ref_match'] = True

    # Mismatch in orientation:
    if segment_one.is_reverse != segment_two.is_reverse:
        aln_diff['dir_match'] = False
        return aln_diff
    else:
        aln_diff['dir_match'] = True

    # Perform strict flag checking:
    if strict_flags:
        if segment_one.flag != segment_two.flag:
            aln_diff['flag_match'] = False
            return aln_diff
        else:
            aln_diff['flag_match'] = True

    # Check if read sequences differ: Now should work with hard clipping.
    if segment_one.query_sequence != segment_two.query_sequence:
        aln_diff['seq_match'] = False
    else:
        aln_diff['seq_match'] = True

    # Register number of bases:
    aln_diff['bases'] = max(segment_one.infer_query_length(), segment_two.infer_query_length())

    # Register reference positions:
    aln_diff['start_pos'] = (segment_one.reference_start, segment_two.reference_start)
    aln_diff['end_pos'] = (segment_one.reference_end, segment_two.reference_end)

    # Figure out hard clipping offsets:
    offset_one = get_hard_clip_offset(segment_one)
    offset_two = get_hard_clip_offset(segment_two)

    # Register number of bases aligned the same way:
    aln_diff['cons_score'] = calc_consistency_score(
        segment_one, segment_two, offset_one, offset_two)

    return aln_diff
