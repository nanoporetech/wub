# -*- coding: utf-8 -*-
"""Filter SAM/BAM records by various criteria."""

import itertools


def get_alignment_score(segement):
    """Get alignment score from pysam segment.

    :param segment: Pysam aligned segment.
    :returns: Alignment score.
    :rtype: int
    """

    score = 0
    try:
        score = segement.get_tag('AS')
    except:
        pass
    return score


def filter_top_per_query(records_iter):
    """Filter pysam records keeping top scoring per query. Assumes
    records are sorted by name.

    :param records_iter: Iterator of pysam aligned segments.
    :returns: Generator of filtered records.
    :rtype: generator
    """
    buff = []
    for rec in itertools.chain(records_iter, [None]):
        if len(buff) == 0:
            buff.append(rec)
        elif rec is None or buff[-1].query_name != rec.query_name:
            sorted_buff = sorted(buff, key=get_alignment_score, reverse=True)
            buff = [rec]
            yield sorted_buff[0]
        else:
            buff.append(rec)


def filter_query_coverage(records_iter, minimum_coverage):
    """Filter pysam records keeping the ones with sufficient query coverage.

    :param records_iter: Iterator of pysam aligned segments.
    :param minimum_coverage: Minimum fraction of covered query.
    :returns: Generator of filtered records.
    :rtype: generator
    """
    for rec in records_iter:
        if rec.is_unmapped:
            yield rec
        elif (float(rec.query_alignment_length) / rec.infer_query_length()) >= minimum_coverage:
            yield rec


def filter_ref_coverage(records_iter, minimum_coverage, header):
    """Filter pysam records keeping the ones with sufficient reference coverage.

    :param records_iter: Iterator of pysam aligned segments.
    :param minimum_coverage: Minimum fraction of covered reference.
    :param header: SAM header with reference lengths.
    :returns: Generator of filtered records.
    :rtype: generator
    """
    ref_lengths = dict((h['SN'], int(h['LN'])) for h in header['SQ'])
    for rec in records_iter:
        if rec.is_unmapped:
            yield rec
        elif (float(rec.query_alignment_length) / ref_lengths[rec.reference_name]) >= minimum_coverage:
            yield rec
