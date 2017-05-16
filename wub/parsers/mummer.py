# -*- coding: utf-8 -*-
""" Parser functions for mummer. """

import six


def _parse_coord_line(line):
    """ Parse a line from a mummer coordinate file. """
    fields = line.replace("|", "").split()
    aln_record = {
        'ref_start': int(fields[0]),
        'ref_end': int(fields[1]),
        'query_start': int(fields[2]),
        'query_end': int(fields[3]),
        'ref_len': int(fields[4]),
        'query_len': int(fields[5]),
        'identity': float(fields[6]),
        'ref': fields[7],
        'query': fields[8],
    }
    return aln_record


def parse_coords(input_object):
    """ Parse coordinates file produced by mummer.

    :param input_object: Input path or file hanlder.
    :returns: List of dictionaries with parsed records.
    :rtype: list
    """
    if type(input_object) == str:
        input_object = open(input_object, 'r')
    records = []
    for line in input_object:
        line = line.strip()
        if line.count('/') > 0:
            continue
        if line.count('NUCMER') > 0:
            continue
        if line.count('[') > 0:
            continue
        if line.count('=') > 0:
            continue
        if len(line) == 0:
            continue
        records.append(_parse_coord_line(line))
    return records
