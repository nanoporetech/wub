# -*- coding: utf-8 -*-
""" Parser functions for blastn outfmt 6. """


def _parse_coord_line(line):
    """ Parse a line from a blast outfmt 6 file. """
    fields = line.split()
    aln_record = {
        'query': fields[0],
        'ref': fields[1],
        'identity': float(fields[2]),
        'aln_length': int(fields[3]),
        'mismatch': int(fields[4]),
        'gapopen': int(fields[5]),
        'query_start': int(fields[6]),
        'query_end': int(fields[7]),
        'ref_start': int(fields[8]),
        'ref_end': int(fields[9]),
        'evalue': float(fields[10]),
        'bitscore': float(fields[11]),
        'strand': '+'
    }

    if aln_record['ref_start'] > aln_record['ref_end']:
        aln_record['strand'] = '-'
        aln_record['ref_start'], aln_record['ref_end'] = aln_record[
            'ref_end'], aln_record['ref_start']
    return aln_record


def parse_coords(input_object):
    """ Parse coordinates file produced by blastn outfmt 6.

    :param input_object: Input path or file hanlder.
    :returns: List of dictionaries with parsed records.
    :rtype: list
    """
    if type(input_object) == str:
        input_object = open(input_object, 'r')
    records = []
    for line in input_object:
        line = line.strip()
        records.append(_parse_coord_line(line))
    return records
