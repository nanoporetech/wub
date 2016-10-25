# -*- coding: utf-8 -*-

import pysam


def pysam_open(alignment_file, in_format='BAM'):
    """Open SAM/BAM file using pysam.

    :param alignment_file: Input file.
    :param in_format: Format (SAM or BAM).
    :returns: pysam.AlignmentFile
    :rtype: pysam.AlignmentFile
    """
    if in_format == 'BAM':
        mode = "rb"
    elif in_format == 'SAM':
        mode = "r"
    else:
        raise Exception("Invalid format: {}".format(in_format))

    aln_iter = pysam.AlignmentFile(alignment_file, mode)
    return aln_iter
