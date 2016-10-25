#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from wub.bam import common as bam_common
from wub.bam import sam_writer
from wub.util import seq as seq_util

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Generate SAM records for the reads present in the input fastq but missing from
    the input SAM/BAM.
    """)
parser.add_argument(
    '-f', metavar='format', type=str, help="Input/output format (SAM).", default='SAM')
parser.add_argument(
    '-q', metavar='fastq', type=str, help="Input fastq.", required=True)
parser.add_argument(
    'infile', metavar='input_file', type=str, help="Input file.")
parser.add_argument(
    'outfile', metavar='output_file', type=str, help="Output SAM file.")

if __name__ == '__main__':
    args = parser.parse_args()

    input_iter = bam_common.pysam_open(args.infile, args.f).fetch(until_eof=True)

    # Get SAM record names:
    sam_names = [record.query_name for record in input_iter]

    writer = sam_writer.SamWriter(args.outfile)

    for read in seq_util.read_seq_records(args.q, 'fastq'):
        if read.id not in sam_names:
            qual = seq_util.quality_array_to_string(read.etter_annotations["phred_quality"])
            sam_record = writer.new_sam_record(qname=read.id, flag=4, rname="*", pos=0, mapq=0, cigar="*", rnext="*",
                                               pnext=0, tlen=0, seq=str(read.seq), qual=qual, tags="AS:i:0")
            writer.write(sam_record)

    writer.close()
