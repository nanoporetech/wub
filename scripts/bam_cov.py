#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import tqdm
import sys
from Bio import SeqIO
from wub.util import misc
from wub.bam import common as bam_common

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Produce refrence coverage table.""")
parser.add_argument(
    '-f', metavar='reference', type=str, help="Reference fasta.", required=True)
parser.add_argument(
    '-c', metavar='region', type=str, help="BAM region (None).", required=False, default=None)
parser.add_argument(
    '-t', metavar='tsv', type=str, default="bam_cov.tsv", help="Output TSV (bam_cov.tsv).", required=False)
parser.add_argument(
    '-q', metavar='aqual', type=int, default=0, help="Minimum alignment quality (0).")
parser.add_argument(
    '-Q', action="store_true", help="Be quiet and do not show progress bars.", default=False)
parser.add_argument('bam', metavar='bam', type=str, help="Input BAM file.")


def _process_bam(bam, out_tsv, chrom_lengths, region=None, min_aqual=0, verbose=True):
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

    tsv = open(out_tsv, "w")
    tsv.write(
        "Read\tRef\tStrand\tRefCov\tReadCov\tReadLength\tReadAlnLength\tRefLength\tRefAlnLength\tMapQual\n")

    for r in bam_iter:
        # Skip unmapped reads:
        if r.is_unmapped:
            continue
        # Skip if mapping quality is too low:
        if r.mapq < min_aqual:
            continue
        strand = '-' if r.is_reverse else '+'
        ref = r.reference_name
        ref_cov = r.reference_length / float(chrom_lengths[ref])
        read = r.query_name
        read_length = r.infer_read_length()
        mapq = r.mapping_quality
        read_aln_len = r.query_alignment_length
        read_cov = read_aln_len / float(read_length)
        ref_aln_length = r.reference_length

        tsv.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(read, ref, strand, ref_cov,
                                                                    read_cov, read_length, read_aln_len, chrom_lengths[ref], ref_aln_length, mapq))

    tsv.flush()
    tsv.close()


if __name__ == '__main__':
    args = parser.parse_args()
    verbose = not args.Q

    # Load reference lengths:
    references = SeqIO.index(args.f, format='fasta')
    chrom_lengths = {name: len(so) for name, so in references.items()}

    # Parse fragments:
    _process_bam(args.bam, args.t, chrom_lengths, args.c, args.q, verbose=verbose)
