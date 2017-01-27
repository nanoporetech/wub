# -*- coding: utf-8 -*-

from collections import OrderedDict


class SamWriter:

    """ Simple class to write SAM files. """

    def __init__(self, out_file, header=None):
        """ Initialise SAM writer object """
        self.out_file = out_file
        self.header = header
        self.out_handler = open(out_file, 'w')
        if header is not None:
            self._write_header()

    def _write_header(self):
        """Write SAM header."""
        for record_type, records in self.header.iteritems():
            for record in records:
                self.out_handler.write("@{}".format(record_type))
                for key, value in record.iteritems():
                    self.out_handler.write("\t{}:{}".format(key, value))
                self.out_handler.write("\n")

    def new_sam_record(self, qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags):
        """Create new SAM record structure.

        :param self: object
        :param qname: Read name.
        :param rname: Reference name.
        :param pos: Position in reference.
        :param mapq: Mapping quality.
        :param cigar: CIGAR string.
        :param rnext: Reference of next read.
        :param pnext: Position of next read.
        :param tlen: Template length.
        :param seq: Read sequence.
        :param qual: Base qualities.
        :param tags: Optional tags.
        :returns: SAM record.
        :rtype: OrderedDict
        """
        record = OrderedDict()

        record['QNAME'] = qname
        record['FLAG'] = flag
        record['RNAME'] = rname
        record['POS'] = pos
        record['MAPQ'] = mapq
        record['CIGAR'] = cigar
        record['RNEXT'] = rnext
        record['PNEXT'] = pnext
        record['TLEN'] = tlen
        record['SEQ'] = seq
        record['QUAL'] = qual
        record['TAGS'] = tags

        return record

    def write(self, record):
        """Write SAM record to file.

        :param self: object
        :param record: SAM record.
        :returns: None
        :rtype: object
        """
        self.out_handler.write("{}\n".format("\t".join(map(lambda x: str(x), record.itervalues()))))

    def close(self):
        """Close SAM file.

        :param self: object
        :returns: None
        :rtype: object
        """
        self.out_handler.flush()
        self.out_handler.close()
