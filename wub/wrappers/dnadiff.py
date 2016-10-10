# -*- coding: utf-8 -*-
""" Wrapper for mummer's dnadiff """

import os
import re
from collections import defaultdict
from collections import namedtuple
import subprocess
import tempfile
from subprocess import STDOUT

dnadiff_extensions = (
    '.1coords', '.1delta', '.delta', '.mcoords', '.mdelta', '.qdiff', '.rdiff', '.report', '.snps')

Property = namedtuple('Property', 'ref query')
PropertyWithPerc = namedtuple('PropertyWithPerc', 'ref ref_perc query query_perc')


def dnadiff(reference, query, working_directory=None, cleanup=True):
    """Run dnadiff on reference and query fasta and parse results.

    :param reference: Reference fasta.
    :param query: Query fasta.
    :param working_directory: Write output in this directory if specified.
    :param cleanup: Delete dnadiff output after parsing if True.
    """
    reference = os.path.abspath(reference)
    query = os.path.abspath(query)
    work_dir = working_directory

    if not os.path.exists(reference):
        raise Exception("Reference fasta {} does not exists!".format(reference))
    if not os.path.exists(query):
        raise Exception("Target fasta {} does not exists!".format(query))
    if work_dir is not None and not os.path.exists(work_dir):
        raise Exception("Working directory {} does not exists!".format(work_dir))

    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix='dnadiff_')

    old_dir = os.getcwd()
    os.chdir(work_dir)

    command = ['dnadiff', reference, query]
    try:
        log = subprocess.check_output(command, stderr=STDOUT)
    finally:
        os.chdir(old_dir)

    report_file = os.path.join(work_dir, 'out.report')
    output = open(report_file, 'r').read()

    results = parse_dnadiff_report(report_file)

    if cleanup:
        cleanup_dnadiff_report(work_dir)
        if working_directory is None:
            os.rmdir(work_dir)

    return results, output, log


def cleanup_dnadiff_report(directory, prefix='out'):
    """Cleanup dnadiff output files in the specified directory.

    :param directory: Output directory.
    :param prefix: Output prefix.
    :returns: None
    :rtype: object
    """
    for ext in dnadiff_extensions:
        name = prefix + ext
        path = os.path.join(directory, name)
        if os.path.exists(path):
            os.unlink(path)


def _parse_dnadiff_into_sections(report_file):
    """Parse dnadiff output lines into sections."""
    report_fh = open(report_file, 'r')
    section = "NO_SECTION"
    sections = defaultdict(list)
    for line in report_fh:
        line = line.strip()
        if len(line) == 0:
            continue
        if line.startswith('/') or line.startswith('NUCMER') or line.startswith('[REF]'):
            continue
        if line.startswith('['):
            section = line
            section = section.replace('[', '')
            section = section.replace(']', '')
        else:
            sections[section].append(line)
    return sections


def _parse_percent_field(field):
    """Parse dnadiff field with percent value."""
    tmp = field.split('(')
    perc = tmp[1].replace(')', '')
    perc = perc.replace('%', '')
    return float(tmp[0]), float(perc)


def _parse_simple_section(lines):
    """Parse a simple dnadiff report section."""
    results = {}
    for line in lines:
        tmp = re.split("\s+", line)
        if '%' not in tmp[1] and '%' not in tmp[2]:
            results[tmp[0]] = Property(float(tmp[1]), float(tmp[2]))
        else:
            ref_prop, ref_prop_perc = _parse_percent_field(tmp[1])
            query_prop, query_prop_perc = _parse_percent_field(tmp[2])
            results[tmp[0]] = PropertyWithPerc(ref_prop, ref_prop_perc, query_prop, query_prop_perc)
    return results


def _parse_complex_section(lines):
    """Parse a complex dnadiff report section."""
    section = "NO_SECTION"
    sections = defaultdict(list)
    results = defaultdict(dict)
    # Parse alignment section into subsections:
    for line in lines:
        if len(line) == 0:
            continue
        # FIXME: Very specific to current dnadiff output:
        if line.startswith('1-to-1') or line.startswith('M-to-M') or line.startswith('Total'):
            tmp = re.split("\s+", line)
            section = tmp[0]
            results[section]['Number'] = Property(float(tmp[1]), float(tmp[2]))
        else:
            sections[section].append(line)

    # Parse subsections and update results dictionary:
    for section, lines in sections.iteritems():
        parsed = _parse_simple_section(lines)
        for name, prop in parsed.iteritems():
            results[section][name] = prop
    return results


def parse_dnadiff_report(report_file):
    """Parse dnadiff report file.

    :param report_file: dnadiff report output.
    :returns: Data structure with parsed results.
    :rtype: dict
    """
    sections = _parse_dnadiff_into_sections(report_file)

    results_sequences = _parse_simple_section(sections['Sequences'])
    results_bases = _parse_simple_section(sections['Bases'])
    results_features = _parse_simple_section(sections['Feature Estimates'])
    results_alignments = _parse_complex_section(sections['Alignments'])
    results_snps = _parse_complex_section(sections['SNPs'])

    results = {
        'Sequences': results_sequences,
        'Bases': results_bases,
        'Features': results_features,
        'Alignments': results_alignments,
        'SNPs': results_snps,
    }
    return results
