# -*- coding: utf-8 -*-
""" Wrapper for mummer's dnadiff """

import os
import subprocess
import tempfile
from subprocess import STDOUT

dnadiff_extensions = (
    '.1coords', '.1delta', '.delta', '.mcoords', '.mdelta', '.qdiff', '.rdiff', '.report', '.snps')


def dnadiff(reference, query, working_directory=None, cleanup=True):

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
    print work_dir

    command = ['dnadiff', reference, query]
    try:
        output = subprocess.check_output(command, stderr=STDOUT)
    finally:
        os.chdir(old_dir)

    results = parse_dnadiff_report(os.path.join(work_dir, 'out.report'))

    if cleanup:
        cleanup_dnadiff_report(work_dir)
        if working_directory is None:
            os.rmdir(work_dir)

    return results, output


def cleanup_dnadiff_report(directory):
    for ext in dnadiff_extensions:
        name = 'out' + ext
        path = os.path.join(directory, name)
        if os.path.exists(path):
            os.unlink(path)


def parse_dnadiff_report(report_file):
    report_fh = open(report_file, 'r')
    for l in report_fh:
        print l
    pass
