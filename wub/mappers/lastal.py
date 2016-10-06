import os
import re
import pandas as pd
from collections import namedtuple, OrderedDict, defaultdict
import subprocess
from subprocess import Popen, PIPE, STDOUT
from Bio import SeqIO

from wub.util import seq as seq_util

lastdb_suffixes = ['bck', 'des', 'prj', 'sds', 'ssp', 'suf', 'tis']

record_template = 'score r_name r_start r_aln_len r_strand r_len r_aln '
record_template += 'q_name q_start q_aln_len q_strand q_len q_aln'
LastRecord = namedtuple('LastRecord', record_template)


def lastdb(ref_dir, ref_name, ref, executable='lastdb', **kwargs):
    """
    Runs lastdb on ref within ref_dir using the label ref_name if
    any errors thrown during runtime, files are checked for existence
    if all files accounted for, successful=False but no errors thrown.
    Otherwise, IOError or CalledProcessError thrown.


    :param ref_dir: directory you will find lastdb files in
    :param ref_name: name of the lastdb files e.g. a for a.prj..
    :param ref: filepath for reference file
    :param executable: path/executable for lastdb e.g. ont_lastdb
    :param kwargs: any -[arg] wanted see lastdb -h for details
    :return: True/False is successful with no errors and command run
    :raises: `IOError` if files don't exist
    :raises: `subprocess.CalledProcessError` for errors during runtime
    """
    # taken from metrichor-bio/alignment:

    ref_dir = os.path.abspath(ref_dir)
    ref = os.path.abspath(ref)
    chdir = 'cd {}'.format(ref_dir)
    kwargs = ['-{} {}'.format(a, b) for a, b in kwargs.items()]
    kwargs_str = ' '.join(kwargs)
    cmd = '{} {} {} {}'.format(executable, ref_name, ref, kwargs_str)
    command = '{}; {}'.format(chdir, cmd)

    # Because a cd is done in subprocess, if that fails cwd is used.
    if not os.path.exists(ref_dir):
        raise IOError('Directory not found: {}'.format(ref_dir))

    try:
        _ = subprocess.check_output(command, shell=True,
                                    stderr=subprocess.STDOUT)
        successful = True
    except subprocess.CalledProcessError as e:
        successful = False
        if not os.path.exists(ref):
            raise IOError('Reference not found: {}'.format(ref))
        elif not utilities.find_executable(executable):
            raise IOError('Executable not found: {}'.format(executable))
        elif self.check_lastdb_file(ref_dir, ref_name):
            raise e

    return os.path.join(ref_dir, ref_name)


def check_lastdb_files(ref_dir, name):
    """
    Check that all lastdb files with `name` label exist within directory

    :param ref_dir: directory to check for lastdb files
    :param name: label to search for e.g. 'a' for a.prj
    :return: list of missing extensions, [] if none missing
    """
    # taken from metrichor-bio/alignment:
    files = 0
    missing = []
    for suffix in lastdb_suffixes:
        lastdb_file = '.'.join([name, suffix])
        lastdb_path = os.path.join(ref_dir, lastdb_file)
        files += os.path.exists(lastdb_path)
        if not os.path.exists(lastdb_path):
            missing.append(suffix)
    return missing


def lastal_align(database, query, executable='lastal', **kwargs):
    """
    Runs lastal via subprocess

    :param database: database prefix
    :param query: filepath for the query file
    :param kwargs: -[args] wanted for lastal e.g. v='' for verbosity
    :return: alignment output
    """
    # some code taken from metrichor-bio/alignment:
    ref_dir, ref_name = os.path.split(database)
    chdir_cmd = 'cd {}'.format(ref_dir)
    kwargs = ['-{} {}'.format(a, b) for a, b in kwargs.items()]
    kwargs_str = ' '.join(kwargs)

    cmd = ' '.join([executable, kwargs_str, ref_name, query])
    command = '{}; {}'.format(chdir_cmd, cmd)

    p = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)

    # FIXME: reading whole alignment in memory, might not always be a good idea:
    return p.stdout.read()


def parse_lastal(res):
    """
    Parse raw lastal output records.
    :param res: Raw lastal results.
    :returns: Generator of lastal alignment records.
    :rtype: generator
    """
    lines = [l for l in res.split('\n') if not l.startswith('#')]
    re_score = re.compile('\s|=')
    re_aln = re.compile('\s+')
    for i, line in enumerate(lines):
        if line.startswith('a '):
            # Parse out score:
            score = int(re_score.split(line)[2])
            # Parse reference record:
            tmp = re_aln.split(lines[i + 1])[1:7]
            ref_name = tmp[0]
            ref_start = int(tmp[1])
            ref_aln_len = int(tmp[2])
            ref_strand = tmp[3]
            ref_len = int(tmp[4])
            ref_aln = tmp[5]
            # Parse query record:
            tmp = re_aln.split(lines[i + 2])[1:7]
            q_name = tmp[0]
            q_start = int(tmp[1])
            q_aln_len = int(tmp[2])
            q_strand = tmp[3]
            q_len = int(tmp[4])
            q_aln = tmp[5]
            # Yield record:
            record = LastRecord(score, ref_name, ref_start, ref_aln_len, ref_strand, ref_len,
                                ref_aln, q_name, q_start, q_aln_len, q_strand, q_len, q_aln)
            yield record


def filter_top_per_query(records):
    """Filter lastal alignment records keeping the best scoring one per query.

    :param records: A collection of LastRecord named tuples.
    :returns: A list of LastRecord named tuples.
    :rtype: list
    """
    top = OrderedDict()
    for alignment in records:
        if alignment.q_name not in top:
            top[alignment.q_name] = alignment
        elif alignment.score > top[alignment.q_name].score:
            top[alignment.q_name] = alignment
    return top.values()


def compare_genomes_lastal(ref_fasta, target_fasta, lastal_options, cleanup=True):
    """Compare a refrence set of sequences to a target set os sequences using lastal alignment.

    :param ref_fasta: Reference sequence set in fasta format.
    :param target_fasta: Target sequence set in fasta format.
    :param lastal_options: Options passed to lastal in a dictionary.
    :param cleanup: If True then lastal database files will be deleted.
    :returns: A pandas data frame with various per-alignment statistics.
    :rtype: DataFrame
    """
    ref_fasta = os.path.abspath(ref_fasta)
    ref_dir, ref_name = os.path.split(ref_fasta)
    db = lastdb(ref_dir, ref_name, ref_fasta)
    alignments = parse_lastal(lastal_align(db, target_fasta, **lastal_options))
    top_alignments = filter_top_per_query(alignments)

    stats = defaultdict(list)
    for aln in top_alignments:
        coverage = aln.r_aln_len / float(aln.r_len)
        aln_stats = seq_util.alignment_stats(aln.r_aln, aln.q_aln)
        stats['reference'].append(aln.r_name)
        stats['target'].append(aln.q_name)
        stats['accuracy'].append(aln_stats.accuracy)
        stats['coverage'].append(coverage)
        stats['ref_len'].append(aln.r_len)
        stats['ref_aln_len'].append(aln.r_aln_len)
        stats['aln_length'].append(len(aln.r_aln))
        stats['substitutions'].append(aln_stats.substitutions)
        stats['deletions'].append(aln_stats.deletions)
        stats['insertions'].append(aln_stats.insertions)
    column_order = ['reference', 'target', 'accuracy',
                    'coverage', 'ref_len', 'ref_aln_len', 'aln_length', 'substitutions', 'deletions', 'insertions']
    stats_frame = pd.DataFrame(stats)
    stats_frame = stats_frame[column_order]
    return stats_frame
