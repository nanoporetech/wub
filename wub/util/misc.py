# -*- coding: utf-8 -*-
"""Yet uncategorised utility functions."""

import pickle
import os.path


def get_fname(fname):
    '''
    get the file name without extension

    :param fname: file name
    :return: file name
    :rtype: str
    '''
    return os.path.splitext(os.path.basename(fname))[0]


def get_extension(fname):
    '''
    get the file extension

    :param fname: file name
    :return: file extention
    :rtype: str format '.*'
    '''
    return os.path.splitext(os.path.basename(fname))[1]


def _getextension(fast):
    '''
    finds and check for the correct extension

    :param fast: fastq or fasta file
    :return: "fastq" or "fasta"
    :rtype: str
    '''

    extension = get_extension(fast)
    if extension in ('.fa', '.fasta'):
        extension = "fasta"
    elif extension in ('.fq', '.fastq'):
        extension = "fastq"
    else:
        raise Exception('Incorrect file format')
        exit()
        # print >> sys.stderr, "Incorrect file format"
    return extension


def mkdir(path):
    '''
    if the dir does not exists it create it

    :param path: dir path
    :return: path
    :rtype: str
    '''
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def pickle_load(fname):
    """ Load object from pickle.

    :param fname: Input pickle file name.
    :returns: Object loaded from pickle file.
    :rtype: object

    """
    return pickle.load(file(fname))


def pickle_dump(obj, fname):
    """Pickle object to file.

    :param obj: Object to be pickled.
    :fname: Output file name.
    :returns: The name of output file.
    :rtype: str

    """
    pickle.dump(obj, file(fname, 'w'))
    return fname
