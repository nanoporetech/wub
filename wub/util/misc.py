# -*- coding: utf-8 -*-
"""Yet uncategorised utility functions."""

import cPickle


def pickle_load(fname):
    """ Load object from pickle.

    :param fname: Input pickle file name.
    :returns: Object loaded from pickle file.
    :rtype: object

    """
    return cPickle.load(file(fname))


def pickle_dump(obj, fname):
    """Pickle object to file.

    :param obj: Object to be pickled.
    :fname: Output file name.
    :returns: The name of output file.
    :rtype: str

    """
    fh = open(fname, "w")
    cPickle.dump(obj, fh)
    fh.flush()
    fh.close()
    return fname
