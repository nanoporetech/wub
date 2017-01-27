# -*- coding: utf-8 -*-

import numpy as np
from collections import OrderedDict

"""Utilities to parse strings into various data sructures."""


def separated_list_to_floats(separated_list, separator=","):
    """ Convert a separated list into a list of floats.

    :param separated_list: A separated list as string.
    :param separator: List separator.
    :returns: List of floats.
    :rtype: list
    """
    return [float(element) for element in separated_list.split(separator)]


def args_string_to_dict(args_string, elements_separator=",", keyvalue_separator=":"):
    """ Convert a two-level separated list into a dictionary.

    :param args_string: Two-level separated string.
    :param elements_separator: Separator between elements.
    :param keyvalue_separator: Separator between key/value pairs.
    :returns: dict
    :rtype: dict
    """
    if len(args_string) == 0:
        return {}
    pairs = [pair.strip() for pair in args_string.split(elements_separator)]
    elements = OrderedDict(pair.split(keyvalue_separator) for pair in pairs)
    parsed = OrderedDict((k.strip(), v.strip()) for k, v in elements.iteritems())
    return parsed


def interval_string_to_tuples(interval_string, elements_separator="|", interval_separator=","):
    """ Convert a two-level separated list into a dictionary.

    :param interval_string: Two-level separated string.
    :param elements_separator: Separator between elements.
    :param keyvalue_separator: Separator between interval boundaries.
    :returns: tuple
    :rtype: tuple
    """
    if len(interval_string) == 0:
        return tuple()
    pairs = [pair.strip() for pair in interval_string.split(elements_separator)]
    elements = OrderedDict(pair.split(interval_separator) for pair in pairs)
    parsed = tuple((int(k.strip()), int(v.strip())) for k, v in elements.iteritems())
    return parsed


def normalise_array(array):
    """ Normalise numpy array so the elments sum to 1.0.

    :param array: Input array.
    :returns: Normalised array.
    :rtype: numpy.array
    """
    temporary_array = array.astype(float)
    return temporary_array / np.sum(temporary_array)
