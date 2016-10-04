# -*- coding: utf-8 -*-

import numpy as np

"""Utilities to parse strings into various data sructures."""


def separated_list_to_floats(separated_list, separator=","):
    """ Convert a separated list into a list of floats.

    :param separated_list: A separated list as string.
    :param separator: List separator.
    :returns: List of floats.
    :rtype: list
    """
    return [float(element) for element in separated_list.split(separator)]


def normalise_array(array):
    """ Normalise numpy array so the elments sum to 1.0.

    :param array: Input array.
    :returns: Normalised array.
    :rtype: numpy.array
    """
    temporary_array = array.astype(float)
    return temporary_array / np.sum(temporary_array)
