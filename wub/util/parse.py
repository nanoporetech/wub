# -*- coding: utf-8 -*-

"""Utilities to parse strings into various data sructures."""

def separated_list_to_floats(separated_list, separator=","):
    """ Convert a separated list into a list of floats.

    :param separated_list: A separated list as string.
    :param separator: List separator.
    :returns: List of floats.
    :rtype: list
    """
    return [float(element) for element in separated_list.split(separator)]
