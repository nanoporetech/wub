# -*- coding: utf-8 -*-
"""Utilities related to running external commands."""

from distutils.spawn import find_executable as exefind_distutils


def find_executable(command):
    """Find executable in path corresponding to a command.

    :param command: Command.
    :returns: Path to executable of False.
    :rtype: str
    """
    # In the future we might want to eliminate the dependency of
    # distutils.
    return exefind_distutils(command)
