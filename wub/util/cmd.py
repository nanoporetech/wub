# -*- coding: utf-8 -*-
"""Utilities related to running external commands."""

from distutils.spawn import find_executable as exefind_distutils
import sys


def find_executable(command):
    """Find executable in path corresponding to a command.

    :param command: Command.
    :returns: Path to executable of False.
    :rtype: str
    """
    # In the future we might want to eliminate the dependency of
    # distutils.
    return exefind_distutils(command)


def ensure_executable(command):
    """Find executable in path corresponding to a command and abort if not found.

    :param command: Command.
    :returns: None
    :rtype: object
    """
    # In the future we might want to eliminate the dependency of
    # distutils.
    if not find_executable(command):
        sys.stderr.write(
            "Required command \"{}\" not found in path! Aborting.!\n".format(command))
        sys.exit(127)
