#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Simulate genome with the specified number of chromosomes, length distribution and base composition.')
parser.add_argument(
    '-i', metavar='input', type=str, help="Input.")


if __name__ == '__main__':
    args = parser.parse_args()
