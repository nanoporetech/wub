#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pprint

from wub.util import misc

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Pretty print the contents of a pickle file.""")
parser.add_argument(
    'pickle', metavar='pickle_file', type=str, help="Input pickle file.")

if __name__ == '__main__':
    args = parser.parse_args()

    pprint.pprint(misc.pickle_load(args.pickle))
