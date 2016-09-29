# -*- coding: utf-8 -*-

import numpy as np

from wub.util import seq as seq_util

uniform = [0.25, 0.25, 0.25, 0.25]

def random_base(probs=uniform):
    return np.random.choice(seq_util.bases, p=probs)

def random_base_except(excluded, probs=uniform):
    # In progress. Must re-normalise probabilities.
    filtered_bases = seq_util.bases[:].remove(excluded)
    return np.random.choice(filtered_bases)
