from prism.hmm import HMMModel
from collections import Counter

import numpy as np
import multiprocessing as mp
import prism.util as util

import random
import cleanlog
import glob
import os
import sys

logger = cleanlog.ColoredLogger('insillico-proofreading')
PSEUDO_COUNT = 1

def hmm_proba(base_pattern, target_patterns, num_base_pattern, hmm_params):
    l = []
    for b in np.split(base_pattern, num_base_pattern):
        model = HMMModel(b, **hmm_params)
        l.append([model.proba(pattern) for pattern in target_patterns])
    return np.array(l)

def proofread_given_n_template(header, patterns, hmm_params, n_t=2, num_iter=100):
    # Abbreviations.
    #
    # l_p: Length of patterns.
    # n_p: Number of patterns.
    # n_t: Number of base patterns.
    # t: Template pattern.
    # l: Likelihood.
    # wl: Weighted likelihood.
    # post: Posterior probabilities.

    l_p, n_p = len(patterns[0]), len(patterns)
    # Initially we choose totally random base patterns.
    # We hope that the EM algorithm will find optimal template pattern!
    t = util.random_pattern(len(patterns[0]), n_t)
    # The uniform prior will be used.
    prior = (np.ones(n_t) / n_t).reshape([n_t, 1])
    
    prev_assignment = np.array([-1] * n_p)
    for i in range(num_iter):
        l = hmm_proba(t, patterns, n_t, hmm_params)

        # Posterior can be expressed as a proportion of weighted likelihoods (wl).
        wl = prior * l
        post = wl / wl.sum(axis=0)

        # Each pattern will be assigned to the template pattern with maximum posterior probability.
        assignment = np.argmax(post, axis=0)

        # If assignments converged, break the loop.
        if np.all(assignment == prev_assignment):
            break
        prev_assignment = assignment.copy()

        new_t = []
        # If there are some base patterns with no assigned patterns,
        # add new base patterns.
        n_assigned_template = len(set(assignment))
        if n_assigned_template != n_t:
            for _ in range(n_t - n_assigned_template):
                new_t.append(util.random_pattern(l_p, 1))

        # Compute new base patterns based with majority voting of the patterns assigned to each of them.
        for j in set(assignment):
            new = np.array((patterns[assignment == j].mean(axis=0) >= 0.5).astype(np.int32))
            new_t.append(new)
        t = np.array(new_t).flatten()

        prior = post.sum(axis=1).reshape([n_t, 1])
        # If it's not the last iteration, add pseudocount to the prior.
        if i != num_iter - 1:
            prior = prior + PSEUDO_COUNT * (np.ones(n_t) / n_t).reshape([n_t, 1])
        prior = prior / prior.sum()

    base_patterns = np.split(t, n_t)
    l = hmm_proba(t, patterns, n_t, hmm_params)
    wl = prior * l
    post = wl / wl.sum(axis=0)
    assignment = np.argmax(post, axis=0)

    sum_loglikelihood = 0
    for i, base_pattern in enumerate(base_patterns):
        model = HMMModel(base_pattern, **hmm_params)
        for j in range(len(patterns)):
            if assignment[j] == i:
                sum_loglikelihood += np.log(model.proba(patterns[j]))

    # Number of parameters.
    k = n_t * l_p  

    # Compute Bayesian Information Criterion for model selection.
    bic = np.log(n_p * l_p) * k - 2 * sum_loglikelihood

    error_corrected_patterns = np.array([base_patterns[a] for a in assignment])
    return bic, header, error_corrected_patterns

def proofread(params, num_iter=100):
    header, patterns, hmm_params = params
    logger.debug(f"Processing {header.split(';')[0]}")

    len_pattern = len(patterns[0])
    num_pattern = len(patterns)

    max_cluster = min(int(num_pattern / len_pattern), 10)

    if max_cluster == 0:
        _, header, error_corrected_patterns = proofread_given_n_template(header, patterns, hmm_params, 10, num_iter)
        return header, error_corrected_patterns

    _, header, error_corrected_patterns = min([proofread_given_n_template(header, patterns, hmm_params, n, num_iter) for n in range(1, max_cluster + 1)], key=lambda tup: tup[0])

    return header, error_corrected_patterns
