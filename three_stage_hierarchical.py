import numpy as np
from math import ceil

def sim_3stage(pool1_size = 100, pool2_size = 10, prevalence = 0.001):
    # INPUTS
    # pool1_size: number of people tested in the first pool
    # pool2_size: number of people tested in each of the sub-pools, if the first pool is positive
    # prevalence: fraction of the population that is infected

    # RETURNS: number of tests per person, in one simulation

    # x is a vector that represents infection status
    # An entry is positive if that individual is infected
    # We will simulate testing for everyone who is in the first pool
    x = np.random.rand(pool1_size) < prevalence

    n_tests = 1

    # Simulate whether the first pool is positive
    positive = np.any(x[1:pool1_size])
    if not(positive):
        return n_tests / pool1_size # No more tests need to be run

    n_subpools = ceil(pool1_size / pool2_size)
    n_tests += n_subpools # need to test all of the sub-pools
    for g in range(n_subpools):
        first = g * pool2_size # index of first individual in the pool
        last = min(first+pool2_size, pool1_size) - 1 # index of last individual in the pool
        positive = np.any(x[first:last])
        if positive:
            # Need to test everyone in this pool individually
            n_tests += last-first + 1

    return n_tests / pool1_size
