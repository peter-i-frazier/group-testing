import numpy as np

def sim_2stage(pool_size = 32, prevalence = 0.001):
    # INPUTS
    # pool_size: number of people tested in the pool
    # prevalence: fraction of the population that is infected
    # RETURNS: number of tests per person, in one simulation

    # x is a vector that represents infection status
    # An entry is positive if that individual is infected
    x = np.random.rand(pool_size) < prevalence

    n_tests = 1
    positive = x[1:pool_size].max()
    if positive:
        n_tests += pool_size
    return n_tests / pool_size
