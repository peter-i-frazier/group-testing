import numpy as np
from math import ceil

def sim_dyadic(n = 10000, pools_per_person = 4, people_per_pool = 256, prevalence = 0.001):
    # INPUTS
    # n: Number of individuals from which pools are formed
    # pools_per_person: Number of pools in which each individual participates
    # people_per_pool: Number of people in each pool
    # prevalence: 0.001

    # Returns the number of tests done per person in this simulation

    n_pools = ceil(n * pools_per_person / people_per_pool) # Number of pools we will need
    pool_assignments = np.full([n, n_pools], False) # Boolean matrix telling us who has been assigned to what pool
    pool_counts = np.zeros(n_pools,dtype=np.int32) # counts of how many people have been assigned to each pool

    # Assign individuals to pools
    for i in range(n):

        # Find eligible pools for the assignment of our next individual, i.e., the pools that are not full yet
        eligible = np.argwhere(pool_counts<people_per_pool)

        # eligible is returned as an ndarray with shape (n_pools,1).
        # To convert it to a 1-d array, we use eligible[:,0]
        eligible = eligible[:,0]

        # Ideally we would place each person into pools_per_person.
        # Indeed, we have enough total testing slots (n_pools * people_per_pool) to do this.
        # However, we can fill in the pools unevenly --- when we have a small number of testing slots left, this can
        # be a problem and cause us to not have enough eligible pools.  For example, pools_per_person is 4, and we
        # only have 4 slots left, we might have 2 slots left in one pool and 1 each in 2 other pools.
        n_to_join = pools_per_person
        if len(eligible) < pools_per_person:
            n_to_join = len(eligible)
            print('Warning: ran out of eligible pools, using less than the desired number')

        # We should always have at least one eligible pool
        assert(n_to_join > 0)

        # Select pools for this individual to join, without replacement
        to_join = np.random.choice(eligible, n_to_join, replace=False)

        for g in to_join:
            pool_assignments[i][g] = True
            pool_counts[g] += 1

    # for each individual, randomize whether they are positive or not
    infection_status = np.random.rand(n) < prevalence

    # Run the tests
    n_tests = n_pools # This will track total number of tests done
    pool_results = np.full(n_pools, False)
    for g in range(n_pools):
        # A pool tests positive if anyone in the pool had a positive infection status
        pool_results[g] = np.any([infection_status[i] and pool_assignments[i][g] for i in range(n)])

    # Figure out how many people were in only positive pools
    for i in range(n):
        # We need to follow up if all pools in which the individual participated were positive
        # I.e., if, for all pools, the pool was either positive or the individual did not participate
        if np.all([pool_results[g] or not(pool_assignments[i][g]) for g in range(n_pools)]):
            n_tests += 1

    return(n_tests / n)
