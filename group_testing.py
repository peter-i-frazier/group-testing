from math import ceil
from random import shuffle
import numpy as np

class GollierGroupTest:
    def __init__(self, infection_pct_blf=0.01, FNR=0):
        self.infection_pct_blf = infection_pct_blf

    def test(self, population):
        # TODO: need to update infection_pct_blf each round
        pop_size = population.get_population_size()
        grp_size_individuals = -1 / np.log(1-self.infection_pct_blf)

        # for now ignore household correlation
        # grp_size_households = 

        num_grps = int(ceil(pop_size / grp_size))
        
        groups = {}
        for grp_num in range(num_grps):
            groups[grp_num] = []

        individuals = list(population.iter_individuals())
        shuffle(individuals)


        for individual, counter in zip(individuals, range(len(individuals))):
            grp = counter % num_grps
            groups[grp].append(individual)

        test_results = {}
        for group_idx, grp_individuals in groups.iteritems():
            result = population.any_infectious(grp_individuals)
            test_results[group_idx] = result
        
        return test_results, groups 


