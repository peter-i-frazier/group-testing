from math import ceil
from random import shuffle
import numpy as np

class EmptyTest:
    def __init__(self, infection_pct_blf=0.01):
        pass
    def test(self, population):
        grp_test_data = {'num_grps': 0,
                        'grp_size': 0,
                        'num_negative_grps': 0}
        test_results = {}
        groups = {}
        return test_results, groups, grp_test_data

class GollierGroupTest:
    def __init__(self, infection_pct_blf=0.01, FNR=0):
        self.infection_pct_blf = infection_pct_blf

    def test(self, population):
        # TODO: need to update infection_pct_blf each round
        pop_size = population.get_population_size()
        grp_size = -1 / np.log(1-self.infection_pct_blf)

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
        num_negative_grps = 0
        for group_idx, grp_individuals in groups.items():
            result = population.any_infectious(grp_individuals)
            test_results[group_idx] = result
            if not result:
                num_negative_grps += 1

        grp_test_data = {'num_grps':num_grps, 'grp_size': grp_size,
                        'num_negative_grps': num_negative_grps}
        
        return test_results, groups, grp_test_data 


