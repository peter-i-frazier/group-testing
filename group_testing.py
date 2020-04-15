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
    def __init__(self, infection_pct_blf=0.01):
        self.infection_pct_blf = infection_pct_blf

    def update_infection_blf(self, infection_blf):
        self.infection_pct_blf = infection_blf

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
            result = population.test_group(grp_individuals)
            test_results[group_idx] = result
            if not result:
                num_negative_grps += 1

        grp_test_data = {'num_grps':num_grps, 'grp_size': grp_size,
                        'num_negative_grps': num_negative_grps}
        
        return test_results, groups, grp_test_data 

class GollierHouseholdGroupTest:
    def __init__(self, infection_pct_blf=0.01, FNR=0):
        self.infection_pct_blf=0.01

    def update_infection_blf(self, infection_blf):
        self.infection_pct_blf = infection_blf

    def test(self, population):
        n_households = population.get_num_households()
        avg_house_size = population.get_avg_household_size()
        # unjustified heuristic to estimate household size
        household_infection_blf = 1 - (1 - self.infection_pct_blf)** avg_house_size
        grp_size = -1/np.log(1-household_infection_blf)

        num_grps = int(ceil(n_households / grp_size))
        
        groups = {}
        for grp_num in range(num_grps):
            groups[grp_num] = []

        households = list(population.iter_households())
        shuffle(households)


        for household, counter in zip(households, range(len(households))):
            grp = counter % num_grps
            for individual in population.iter_household_individuals(household):
                groups[grp].append(individual)

        test_results = {}
        num_negative_grps = 0
        for group_idx, grp_individuals in groups.items():
            result = population.test_group(grp_individuals)
            test_results[group_idx] = result
            if not result:
                num_negative_grps += 1

        grp_test_data = {'num_grps':num_grps, 'grp_size': grp_size,
                        'num_negative_grps': num_negative_grps}
        
        return test_results, groups, grp_test_data 

