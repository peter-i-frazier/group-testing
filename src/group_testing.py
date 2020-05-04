from math import ceil, sqrt, floor
import random
import numpy as np

def fnr_from_dilution(prevalence_in_group, group_size):
    i = prevalence_in_group
    n = group_size
    # 5% error rate for i = n = 1
    # 10% error rate for i = 1; n = 32
    fnr = 0.05 **(i/(n**0.076))
    return  fnr

class MatrixGroupTest:
    def __init__(self, group_size, false_negative_rate, false_positive_rate, fnr_at_swab_level=False):
        # group_size: size of group on which each test is run.  That is, the method
        #            forms 'matrix groups' of dimension group_size * group_size
        # false_negative_rate: how frequently does a test return false for a positive individual
        # fnr_at_swab_level: boolean. are false negatives at the swab level or the test level?
        #                    If they are at the swab level, then an error will cause all tests on an individual's swab
        #                    to come back negative.  If they are at the test level, then each test has an error
        #                    independently.
        # group_by_household: boolean. form groups by household or by individuals?
        self.group_size = int(group_size)
        self.false_negative_rate = false_negative_rate
        self.false_positive_rate = false_positive_rate
        self.fnr_at_swab_level = fnr_at_swab_level


    def _test_group(self, population, group):
        # test a particular group for the presence of the disease
        # population: the relevant population object
        # group: the list of households in the group
        for i in group:
            for j in range(population.household_sizes[i]):
                # For each swab that was tested in another group, we cache whether it was successful or not
                if (i,j) in population.infected_individuals and (i,j) not in self.was_swab_successful:
                    if self.fnr_at_swab_level:
                        # Randomize whether this new swab was successful.
                        self.was_swab_successful[(i,j)] = (random.random() > self.false_negative_rate)
                    else:
                        # Since we aren't doing testing at the swab level, all swabs are successful
                        self.was_swab_successful[(i,j)] = True

        # Number that are infected and have a successful swab
        num_in_group_infected = len([(i,j) for i in group for j in range(population.household_sizes[i])
                                    if (i,j) in population.infected_individuals 
                                    and self.was_swab_successful[(i,j)]])

        if self.fnr_at_swab_level and num_in_group_infected > 0:
            false_negative_pct = fnr_from_dilution(num_in_group_infected, len(group))
        else:
            # When false negatives are at the test level, we need to have an error for all infected individuals in the
            # test for it come back negative
            false_negative_pct = self.false_negative_rate ** (num_in_group_infected)

        if num_in_group_infected == 0:
            if random.random() < self.false_positive_rate:
                test_detected_presence = True
            else:
                # False negatives from swab failures will cause num_in_group_infected == 0 and will show up here
                test_detected_presence = False
        else:
            if random.random() < false_negative_pct:
                test_detected_presence = False
            else:
                test_detected_presence = True

        return test_detected_presence



    def test(self, population):
        # Test the entire population.  Runs _test_group several times and aggregates results.
        avg_household_size = population.get_avg_household_size()
        household_group_size = ceil(self.group_size / float(avg_household_size))
        household_matrix_size = household_group_size ** 2

        self.was_swab_successful = {}

        test_results = {}

        remaining_households = list(population.households)
        random.shuffle(remaining_households)


        n_tests = 0
        while len(remaining_households) > 0:
            matrix_size = int(min(household_matrix_size, len(remaining_households)))
            matrix_households = remaining_households[0:matrix_size]
            remaining_households = remaining_households[matrix_size:]

            group_size = int(ceil(sqrt(matrix_size)))

            row_groups = {i:[] for i in range(group_size)}
            column_groups = {j:[] for j in range(group_size)}
           
            # Form the row group assignments and the column group assignments.
            # All members of a household are assigned together.
            for idx, household in enumerate(matrix_households):
                row = floor(idx / float(group_size))
                col = idx % group_size
                row_groups[row].append(household)
                column_groups[col].append(household)

             
            # Perform the tests for both row and column groups.
            row_test_results = {}
            col_test_results = {}
            for idx in range(group_size):
                row_test_results[idx] = self._test_group(population, row_groups[idx])
                col_test_results[idx] = self._test_group(population, column_groups[idx])
                n_tests += 2


            # Look at the intersection of all positive tests.
            # All members of a household are positive if their row and column both come back positive.
            for idx, household in enumerate(matrix_households):
                row = floor(idx / float(group_size))
                col = idx % group_size
                test_detected_presence = (row_test_results[row] and col_test_results[col])
                for j in range(population.household_sizes[household]):
                    test_results[(household,j)] = test_detected_presence

        return test_results, n_tests



class SymptomaticIndividualTest:

    def __init__(self, false_negative_rate):
        self.false_negative_rate = false_negative_rate

    def test(self, population):
        symptomatic_individuals = population.get_symptomatic_individuals()
        quarantined_individuals = population.quarantined_individuals

        test_individuals = symptomatic_individuals.union(quarantined_individuals)

        test_results = {}
        num_tests = len(test_individuals)
        for (i,j) in test_individuals:
            if (i,j) not in population.infected_individuals:
                test_results[(i,j)] = False
            elif random.random() < self.false_negative_rate:
                test_results[(i,j)] = False
            else:
                test_results[(i,j)] = True

        return test_results, num_tests

class ThreeStageHierarchicalTest:

    def __init__(self,
            large_group_size,
            small_group_size,
            group_test_participation_rate,
            outer_false_negative_rate,
            inner_false_negative_rate):
        self.large_group_size = large_group_size
        self.small_group_size = small_group_size
        assert(self.small_group_size > 1)
        self.group_test_participation_rate = group_test_participation_rate
        self.inner_false_negative_rate = inner_false_negative_rate

        self.outer_group_test = HouseholdGroupTest(large_group_size,
                                                    group_test_participation_rate,
                                                    outer_false_negative_rate)
        self.inner_group_test = HouseholdGroupTest(small_group_size,
                                                    group_test_participation_rate,
                                                    inner_false_negative_rate)

        self.previous_positive_individuals = set()

    def get_group_size(self):
        return self.large_group_size

    def set_group_size(self, group_size):
        self.large_group_size = group_size


    def test(self, population):

        outer_test_results, num_outer_tests = self.outer_group_test.test(population, 
                                            non_participating_individuals=self.previous_positive_individuals)
        positive_outer_individuals = set([(i,j) for (i,j) in outer_test_results if outer_test_results[(i,j)]])
        
        inner_test_results, num_inner_tests = self.inner_group_test.test(population, 
                                                                test_individuals=positive_outer_individuals)

        positive_inner_individuals = self.previous_positive_individuals.union(set([(i,j) for (i,j) in inner_test_results if inner_test_results[(i,j)]]))
        
        positive_individuals = set()
        num_individual_tests = 0

        for (i,j) in positive_inner_individuals:
            num_individual_tests += 1

            if (i,j) in population.infected_individuals:

                if random.random() >  self.inner_false_negative_rate:
                    positive_individuals.add((i,j))

        total_tests = num_outer_tests + num_inner_tests + num_individual_tests
        test_results = {(i,j):False for (i,j) in outer_test_results}
        for (i,j) in positive_inner_individuals:
            if (i,j) in positive_individuals:
                test_results[(i,j)] = True
            else:
                test_results[(i,j)] = False

        return test_results, total_tests


class HouseholdGroupTest:
    def __init__(self, 
            group_test_size,
            group_test_participation_rate,
            false_negative_rate,
            false_positive_rate = 0):
        self.group_test_size = group_test_size
        self.group_test_participation_rate = group_test_participation_rate
        self.false_negative_rate = false_negative_rate
        self.false_positive_rate = false_positive_rate

    def get_group_size(self):
        return self.group_test_size

    def set_group_size(self, group_size):
        self.group_test_size = group_size

    def test(self, population, test_individuals=None, non_participating_individuals=None):

        # Make the test groups
        
        if test_individuals == None:
            test_individuals = set()

            if non_participating_individuals == None:
                non_participating_individuals = set()

            for (i,j) in population.population - population.fatality_individuals - non_participating_individuals:
                if (i,j) in population.quarantined_individuals or random.random() <= self.group_test_participation_rate:
                    test_individuals.add((i,j))

        number_of_groups = int(ceil(len(test_individuals) / self.group_test_size))
        ordered_test_individuals = list(test_individuals)
        random.shuffle(ordered_test_individuals)
        
        test_groups = {idx:set() for idx in range(number_of_groups)}
        group_idx_counter = 0
        already_grouped_individuals = set()

        for (i,j) in ordered_test_individuals:

            group_idx = group_idx_counter % number_of_groups
            group_idx_counter += 1
            
            if (i,j) not in already_grouped_individuals:

                test_groups[group_idx].add((i,j))
                already_grouped_individuals.add((i,j))

                for j_new in range(population.household_sizes[i]):
                    if (i,j_new) in test_individuals and j_new != j:
                        assert((i,j_new) not in already_grouped_individuals)
                        test_groups[group_idx].add((i,j_new))
                        already_grouped_individuals.add((i,j_new))
        
        # Test the groups
        #import pdb ; pdb.set_trace() 
        test_results = {}

        for group_idx in test_groups:
            num_in_group_infected = len([(i,j) for (i,j) in test_groups[group_idx] 
                                    if (i,j) in population.infected_individuals])
            if num_in_group_infected == 0:
                if random.random() < self.false_positive_rate:
                    test_detected_presence = True
                else:
                    test_detected_presence = False
            else:
                false_negative_pct = self.false_negative_rate ** num_in_group_infected
                if random.random() < false_negative_pct:
                    test_detected_presence = False
                else:
                    test_detected_presence = True

            for (i,j) in test_groups[group_idx]:
                test_results[(i,j)] = test_detected_presence

        return test_results, number_of_groups

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
        self.grp_size = -1 / np.log(1-self.infection_pct_blf)

    def set_group_size(self, group_size):
        self.grp_size = group_size

    def get_group_size(self):
        return self.grp_size

    def update_infection_blf(self, infection_blf):
        self.infection_pct_blf = infection_blf
        self.grp_size = -1 / np.log(1-self.infection_pct_blf)

    def test(self, population):
        # TODO: need to update infection_pct_blf each round
        pop_size = population.get_population_size()
        grp_size = self.grp_size
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


