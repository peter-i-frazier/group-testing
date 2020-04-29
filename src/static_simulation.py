from population import Population
from group_testing import SymptomaticIndividualTest
import numpy as np

def _safe_divide(a,b):
    if b == 0:
        return 0
    else:
        return a / float(b)

class StaticSimulation:
    # Runs a simulation over the initial population to estimate number of quarantined, number of people that should
    # have been quarantined but weren't, and number of tests

    def __init__(self, population, group_test):
        # Here we give some guidance on the population passed. In particular, we only depend on the parameters that
        # influence the infection status of the initial population.
        #
        # The size of the population (determined by n_households, household_size)
        # can be set to be the group_test.get_group_size(), or to a small multiple of that.
        # If that is a pain to do, you can instead set it to a large number.
        #
        # The initial population also depends on
        #   initial_prevalence,
        #   daily_secondary_attack_rate (the "daily" part of this is not important for this simulation)
        # It does not depend on:
        #   disease_length, time_until_symptomatic, non_quarantine_alpha,
        #   fatality_pct, daily_outside_infection_pct, outside_symptomatic_prob, initial_quarantine
        self.population = population
        self.group_test = group_test

    def __rep__(self):
        # Run one simulation of the group testing protocol on the population
        # Runs self.population.reset() when it starts
        # Returns   number of people quarantined,
        #           quarantine false negatives (number of people that should have been quarantined, but weren't),
        #           number_of_tests

        self.population.reset()
        attempts = 0
        while self.population.get_num_infected() == 0 and attempts < 400:
            self.population.reset()
            attempts += 1
        if self.population.get_num_infected() == 0:
            raise(Exception("could not infect population"))
        self.test_results, number_of_tests = self.group_test.test(self.population)
        #import pdb; pdb.set_trace()

        quarantines = 0
        quarantine_false_negatives = 0
        quarantine_false_positives = 0
        for (i,j), test_detected_presence in self.test_results.items():
            if test_detected_presence:
                quarantines += 1
                if (i,j) not in self.population.infected_individuals:
                    quarantine_false_positives += 1
            elif (i,j) in self.population.infected_individuals:
                # no virus detected but individual is positive
                quarantine_false_negatives += 1


        return { 'QFN' : quarantine_false_negatives,
                 'QFP' : quarantine_false_positives,
                 'tests' : number_of_tests,
                 'quarantines' : quarantines,
                 'positives' : self.population.get_num_infected(),
                 }


    def sim(self, nreps):
        # Runs the simulation for many replications and aggregates the output
        quarantine_false_negatives = [0]*nreps
        quarantine_false_positives = [0]*nreps
        number_of_tests = [0]*nreps
        quarantines = [0]*nreps
        positives = [0]*nreps
        for m in range(nreps):
            r = self.__rep__()
            quarantine_false_negatives[m]=r['QFN']
            quarantine_false_positives[m]=r['QFP']
            number_of_tests[m] = r['tests']
            quarantines[m] = r['quarantines']
            positives[m] = r['positives']

        n = self.population.get_num_individuals()
        # quarantine false negative rate is the number of quarantine false negatives,
        # i.e., positives that were incorrectly reported, divided by the overall number of positives
        #QFNR = np.mean(quarantine_false_negatives) / np.mean(positives)
        QFNR = np.mean([_safe_divide(QFN,P) for QFN, P in zip(quarantine_false_negatives, positives)])
        # similarly, quarantine false positive rate is the number of negatives that were incorrectly reported divided
        # by the overall number of negatives
        #QFPR = np.mean(quarantine_false_positives) / (n - np.mean(positives))
        QFPR = np.mean([_safe_divide(QFP,n - P) for QFP, P in zip(quarantine_false_positives, positives)])
        tests_per_person = np.mean(number_of_tests) / n
        quarantines_per_person = np.mean(quarantines) / n

        # TODO: should also return standard errors for these
        return QFNR, QFPR, tests_per_person, quarantines_per_person


    @staticmethod
    def __test__():
        pop = Population(n_households=100,
                         household_size=1,
                         initial_prevalence=.01,
                         disease_length=0,
                         time_until_symptomatic=0,
                         non_quarantine_alpha=0,
                         daily_secondary_attack_rate=.34,
                         fatality_pct=0,
                         daily_outside_infection_pct=0,
                         outside_symptomatic_prob=0,
                         initial_quarantine=0)

        group_test = SymptomaticIndividualTest(false_negative_rate=0.3)
        s = StaticSimulation(pop,group_test)
        print(s.__rep__())
        print(s.sim(10))


if __name__ == '__main__':
    StaticSimulation.__test__()









