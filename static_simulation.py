from population import Population
from group_testing import SymptomaticIndividualTest

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

    def rep(self):
        # Run one simulation of the group testing protocol on the population
        # Runs self.population.reset() when it starts
        # Returns   number of people quarantined,
        #           quarantine false negatives (number of people that should have been quarantined, but weren't),
        #           number_of_tests
        self.population.reset()
        self.test_results, number_of_tests = self.group_test.test(self.population)

        quarantined = 0
        quarantine_false_negatives = 0
        for (i,j), test_detected_presence in self.test_results.items():
            if test_detected_presence:
                quarantined += 1
            elif (i,j) in self.population.infected_individuals:
                # no virus detected but individual is positive
                quarantine_false_negatives += 1

        return quarantined, quarantine_false_negatives, number_of_tests

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
        print(s.rep())


if __name__ == '__main__':
    StaticSimulation.__test__()









