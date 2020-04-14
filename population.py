import numpy as np

class population:

    """ This class describes a population of people and their infection status, optionally including information about their
    household organization and the ability to simulate infection status in an correlated way across households.
    It has methods for simulating infection status forward in time """

    def __init__(self, n, prevalence):
        # Create a population with n individuals in which the household sizes are all 1 with the given initial disease
        # prevalence

        self.infectious = np.random.rand(n) < prevalence[0]  # True if infectious
        self.quarantined = [False]*n    # True if quarantined

        # True if susceptible to the disease, i.e., has not died or developed immunity, and not currently infected
        self.susceptible = ~self.infectious

    def __init__(self, n_households, household_size_dist, prevalence, SAR, R0, d0):
        # Initialize a population with non-trivial households
        # n_households:         the number of households in the population
        # household_size_dist:  a numpy array that should sum to 1, where household_size_dist[i] gives the fraction of
        #                       households with size i+1
        # prevalence:           prevalence of the disease in the population, used to simulate initial infection status
        # SAR:                  secondary attack rate, used for simulating initial infection status and for simulating
        #                       forward in time
        #
        assert np.isclose(np.sum(household_size_dist),1.)

        # This code is for Yujia to write

        # simulate the household sizes

        # simulate their infection statuses

    def __check_indices(self,x):
        # Make sure that a passed set of individual array indices are valid for our population
        assert(max(x) < self.n)
        assert(min(x) >= 0)

    def infected(self, x):
        # Given a set of indices x, return the infection status of those individuals
        # Used by tests based on this class
        __check_indices(x)
        return self.infected[x]

    def step(self):
        # Simulate one step forward in time
        # Simulate how infectious individuals infect each other

        # For Yujia to fill in

    def quarantine(self, x):
        # Put into quarantine all individuals whose indices are in the list x
        __check_indices(x)
        self.quarantined[x] = True

    def unquarantine(self, x):
        # Remove from quarantine all individuals whose indices are in the list x
        __check_indices(x)
        self.quarantined[x] = False

    def num_infected(self):
        

    def num_recovered_dead(self):
