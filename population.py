import numpy as np

class Population:

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

    def __init__(self, n_households, household_size_dist, prevalence, SAR, R0, d0, fatality_pct):
        # Initialize a population with non-trivial households
        # n_households:         the number of households in the population
        # household_size_dist:  a numpy array that should sum to 1, where household_size_dist[i] gives the fraction of
        #                       households with size i+1
        # prevalence:           prevalence of the disease in the population, used to simulate initial infection status
        # SAR:                  secondary attack rate, used for simulating initial infection status and for simulating
        #                       forward in time; SAR is defined as the probability that an infection occurs among
        #                       susceptible people within a household related to a confirmed case
        #
        assert np.isclose(np.sum(household_size_dist), 1.)

        self.infectious = []
        self.quarantined = []
        self.susceptible = []
        self.prevalence = prevalence
        self.n_households = n_households
        self.R0 = R0
        self.d0 = d0
        self.SAR = SAR
        self.fatality_pct = fatality_pct

        self.household_risk_status = [True] * self.n_households

        self.remaining_days_infected = {}
        self.deaths = set([])
        self.recoveries = set([])

        self.total_pop = 0


        for i in range(n_households):
            # generate household size h from household_size_dist
            h = int(np.random.choice(np.arange(1, len(household_size_dist)+1), 1, p=household_size_dist))
            self.total_pop += h
            # compute primary case probability = p*h/(1+SAR*(h-1))
            prob_prim = prevalence*h/(1+SAR*(h-1))
            infectious = [np.random.uniform() < prob_prim]
            quarantined = [True]        
            if h > 1:
                # if there are >1 members in the household, and there is a primary case,
                # generate secondary cases from Bin(h-1, SAR); otherwise, set everyone to be uninfected
                infectious.extend(np.random.binomial(1, SAR, h-1) * infectious[0] == True)
                quarantined.extend([True]*(h-1))
            
            susceptible = [not infected for infected in infectious]
            
            self.susceptible.append(np.array(susceptible))
            self.infectious.append(np.array(infectious))
            self.quarantined.append(np.array(quarantined))

        self.update_infection_days()

    def update_infection_days(self):
        for i in range(self.n_households):
            for j in range(len(self.infectious[i])):
                if self.infectious[i][j]:
                    if (i,j) not in self.remaining_days_infected:
                        self.remaining_days_infected[(i,j)] = self.d0
                    elif self.remaining_days_infected[(i,j)] > 0:
                        self.remaining_days_infected[(i,j)] -= 1
                        if self.remaining_days_infected[(i,j)] == 0:
                            if np.random.uniform() < self.fatality_pct:
                                self.deaths.add((i,j))
                            else:
                                self.recoveries.add((i,j))
                            self.susceptible[i][j] = False
                            self.quarantined[i][j] = False
                            self.infectious[i][j] = False



    def __check_indices(self,x):
        # Make sure that a passed set of individual array indices are valid for our population
        # TODO: This should be adapted for households?
        assert(max(x) < self.n)
        assert(min(x) >= 0)

    def infected(self, x):
        # Given a set of indices x, return the infection status of those individuals
        # Used by tests based on this class
        self.__check_indices(x)
        return self.infectious[x]

    def get_prevalence(self):
        # Get current prevalence = (number of infectious unquarantined) / (number of unquarantined)
        infected_counts = 0
        unquarantined_counts = 0

        for i in range(self.n_households):
            infected_counts += np.sum(self.infectious[i] & ~self.quarantined[i])
            unquarantined_counts += np.sum(~self.quarantined[i])
        if unquarantined_counts == 0:
            return 0
        else:
            return infected_counts / unquarantined_counts

    def step(self):
        # Simulate one step forward in time
        # Simulate how infectious individuals infect each other
        # Unquarantined susceptible people become infected w/ probability = alpha*current prevalence
        prevalence = self.get_prevalence()
        for i in range(self.n_households):
            household_infected = any(self.infectious[i])
            for j in range(len(self.quarantined[i])):
                if ~self.quarantined[i][j] and self.susceptible[i][j] and (not self.infectious[i][j]):
                    if household_infected:
                        secondary_prob = self.SAR # maybe should depend on total # of infected in household, but this is a start
                        new_secondary = np.random.uniform() < secondary_prob
                    else:
                        new_secondary = False
                    primary_prob = np.log(self.R0) * prevalence / self.d0
                    new_primary = np.random.uniform() < primary_prob
                    self.infectous[i][j] = new_primary or new_secondary

        self.update_infection_days()

    def get_population_size(self):
        return sum([len(household) for household in self.infectious])
    
    def get_num_households(self):
        return self.n_households

    def get_avg_household_size(self):
        return np.mean([len(household) for household in self.infectious])
    
    def iter_individuals(self):
        for i in range(self.n_households):
            for j in range(len(self.infectious[i])):
                yield (i,j)

    def any_infectious(self, grp_individuals):
        return any([self.infectious[i][j] for (i,j) in grp_individuals])

    def update_qurantine_status(self, test_results, groups):
        # unquarantine all confirmed negative results who do not also have
        # an infected family member
        self.household_risk_status = [False] * self.n_households
        for group_idx, result in test_results.iteritems():
            # for any positive group, set all household risk statuses to True
            # within that group
            if result:
                for i,_ in groups[group_idx]:
                    self.household_risk_status[i] = True

        # unquarantine all households with negative risk status
        for i, risk_status in enumerate(self.household_risk_status):
            if not risk_status:
                for j in range(len(self.quarantined[i])):
                    self.quarantined[i][j] = False


    def quarantine(self, x):
        # Put into quarantine all individuals whose indices are in the list x
        self.__check_indices(x)
        self.quarantined[x] = True

    def unquarantine(self, x):
        # Remove from quarantine all individuals whose indices are in the list x
        self.__check_indices(x)
        self.quarantined[x] = False

    def num_infected(self):
        infected_counts = 0
        for i in range(self.n_households):
            infected_counts += np.sum(self.infectious[i] and ~self.quarantined[i])
        return infected_counts

    #def num_recovered_dead(self):
    # for now we haven't taken into consideration the duration of infection and recovery yet
