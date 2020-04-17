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

    def __init__(self, n_households, household_size_dist, prevalence, SAR, R0, R0_social_dist, d0, fatality_pct, FNR=0):
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

        self.infectious = {}
        self.susceptible = {}
        self.prevalence = prevalence
        self.n_households = n_households
        self.R0 = R0
        self.R0_social_dist = R0_social_dist
        self.d0 = d0
        self.SAR = SAR
        self.fatality_pct = fatality_pct

        self.FNR= FNR

        self.remaining_days_infected = {}
        self.deaths = set([])
        self.recoveries = set([])
        self.active_cases = set([])
        self.quarantined_households = set([])
        self.quarantined_individuals = set()

        self.socially_distant_households = set()
        self.socially_distant_individuals = set()

        self.total_pop = 0


        for i in range(n_households):
            # generate household size h from household_size_dist
            randomizer = np.random.uniform()
            h = 0
            while randomizer >= 0:
                randomizer -= household_size_dist[h]
                h += 1
            self.total_pop += h
            # compute primary case probability = p*h/(1+SAR*(h-1))
            prob_prim = prevalence*h/(1+SAR*(h-1))
            infectious = [np.random.uniform() < prob_prim]
            if h > 1:
                # if there are >1 members in the household, and there is a primary case,
                # generate secondary cases from Bin(h-1, SAR); otherwise, set everyone to be uninfected
                infectious.extend(np.random.binomial(1, SAR, h-1) * infectious[0] == True)
            
            susceptible = [not infected for infected in infectious]
            
            self.susceptible[i]=(np.array(susceptible))
            self.infectious[i] =(np.array(infectious))

        self.update_infection_days()

    def begin_social_distancing(self):
        for i in range(self.n_households):
            self.socially_distant_households.add(i)
            for j in range(len(self.infectious[i])):
                self.socially_distant_individuals.add((i,j))

    def end_social_distancing(self):
        for i in range(self.n_households):
            self.socially_distant_households.discard(i)
            for j in range(len(self.infectious[i])):
                self.socially_distant_individuals.discard((i,j))

    def un_distance_household(self, i):
        self.socially_distant_households.discard(i)
        for j in range(len(self.infectious[i])):
            self.socially_distant_individuals.discard((i,j))

    def distance_household(self, i):
        self.socially_distant_households.discard(i)
        for j in range(len(self.infectious[i])):
            self.socially_distant_individuals.discard((i,j))
    
    def get_num_social_dist_households(self):
        return len(self.socially_distant_households)

    def get_num_social_dist(self):
        return len(self.socially_distant_individuals)


    def set_FNR(self, FNR):
        self.FNR=FNR

    def update_infection_days(self):
        for i in range(self.n_households):
            for j in range(len(self.infectious[i])):
                if self.infectious[i][j]:
                    if (i,j) not in self.remaining_days_infected:
                        self.remaining_days_infected[(i,j)] = self.d0
                        self.active_cases.add((i,j))
                    elif self.remaining_days_infected[(i,j)] > 0:
                        self.remaining_days_infected[(i,j)] -= 1
                        if self.remaining_days_infected[(i,j)] == 0:
                            if np.random.uniform() < self.fatality_pct:
                                self.deaths.add((i,j))
                            else:
                                self.recoveries.add((i,j))
                            self.active_cases.discard((i,j))
                            self.susceptible[i][j] = False
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
        infectious_unquarantined = len(self.active_cases - self.quarantined_individuals)
        total_unquarantined = self.total_pop - len(self.quarantined_individuals)

        if total_unquarantined == 0:
            return 0
        else:
            return infectious_unquarantined / total_unquarantined

    def step(self):
        # Simulate one step forward in time
        # Simulate how infectious individuals infect each other
        # Unquarantined susceptible people become infected w/ probability = alpha*current prevalence
        prevalence = self.get_prevalence()
        for i in range(self.n_households):
            household_infected = any(self.infectious[i])
            household_socially_distant = (i in self.socially_distant_households)

            for j in range(len(self.infectious[i])):
                if household_infected and self.susceptible[i][j]:
                    secondary_prob = self.SAR # maybe should depend on total # of infected in household, but this is a start
                    new_secondary = np.random.uniform() < secondary_prob
                else:
                    new_secondary = False

                if not (i,j) in self.quarantined_individuals and self.susceptible[i][j]:
                    if household_socially_distant:
                        primary_prob = np.log(self.R0) * prevalence / self.d0
                    else:
                        primary_prob = np.log(self.R0) * prevalence / self.d0
                    new_primary = np.random.uniform() < primary_prob
                else:
                    new_primary = False

                self.infectious[i][j] = new_primary or new_secondary or self.infectious[i][j]

        self.update_infection_days()

    def get_population_size(self):
        return self.total_pop
    
    def get_num_households(self):
        return self.n_households

    def get_avg_household_size(self):
        return np.mean([len(household) for household in self.infectious.values()])
    
    def iter_individuals(self):
        for i in range(self.n_households):
            for j in range(len(self.infectious[i])):
                yield (i,j)

    def iter_households(self):
        for i in range(self.n_households):
            yield i

    def iter_household_individuals(self, household):
        for j in range(len(self.infectious[household])):
            yield (household, j)

    def test_group(self, grp_individuals):
        infected_subgrp = self.active_cases.intersection(set(grp_individuals))
        if len(infected_subgrp) == 0:
            return False
        else:
            false_neg_prob = self.FNR ** len(infected_subgrp)
            if np.random.uniform() < false_neg_prob:
                return False
            else:
                return True

    def any_infectious(self, grp_individuals):
        return any([self.infectious[i][j] for (i,j) in grp_individuals])

    def react_to_test(self, test_results, groups):
        # unquarantine all confirmed negative results who do not also have
        # an infected family member
        household_risk_status = [False] * self.n_households
        for group_idx, result in test_results.items():
            # for any positive group, set all household risk statuses to True
            # within that group
            if result:
                for i,_ in groups[group_idx]:
                    household_risk_status[i] = True

        # unquarantine all households with negative risk status
        for i, risk_status in enumerate(household_risk_status):
            if risk_status:
                self.quarantine_household(i)
            else:
                self.unquarantine_household(i)

    def quarantine_household(self, i):
        for j in range(len(self.infectious[i])):
            self.quarantined_individuals.add((i,j))
        self.quarantined_households.add(i)

    def unquarantine_household(self, i):
        for j in range(len(self.infectious[i])):
            self.quarantined_individuals.discard((i,j))
        self.quarantined_households.discard(i)


    def quarantine(self, i, j):
        # Put into quarantine all individuals whose indices are in the list x
        if (i,j) not in self.deaths:
            self.quarantined[i][j] = True

    def unquarantine(self, i, j):
        # Remove from quarantine all individuals whose indices are in the list x
        self.quarantined[i][j] = False

    def get_num_infected(self):
        return len(self.active_cases)

    def get_total_dead(self):
        return len(self.deaths)

    def get_total_recovered(self):
        return len(self.recoveries)

    def get_num_quarantined(self):
        return len(self.quarantined_individuals)
    
    def get_num_households_quarantined(self):
        return len(self.quarantined_households)
    #def num_recovered_dead(self):
    # for now we haven't taken into consideration the duration of infection and recovery yet
