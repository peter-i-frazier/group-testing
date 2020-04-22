import numpy as np
import random
from math import log
import copy

class Population:

    """ This class describes a population of people and their infection status, optionally including information about their
    household organization and the ability to simulate infection status in an correlated way across households.
    It has methods for simulating infection status forward in time """
    
    def __init__(self, n_households, 
                    household_size_dist, 
                    target_prevalence, 
                    disease_length,
                    time_until_symptomatic,
                    non_quarantine_alpha,
                    daily_secondary_attack_rate,
                    fatality_pct,
                    daily_outside_infection_pct,
                    outside_symptomatic_prob,
                    initial_quarantine,
                    initial_prevalence=None,
                    eps = 0.01
                   ):
        # Initialize a population with non-trivial households
        # n_households:         the number of households in the population
        # household_size_dist:  a numpy array that should sum to 1, where household_size_dist[i] gives the fraction of
        #                       households with size i+1
        # prevalence:           prevalence of the disease in the population, used to simulate initial infection status
        # SAR:                  secondary attack rate, used for simulating initial infection status and for simulating
        #                       forward in time; SAR is defined as the probability that an infection occurs among
        #                       susceptible people within a household related to a confirmed case
        #

        self.n_households = n_households
        self.household_size_dist = household_size_dist
        self.target_prevalence = target_prevalence
        self.daily_secondary_attack_rate = daily_secondary_attack_rate
        self.non_quarantine_alpha = non_quarantine_alpha
        self.disease_length = disease_length
        self.time_until_symptomatic = time_until_symptomatic
        self.fatality_pct = fatality_pct
        self.daily_outside_infection_pct = daily_outside_infection_pct
        self.outside_symptomatic_prob = outside_symptomatic_prob
        self.initial_quarantine = initial_quarantine

        assert(np.isclose(sum(self.household_size_dist), 1))
        
        if initial_prevalence == None:
            self._compute_initial_infection_probability()
        else:
            self.initial_prevalence = initial_prevalence
            
        # Reset the population and create an initial infection
        self.reset()

    def sanity_check_prevalences(self):
        prevalences = []
        for _ in range(100):
            self.reset()
            prevalences.append(self.get_num_infected() / float(self.get_num_individuals()))
        if abs(np.mean(prevalences)-self.target_prevalence) > eps:
            print("Warning: Monte-carlo estimate of gap between desired target_prevalence and \
                    Expected target prevalence is large: {}".format(
                        abs(np.mean(prevalences)-self.target_prevalence) ))
            return False

    def _evaluate_r_monte_carlo(self, r, nreps=100):
        prevalences = []
        self.initial_prevalence = r
        for _ in range(nreps):
            self.reset()
            prevalences.append(self.get_num_infected() / float(self.get_num_individuals()))
        return np.mean(prevalences)

    def _evaluate_r_analytic(self, r):

        expected_infections = 0
        expected_population = 0
        SAR = self.daily_secondary_attack_rate
        for i in range(len(self.household_size_dist)):
            house_size = i+1
            prob = self.household_size_dist[i]
            expected_infections = house_size * (r + (1 - r) * SAR - SAR * (1 - r) ** house_size)

            n_houses_of_size = prob * self.n_households
            expected_infections += n_houses_of_size * expected_infections

            expected_population += n_houses_of_size * house_size

        return expected_infections / expected_population



    def _compute_initial_infection_probability(self, eps=0.01):
        # Find initial infection probability r under which E[infection prevalence(r,SAR)] = target_prevalence
        # do binary search over r using Monte Carlo estimates of E[infection prevalence(r,SAR)]
        # nreps: number of samples to use for Monte Carlo estimates
        # eps: return r lying in range r* +/- eps
        r_min = 0
        r_max = self.target_prevalence

        r_prev = 0
        r = 0.5 * (r_min + r_max)
        while abs(r - r_prev) > eps * self.target_prevalence:
            r_prev = r

            value = self._evaluate_r_monte_carlo(r) - self.target_prevalence
            #value = self._evaluate_r_analytic(r) - self.target_prevalence

            # if E[infection_prevalence[(r,SAR)] >= target_prevalence, we have to decrease r
            if value >= 0:
                r_max = r
                r = 0.5 * (r_min + r)
            else:
                r_min = r
                r = 0.5 * (r + r_max)

        self.initial_prevalence = r

    def _sample_house_size(self):
        # Randomly select a house size from the self.household_size_dist distribution
        selector = random.random()
        house_size = 0
        while selector >= 0:
            selector -= self.household_size_dist[house_size]
            house_size += 1
        return house_size

    def reset(self):
        # Resets the population to its initial state, including creation of an initial infection
        n_households = self.n_households

        self.households = set(range(n_households))
        self.household_sizes = {i:self._sample_house_size() for i in range(n_households)}

        self.population = set([(i,j) for i in range(n_households) for j in range(self.household_sizes[i])]) 

        if self.initial_quarantine:
            self.quarantined_individuals = self.population.copy()
            self.unquarantined_individuals = set()
        else:
            self.unquarantined_individuals = self.population.copy()
            self.quarantined_individuals = set()

        self.fatality_individuals = set()
        self.recovered_individuals = set()

        self.days_infected_so_far = {(i,j):0 for (i,j) in self.population}
        self.infected_individuals = set([])

        self.cumulative_infected_individuals = set()
        self.infections_from_inside = 0

        self.days_halted = 0
        self.currently_operating = True

        # infect initial population
        for (i,j) in self.population:
            if random.random() < self.initial_prevalence:
                self.infected_individuals.add((i,j))
                self.cumulative_infected_individuals.add((i,j))

        # Simulate secondary infections
        for i in self.households:
            if any([(i,j) in self.infected_individuals for j in range(self.household_sizes[i])]):
                for j in range(self.household_sizes[i]): 
                    if (i,j) not in self.cumulative_infected_individuals:
                        if random.random() < self.daily_secondary_attack_rate:
                            self.infected_individuals.add((i,j))
                            self.cumulative_infected_individuals.add((i,j))

    def halt_operations(self):
        self.currently_operating = False

    def resume_operations(self):
        self.currently_operating = True

    def is_symptomatic(self, individual):
        return self.days_infected_so_far[individual] >= self.time_until_symptomatic or random.random() < self.outside_symptomatic_prob

    def get_symptomatic_individuals(self):
        symptomatic = set([(i,j) for (i,j) in self.population if self.is_symptomatic((i,j))])
        return symptomatic

    def step(self):
        # Simulate one step forward in time
        # Simulate how infectious individuals infect each other
        # Unquarantined susceptible people become infected w/ probability = alpha*current prevalence
        # TODO: is this the right use of alpha?  Shouldn't it be (1-alpha)*current prevalence?
        # TODO: That way, prevalence(day t+1) = alpha * prevalence(day t)


        # First simulate new primary cases
        current_prevalence = self.get_current_prevalence()

        if self.currently_operating:
            probability_new_infection = log(self.non_quarantine_alpha) * current_prevalence
        else:
            probability_new_infection = 0
            self.days_halted += 1

        for (i,j) in self.unquarantined_individuals:
            if (i,j) not in self.cumulative_infected_individuals:
                if random.random() < probability_new_infection:
                    self.infections_from_inside += 1
                    self.infected_individuals.add((i,j))
                    self.cumulative_infected_individuals.add((i,j))

                elif random.random() < self.daily_outside_infection_pct:
                    self.infected_individuals.add((i,j))
                    self.cumulative_infected_individuals.add((i,j))

        # Next simulate secondary cases
        for i in self.households:
            if any([(i,j) in self.infected_individuals for j in range(self.household_sizes[i])]):
                for j in range(self.household_sizes[i]): 
                    if (i,j) not in self.cumulative_infected_individuals:
                        if random.random() < self.daily_secondary_attack_rate:
                            self.infected_individuals.add((i,j))
                            self.cumulative_infected_individuals.add((i,j))

        individuals_to_resolve = set()
        # update infection counts
        for (i,j) in self.infected_individuals:
            self.days_infected_so_far[(i,j)] += 1

            # see if the disease has lasted its course 
            if self.days_infected_so_far[(i,j)] >= self.disease_length:
                individuals_to_resolve.add((i,j))

        for (i,j) in individuals_to_resolve:

                self.infected_individuals.discard((i,j)) 

                if random.random() < self.fatality_pct:
                    self.fatality_individuals.add((i,j))
                else:
                    self.recovered_individuals.add((i,j))


    def get_avg_household_size(self):
        return np.mean([self.household_sizes[i] for i in range(self.n_households)])


    def get_num_unquarantined(self):
        # Returns the number of people not in quarantine
        return len(self.unquarantined_individuals)

    def get_num_infected(self):
        return len(self.infected_individuals)

    def get_num_individuals(self):
        return len(self.population)

    def get_current_prevalence(self):
        # This returns prevalence in the unquarantined population
        total_unquarantined_infected = len(self.infected_individuals.intersection(self.unquarantined_individuals))
        total_unquarantined = len(self.unquarantined_individuals)

        if total_unquarantined == 0:
            return 0
        else:
            return total_unquarantined_infected / total_unquarantined

    def quarantine_household(self, i):
        for j in range(self.household_sizes[i]):
            self.unquarantined_individuals.discard((i,j))
            self.quarantined_individuals.add((i,j))


    def unquarantine_household(self, i):
        for j in range(self.household_sizes[i]):
            self.quarantined_individuals.discard((i,j))
            self.unquarantined_individuals.add((i,j))

   
