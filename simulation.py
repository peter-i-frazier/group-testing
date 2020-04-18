class Simulation:

    def __init__(self, population, 
                        group_test, 
                        test_frequency):
        self.population = population
        self.group_test = group_test
        self.test_frequency = test_frequency

        self.current_day = 0
        self.last_test_day = 0

        self.recorded_data = {}

        self.cumulative_tests_to_date = 0

    def step(self):
        if self.current_day == 0 or self.current_day - self.last_test_day >= self.test_frequency:


            test_results, number_of_tests = self.group_test.test(self.population)

            self.last_test_day = self.current_day

            # enact test results. first unquarantine all negative households
            for (i,j), test_detected_presence in test_results.items():
                if not test_detected_presence:
                    self.population.unquarantine_household(i)

            # next quarantine all positive homes (quarantine takes precedence)
            for (i,j), test_detected_presence in test_results.items():
                if test_detected_presence:
                    self.population.quarantine_household(i)
        else:
            number_of_tests  = 0

        self.cumulative_tests_to_date += number_of_tests

        self.population.step()

        population_size = len(self.population.population)

        sim_data = {
                'in_quarantine_fraction': len(self.population.quarantined_individuals) / population_size,
                'fatality_fraction': len(self.population.fatality_individuals) / population_size,
                'cumulative_infected_fraction': len(self.population.cumulative_infected_individuals) / 
                                                                                    population_size,
                'cumulative_tests_to_date': self.cumulative_tests_to_date
        }

        self.recorded_data[self.current_day] = sim_data

        self.current_day += 1

        return sim_data



