class Simulation:

    def __init__(self, population, 
                        group_test, 
                        test_frequency,
                        test_latency,
                        halt_operations_if_case_detected):
        self.population = population
        self.group_test = group_test
        self.test_frequency = test_frequency

        self.current_day = 0
        self.last_test_day = 0

        self.recorded_data = {}

        self.cumulative_tests_to_date = 0
        self.test_latency = test_latency
        self.has_recent_test_been_reported = False
        self.test_results = None
        self.halt_operations_if_case_detected = halt_operations_if_case_detected

    def report_test(self):
        assert(self.test_results != None and not self.has_recent_test_been_reported)
        self.has_recent_test_been_reported = True

        if self.halt_operations_if_case_detected:
            any_detected = False
            for (i,j) in self.test_results:
                if (i,j) in self.population.unquarantined_individuals:
                    if self.test_results[(i,j)]:
                        any_detected = True
            if any_detected:
                self.population.halt_operations()
            else:
                self.population.resume_operations()

         # enact test results. first unquarantine all negative households
        for (i,j), test_detected_presence in self.test_results.items():
            if not test_detected_presence:
                self.population.unquarantine_household(i)

        # next quarantine all positive homes (quarantine takes precedence)
        for (i,j), test_detected_presence in self.test_results.items():
            if test_detected_presence:
                self.population.quarantine_household(i)

    def step(self):
        if self.current_day == 0 or self.current_day - self.last_test_day >= self.test_frequency:


            self.test_results, number_of_tests = self.group_test.test(self.population)
            self.has_recent_test_been_reported = False


            self.last_test_day = self.current_day

        else:
            number_of_tests  = 0

        if self.current_day - self.last_test_day >= self.test_latency and not self.has_recent_test_been_reported:
            self.report_test()
            self.has_recent_test_been_reported = True

        self.cumulative_tests_to_date += number_of_tests

        self.population.step()

        population_size = len(self.population.population)

        sim_data = {
                'in_quarantine_fraction': len(self.population.quarantined_individuals) / population_size,
                'infected_unquarantined_fraction': len(self.population.unquarantined_individuals.intersection(self.population.infected_individuals)) / population_size,
                'infected_fraction': len(self.population.infected_individuals) / population_size,
                'fatality_fraction': len(self.population.fatality_individuals) / population_size,
                'cumulative_infected_fraction': len(self.population.cumulative_infected_individuals) / 
                                                                                    population_size,
                'cumulative_infected_within_population': self.population.infections_from_inside / 
                                                                                    population_size,
                'cumulative_tests_to_date': self.cumulative_tests_to_date,
                'cumulative_days_halted': self.population.days_halted
        }

        self.recorded_data[self.current_day] = sim_data

        self.current_day += 1

        return sim_data



