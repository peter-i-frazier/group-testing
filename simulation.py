class Simulation:

    def __init__(self, population, grptest, test_freq):
        self.population = population
        self.grptest = grptest
        self.test_freq = test_freq
        self.current_day = 0
        self.last_test_day = 0

    def step(self):
        if self.current_day == 0 or self.current_day - self.last_test_day >= self.test_freq:
            test_results, groups = self.grptest.test(self.population)
            self.population.update_quarantine_status(test_results, groups)
            self.last_test_day = self.current_day

        self.current_day += 1
        self.population.step()


