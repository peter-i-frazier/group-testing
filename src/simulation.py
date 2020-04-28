"""
Implements generic parent class for running simulations using a Population object
and a TestProtocol object

A Simulation object is responsible for the following:
    * simulating the day-by-day population dynamics by calling the appropriate
      methods in the population object
    * executing the testing protocol at a frequency that is configurable for each instance
    * reporting test results using a test latency that is configurable for each instance
    * recording day-by-day summary data about the disease progression

BaseSimulation has the following attributes:
    * population: an object that inherits BasePopulation
    * test_protocol: an object that inherits BaseTestProtocol
    * test_frequency: an int specifying how frequently we run tests
    * test_latency: an int specifying how many days latency between time-of-test
                    and time-of-test-results
    * test_results: a dict mapping day_test_was_run to test_results
    * current_day: integer specifying how many days since the simulation started
    * recorded_summary_data: dict mapping days to summary data obtained from the population
                             get_summary_data() method

BaseSimulation has the following methods:
    * step(): run one day of the simulation
    * run_test_today(): check whether we run a test today
    * report_test_today(day_test_was_run): check whether we report test results
        today, given the current day and the day the test in question was run

"""

class BaseSimulation:

    def __init__(self, population, 
                        test_protocol, 
                        test_frequency,
                        test_latency):
        self.population = population
        self.test_protocol = test_protocol
        self.test_frequency = test_frequency
        self.test_latency = test_latency

        self.current_day = 0
        self.day_of_last_test = -1
        self.recorded_test_results = {}

        self.summary_population_data = {-1: population.get_summary_data()}
        self.summary_test_data = {-1: test_protocol.get_summary_data()}

    def step(self):
        # begin simulation step() by running one iteration of disease progression on population
        self.population.step()

        # check if today is a testing day; if so run a test
        if self.run_test_today():
            test_results = self.test_protocol.run_test(self.population)
            self.recorded_test_results[self.current_day] = test_results
            self.day_of_last_test = self.current_day


        # check if we have outstanding test results which get reported today;
        # if so report the test
        outstanding_test_days = set(self.recorded_test_results.keys())
        for test_day in outstanding_test_days:
            if self.report_test_today(test_day):
                test_results = self.recorded_test_results[test_day]
                self.test_protocol.respond_to_test(self.population, test_results, test_day)

                # remove test_day from list of outstanding test results
                self.recorded_test_results.pop(test_day)

        # track infection status and quarantine status
        self.population.track_infection_status()
        self.population.track_quarantine_status()

        # record data
        self.summary_population_data[self.current_day] = self.population.get_summary_data()
        self.summary_test_data[self.current_day] = self.test_protocol.get_summary_data()

        self.current_day += 1
    
    
    def run_test_today(self):
        return self.current_day == 0 or \
            self.current_day - self.day_of_last_test >= self.test_frequency
    
    def report_test_today(self, day_test_was_run):
        return self.current_day - day_test_was_run >= self.test_latency 
