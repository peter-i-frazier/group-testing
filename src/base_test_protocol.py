"""
Implements a generic glass for a test protocol object

A test protocol is responsible for running a test on a BasePopulation object
and reporting a standardized test result 

BaseTestProtocol has the following attributes:
    * cumulative_num_tests

BaseTestProtocol has the following methods:
    * run_test(population)
    * get_summary_data()

"""


class BaseTestProtocol:


    def __init__(self):
        self.cumulative_num_tests = 0


    def run_test(self, population):
        raise(Exception("run_test() must be implented by child class"))

    def respond_to_test(self, population, test_results, test_day):
        raise(Exception("respond_to_test() must be implemented by child class"))

    def get_summary_data(self):
        return {'cumulative_num_tests':self.cumulative_num_tests}

    def update_test_count(self, num_tests):
        self.cumulative_num_tests += num_tests


class EmptyTestProtocol(BaseTestProtocol):

    def run_test(self, population):
        return {}

    def respond_to_test(self, population, test_results, test_day):
        pass


class QuarantineSymptomaticProtocol(BaseTestProtocol):

    def run_test(self, population):
        test_results = set()

        for agent_id in population.agents:
            if population.is_agent_symptomatic(agent_id):
                test_results.add(agent_id)

        return test_results
    

    def respond_to_test(self, population, test_results, test_day):
        for agent_id in test_results:
            population.quarantine_agent(agent_id)

class ContactTraceProtocol(BaseTestProtocol):

    def __init__(self):
        super().__init__()
        self.followup_contacts = set()

    def run_test(self, population):
        num_tests = 0
        test_results = set()

        for agent_id in population.iter_unquarantined_symptomatic():
            num_tests += 1
            test_results.add(agent_id)


        for agent_id in self.followup_contacts:
            if population.is_agent_infected(agent_id):
                test_results.add(agent_id)
                num_tests +=  1

        recent_contacts = population.get_recent_contacts(test_results)
        new_followup = recent_contacts.difference(test_results)
       
        self.followup_contacts = new_followup

        return test_results

    def respond_to_test(self, population, test_results, test_day):
        for agent_id in test_results:
            population.quarantine_agent(agent_id)

