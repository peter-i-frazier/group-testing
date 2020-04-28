import numpy as np

from group_testing import GroupTesting

class DummyPopulation:
    def __init__(self, pop_size=100000, infection_pct=0.01):
        self.pop_size = pop_size
        self.infection_status = {}
        for idx in range(self.pop_size):
            self.infection_status[idx] = (np.random.uniform() <= infection_pct)

class SampleCollection:
    def __init__(self):
        pass
    
    def collect_sample(self, population, sample_pct):
        sample = {}
        for idx, status in population.infection_status.iteritems():
            if np.random.uniform() <= sample_pct:
                sample[idx] = status
        return sample

if __name__ == '__main__':
    population = DummyPopulation()
    sample_collector = SampleCollection()
    sample = sample_collector.collect_sample(population, 0.9)
    group_tester = GroupTesting(sample)
    test_results, num_grps = group_tester.gollier_exact()
    num_safe = len([result for result in test_results.values() if not result])
    print("Exact Gollier identified {} safe individuals out of a sample size {}".format(num_safe,
                                                                                    len(test_results)))
    print("A total of {} tests were used".format(num_grps)) 
