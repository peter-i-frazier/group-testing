from math import ceil
from random import shuffle

class GroupTesting:
    def __init__(self, pop_sample, infection_pct=0.01):
        self.pop_sample = pop_sample
        self.infection_pct = infection_pct
        self.sample_size = len(self.pop_sample)

    def gollier_exact(self):
        grp_size = 1 / self.infection_pct

        num_grps = int(ceil(self.sample_size / float(grp_size)))
        
        groups = {}
        for grp_num in range(num_grps):
            groups[grp_num] = []

        individuals = self.pop_sample.keys()
        shuffle(individuals)

        for idx, counter in zip(individuals, range(len(individuals))):
            grp = counter % num_grps
            groups[grp].append(idx)

        test_results = {}
        for grp_idxs in groups.values():
            result = any([self.pop_sample[idx] for idx in grp_idxs])
            for idx in grp_idxs:
                test_results[idx] = result
        
        return test_results, num_grps


