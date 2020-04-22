from eval_r import eval_r, match_r, US_dist
from population import Population
from group_testing import HouseholdGroupTest, MatrixGroupTest
from static_simulation import StaticSimulation

def test_properties(prevalence, group_size, FNR = 0.3, FPR = 0.1):
    initial_prevalence = eval_r(match_r, prevalence)

    #doubling_time = 3.0
    #alpha = 2 ** (1 / doubling_time)
    alpha = 0 # I don't think we need this
    SAR = 0.374

    nreps = 100

    pop = Population(n_households=22500, # Should be big relative to the largest group size
                      household_size_dist=US_dist,
                      target_prevalence=prevalence,
                      disease_length=0,
                      time_until_symptomatic=0,
                      non_quarantine_alpha=alpha,
                      daily_secondary_attack_rate=SAR,
                      fatality_pct=0,
                      daily_outside_infection_pct=0,
                      outside_symptomatic_prob=0,
                      initial_quarantine=0,
                      initial_prevalence=initial_prevalence)

    # group_test = HouseholdGroupTest(group_size, 1, FNR, FPR)
    group_test = MatrixGroupTest(group_size, FNR, FPR, fnr_at_swab_level=False)
    QFNR, QFPR, tests_per_person, quarantines_per_person = StaticSimulation(pop,group_test).sim(nreps)

    return QFNR, QFPR, tests_per_person, quarantines_per_person


if __name__ == '__main__':
    print(test_properties(0.01, 100))
    print('OK')