from math import ceil
import numpy as np
from group_testing import HouseholdGroupTest, ThreeStageHierarchicalTest
from frequency import frequency
from population import Population
from static_simulation import StaticSimulation

def optimize_gollier_group_size(population, beta, FNR, FPR, tests_per_person_per_week_ub,
                                group_size_min = 10, group_size_max = 1000, grid_size = 100, nreps = 100,
                                group_test_name='Gollier',
                                verbose=True):
    best_group_size = -1
    best_quarantines_per_person = 1.01
    best_days_between_tests = -1
    best_tests_per_person_per_week = -1
    group_size_grid = np.linspace(group_size_min, group_size_max, grid_size)

    for group_size in group_size_grid:
        if group_test_name == 'Gollier':
            group_test = HouseholdGroupTest(group_size, 1, FNR, FPR)
        elif group_test_name == 'ThreeStageHierarchical':
            group_test = ThreeStageHierarchicalTest(group_size, np.sqrt(group_size), 1, FNR, FNR)
        else:
            print('Unknown group test requested {}.  Using Gollier'.format(group_test_name))
            group_test = HouseholdGroupTest(group_size, 1, FNR, FPR)

        # We want (population size % group size) / (population size) to be small

        QFNR, QFPR, tests_per_person, quarantines_per_person = StaticSimulation(population,group_test).sim(nreps)

        # Number of days between screenings
        days_between_tests = frequency(QFNR, population.non_quarantine_alpha, beta)

        weeks_between_tests = days_between_tests / 7

        tests_per_person_per_week = tests_per_person / weeks_between_tests

        if verbose:
            print('group size={:.2f} QFNR={} QFPR={:.2f} tests/person/wk={:.3f} quarantines/person={:.2f} days between tests={:2f}'.format(
                group_size,QFNR,QFPR,tests_per_person_per_week,quarantines_per_person,days_between_tests))
        
        if tests_per_person_per_week <= tests_per_person_per_week_ub:
            if quarantines_per_person < best_quarantines_per_person:
                best_quarantines_per_person = quarantines_per_person
                best_group_size = group_size
                best_days_between_tests = days_between_tests
                best_tests_per_person_per_week = tests_per_person_per_week

    return best_group_size, best_quarantines_per_person, best_tests_per_person_per_week, best_days_between_tests


def larry_analysis(prevalence=0.01,household_size_dist=[1],group_test_name='Gollier',
                   group_size_min = 10, group_size_max = 100, grid_size = 9, nreps = 100, verbose=True):
    doubling_time = 3. # number of days for epidemic to double in the absence of distancing
    alpha = 2**(1/doubling_time)
    SAR = 0.374

    assert(np.isclose(sum(household_size_dist), 1)) 

    pop = Population(n_households=10000, # Should be big relative to the largest group size
                     household_size_dist=household_size_dist,
                     target_prevalence=prevalence,
                     disease_length=0,
                     time_until_symptomatic=0,
                     non_quarantine_alpha=alpha,
                     daily_secondary_attack_rate=SAR,
                     fatality_pct=0,
                     daily_outside_infection_pct=0,
                     outside_symptomatic_prob=0,
                     initial_quarantine=0)

    if verbose:
        print("Done initializing population.  Using initial prevalence {} for target prevalence {}".format(
            pop.initial_prevalence, pop.target_prevalence))

    beta = 1
    FNR = 0.3
    FPR = 0.1
    uspop = 328E6
    tests_per_week_ub = 6E6

    tests_per_week_per_person_ub = tests_per_week_ub / uspop
    best_group_size, best_quarantines_per_person, test_per_person_per_week, days_between_tests =\
        optimize_gollier_group_size(pop,beta,FNR, FPR,tests_per_week_per_person_ub,
                                    group_size_min, group_size_max, grid_size, nreps,
                                    group_test_name, True)
    print('{} doubling time {} days, household size dist={} SAR={} test FNR={} FPR={} beta={}'.format(
        group_test_name, doubling_time,household_size_dist, SAR,FNR,FPR,beta))
    print('Optimum @ prevalence={}: group size={} quarantines/person={:.2f} tests/week/person={:.3f} days between tests={:2f}'.format(
        prevalence,best_group_size,best_quarantines_per_person,test_per_person_per_week,days_between_tests))

#https://www.statista.com/statistics/242189/disitribution-of-households-in-the-us-by-household-size/ 
us_census_household_dist = [0.2837, 0.3451, 0.1507, 0.1276, 0.0578, 0.0226, 0.0125] 

if __name__ == '__main__':
    # 0.0183 is the threshold on tests/person/wk

    #larry_analysis(0.005,1,'ThreeStageHierarchical', group_size_min = 30, group_size_max = 80, grid_size = 10, nreps = 1000)

    # Optimum @ prevalence = 0.005: group size = 70.00 quarantines / person = 0.21 tests / week / person = 0.017 days between tests = 6.142673
    #larry_analysis(0.005,1,'Gollier', group_size_min = 70, group_size_max = 80, grid_size = 3, nreps = 1000)

    # Optimum @ prevalence = 0.01: group size = 60.00 quarantines / person = 0.34 tests / week / person = 0.017 days between tests = 6.970806
    # Spreadsheet (Appendix 2) confirms this (fraction of pop without a known case that can exit is 65.61%, days between tests is 7)
    larry_analysis(0.005,us_census_household_dist,'Gollier', group_size_min = 70, group_size_max = 70, grid_size = 1, nreps = 1000)
    # At household size = 2, according to spreadsheet (Appendix 3), can do group size = 45, quarantine fraction = 30.5%, 8.19E05 tests / day, 8.9 days between tests

    larry_analysis(0.01,us_census_household_dist,'Gollier',
                    group_size_min = 45, group_size_max = 45, grid_size = 1, nreps = 1000)
    larry_analysis(0.04,us_census_household_dist,'Gollier', group_size_min = 40, group_size_max = 40, grid_size = 1, nreps = 1000)
    # Optimum @ prevalence = 0.04: group size = 40.00 quarantines / person = 0.68 tests / week / person = 0.017 days between tests = 10.065353
    # larry_analysis(0.04,1,'Gollier', group_size_min = 30, group_size_max = 50, grid_size = 5, nreps = 1000)

def yujia_analysis():
    # iterate over prevalence
    # iterate over group sizes
    for p in np.linspace(prevalence_min, prevalence_max, grid_size)




