from math import ceil
import numpy as np
from group_testing import HouseholdGroupTest
from frequency import frequency

def static_simulation(group_test):
    # FILL IN MEANINGFUL DETAILS
    QFNR = 0.1
    QFPR = 0.1
    tests_per_person_per_run = 0.1
    quarantines_per_person = 0.4
    return QFNR, QFPR, tests_per_person_per_run, quarantines_per_person

def optimize_gollier_group_size(R0, d0, beta, population_size, FNR, tests_per_week_ub):
    group_size_min = 10
    group_size_max = 1000
    grid_size = 100

    best_group_size = -1
    best_quarantines_per_person = 1.01
    group_size_grid = np.linspace(group_size_min, group_size_max, grid_size)

    for group_size in group_size_grid:
        group_test = HouseholdGroupTest(group_size, 1, FNR)
        QFNR, QFPR, tests_per_person, quarantines_per_person = static_simulation(group_test)

        test_freq = frequency(QFNR, R0, d0, beta)

        freq_per_week = test_freq / 7

        tests_per_week = tests_per_person * population_size * freq_per_week
        
        if tests_per_week <= tests_per_week_ub:
            if quarantines_per_person < best_quarantines_per_person:
                best_quarantines_per_person = quarantines_per_person
                best_group_size = group_size

    return best_group_size, best_quarantines_per_person
    

