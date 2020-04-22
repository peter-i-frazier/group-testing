from eval_r import eval_r, match_r, US_dist
from population import Population
from group_testing import HouseholdGroupTest, MatrixGroupTest
from static_simulation import StaticSimulation
import numpy as np

def test_properties(prevalence, group_size, FNR = 0.3, FPR = 0.1):
    initial_prevalence = eval_r(match_r, prevalence)


    alpha = 0 # I don't think we need this
    SAR = 0.374

    nreps = 1
    #nreps = 100
    # n_households=22500

    pop = Population(n_households=10000, # Should be big relative to the largest group size
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

    return QFNR, QFPR, tests_per_person

def __test_test_properties():
    print(test_properties(0.01, 100))
    print('OK')


def simple_simulation(group_sizes=[80,90,100,110], prevalence=0.01, resolved=0.):

    T = len(group_sizes)

    doubling_time = 3.0
    alpha = 2 ** (1 / doubling_time)
    test_period = 7 # days between tests

    state = { 'FI': prevalence,  # free infected
              'FS': 1. - prevalence - resolved, # free susceptible
              'FR': 0.,  #free resolved
              'QI1' : 0., # infected,    in quarantine for 1 test period
              'QS1' : 0., # susceptible, in quarantine for 1 test period
              'QR1' : 0., # resolved,    in quarantine for 1 test period
              'QI2': 0.,  # same thing, but in quarantine for 2 test periods
              'QS2': 0.,
              'QR2': 0.,
              }

    for t in range(T):
        print()
        print('Week {}'.format(t))

        print_state(state, 'Start of Week, Before Test:')

        # Assumes that we don't know that recovered people are recovered, and so we test them
        size_of_tested_pop = state['FI'] + state['FS'] + state['FR']
        prevalence_in_tested_pop = state['FI'] / size_of_tested_pop

        QFNR, QFPR, tests_per_person = test_properties(prevalence_in_tested_pop, group_sizes[t])

        tests_per_person_overall = tests_per_person * size_of_tested_pop

        print('Test Properties:            Group Size={}, PCR/person={:.3f} QFNR={:.1f}% QFPR={:.1f}%'.format(group_sizes[t], tests_per_person_overall,QFNR*100,QFPR*100))

        # Simulate the effect of the test

        old_state = state.copy() # Remember previous state so that when we modify entries in state we can refer back to what the previous state was

        state['FI'] = old_state['FI']*QFNR # the only people who remain free infected are those that are missed in the test
        state['FS'] = old_state['FS']*(1-QFPR) + old_state['QS2']
        state['FR'] = old_state['FR']*(1-QFPR) + old_state['QI2'] + old_state['QR2'] # Assumes all infected leaving quarantine are now resolved

        state['QI1'] = old_state['FI'] * (1 - QFNR)
        state['QS1'] = old_state['FS'] * QFPR
        state['QR1'] = old_state['FR'] * QFPR

        state['QI2'] = old_state['QI1']
        state['QS2'] = old_state['QS1']
        state['QR2'] = old_state['QR1']

        print_state(state, 'Start of Week, After test: ')

        # Simulate the growth of the disease in the free population
        fraction_susceptible = 1
        growth = alpha**(test_period*fraction_susceptible)
        new_infections = state['FI'] * (growth-1)
        state['FI'] = state['FI'] + new_infections
        state['FS'] = state['FS'] - new_infections
        # No change in resolved

        print_state(state, 'End of Week:               ')

def print_state(state,preamble='',verbose=False):

    free = state['FI'] + state['FS'] + state['FR']
    infected  = state['FI'] + state['QI1'] + state['QI2']
    recovered = state['FR'] + state['QR1'] + state['QR2']
    print('{} Free={:.1f}%, Infected={:.1f}%, Infected&Free={:.1f}%, RecoveredOrDead={:.1f}%'.format(
        preamble,free*100,infected*100,state['FI']*100,recovered*100))

    if verbose:
        print(state)

    assert (np.isclose(sum([p[1] for p in state.items()]), 1.))


if __name__ == '__main__':
    # __test_test_properties():
    # simple_simulation()


