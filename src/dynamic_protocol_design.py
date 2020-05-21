from eval_r import eval_r, match_r, US_dist
from population import Population
from group_testing import HouseholdGroupTest, MatrixGroupTest
from static_simulation import StaticSimulation
import numpy as np
import sys

def test_properties(prevalence, group_size, FNR = 0.3, FPR = 0.1, n_households = 1000, nreps = 10):
    # To make the code accurate, use n_households = 22500 and n_reps = 100
    # To make it faster use the default parameters
    initial_prevalence = eval_r(match_r, prevalence)
    pop = Population(n_households=n_households, # Should be big relative to the largest group size
                      household_size_dist=US_dist,
                      target_prevalence=prevalence,
                      disease_length=0,
                      time_until_symptomatic=0,
                      non_quarantine_alpha=0, # not needed for present purposes
                      daily_secondary_attack_rate=0.374,
                      fatality_pct=0,
                      daily_outside_infection_pct=0,
                      outside_symptomatic_prob=0,
                      initial_quarantine=0,
                      initial_prevalence=initial_prevalence)

    group_test = MatrixGroupTest(group_size, FNR, FPR, fnr_at_swab_level=True)
    QFNR, QFPR, tests_per_person, quarantines_per_person = StaticSimulation(pop,group_test).sim(nreps)

    return QFNR, QFPR, tests_per_person

def __test_test_properties():
    print(test_properties(0.01, 100))
    print('OK')


def simple_simulation(group_sizes=[80,90,100,110], prevalence=0.01, resolved=0., seed=1, verbose=True):
    # This code has been modified to support Raul's 2020 NeurIPS submission with Peter
    # From the perspective of his paper, this code takes as input the decision variables for each node,
    # which are the group sizes, i.e., group_sizes[i] is the input for node i.
    # It returns as output, for each node:
    #   the output state, which is a vector of occupancies
    #   the loss
    # The total loss is then the sum of the losses across all of the nodes

    np.random.set_state(seed)

    T = len(group_sizes)


    # state[t] is the input to node t
    # group_sizes[t] is the control for node t
    # state[t+1] is the output from node t
    # loss[t] is the loss that results from node t
    state = [{} for i in range(T+1)]
    loss = [0 for i in range(T)]

    # Construct the initial state
    state[0] = { 'FI': prevalence,  # free infected
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
        if verbose:
            print('Week {}'.format(t))
            print_state(state[t], 'Start of Week, Before Test:')

        state[t+1], loss[t] = node(state[t], group_sizes[t],verbose=verbose)

        if verbose:
            print_state(state[t], 'End of Week:               ')

    # When we return the state, we drop state[T] and only return 0 through T-1 inclusive
    # loss has already been prepared so that it has only has components 0 through T-1
    return state[0:T], loss


def node(state, group_size, costs = {'test' : .001, 'infections' : 1, 'quarantine' : .01}, verbose=True):
    # Takes as input:
    #   the state of the system (which is the output from the previous node, or the initial state)
    #   the group size (which is the control to exercise in that period)
    # Whether to provide verbose output
    # Returns the output state (a dictionary) and a loss
    #
    # Interesting question from Raul: is the loss a linear function of the state?
    # Answer: no, but it's probably pretty simple to model, and some components are linear
    # The loss is a linear combination of the PCR capacity used,
    # (Polymerase chain reaction, chemical test that tells you whether )
    # the number of people in quarantine at the end of the week,
    # and the number of new infections

    loss = 0

    # Set some epidemiological parameters
    doubling_time = 3.0
    alpha = 2 ** (1 / doubling_time)
    test_period = 7 # days between tests

    # Assumes that we don't know that recovered people are recovered, and so we test them
    size_of_tested_pop = state['FI'] + state['FS'] + state['FR']
    prevalence_in_tested_pop = state['FI'] / size_of_tested_pop

    QFNR, QFPR, tests_per_person = test_properties(prevalence_in_tested_pop, group_size)

    tests_per_person_overall = tests_per_person * size_of_tested_pop

    # We add a component to the loss which is the number of tests performed times the cost per test
    loss = tests_per_person_overall * costs['test']

    if verbose:
        print('Test Properties: Group Size={}, PCR/person={:.3f} QFNR={:.1f}% QFPR={:.1f}%'.format(group_size,
            tests_per_person_overall,
            QFNR * 100,
            QFPR * 100))

    # Simulate the effect of the test
    old_state = state.copy() # Remember previous state so we can refer back after modifying entries

    state['FI'] = old_state['FI'] * QFNR  # the only people who remain free infected are those that are missed in the test
    state['FS'] = old_state['FS'] * (1 - QFPR) + old_state['QS2']
    state['FR'] = old_state['FR'] * (1 - QFPR) + old_state['QI2'] + old_state['QR2']  # Assumes all infected leaving quarantine are now resolved

    state['QI1'] = old_state['FI'] * (1 - QFNR)
    state['QS1'] = old_state['FS'] * QFPR
    state['QR1'] = old_state['FR'] * QFPR

    state['QI2'] = old_state['QI1']
    state['QS2'] = old_state['QS1']
    state['QR2'] = old_state['QR1']

    # Add the loss from having people in quarantine
    people_in_quarantine  = state['QI1'] + state['QS1'] + state['QR1']
    people_in_quarantine += state['QI2'] + state['QS2'] + state['QR2']
    loss += people_in_quarantine * costs['quarantine']

    if verbose:
        print_state(state, 'Start of Week, After test: ')


    # Simulate the growth of the disease in the free population
    fraction_susceptible = 1
    growth = alpha ** (test_period * fraction_susceptible)
    new_infections = state['FI'] * (growth - 1)
    state['FI'] = state['FI'] + new_infections
    state['FS'] = state['FS'] - new_infections

    loss += new_infections * costs['infections']

    return state, loss


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
    prevalence = float(sys.argv[1])
    group_sizes = [int(sys.argv[idx]) for idx in range(2,6)]
    print("Beginning simulation with prevalence {} and group sizes {}".format(prevalence, group_sizes))
    state, loss = simple_simulation(group_sizes, prevalence)
    print(state)
    print(loss)
