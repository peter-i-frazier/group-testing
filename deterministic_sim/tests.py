import numpy as np
from sim import sim, well_mixed_infection_rate
from micro import __days_infectious_perfect_sensitivity__, days_infectious


def is_constant_population_size(s):
    """Checks that the population size of each group is contant across sim."""
    pop = s.get_metric("S") + s.get_metric("I") + s.get_metric("R")
    return all(np.isclose(pop, np.max(pop)))


def test_sim_symmetric_populations():
    """Test spread through symmetric populations is spread."""
    K = 3
    T = 20

    # symmetric populations
    pop = 100 * np.ones(K)/K
    R0 = 100 * np.array([.1, .1, .1])
    I0 = 100 * np.array([.01, .01, .01])
    S0 = pop - R0 - I0

    # Here, each group has an overall rate at which they infect others.
    # These infections are distributed proportionally to the overall population.
    # infection_rate[i,j] is the number of new infections that an infected
    # person in i creates in j.
    contact_rates = np.array([0, 1, 2])
    infection_rate = np.outer(contact_rates, pop / 100)

    s = sim(T, S0, I0, R0, infection_rate=infection_rate)
    s.step(T-1)

    assert is_constant_population_size(s)
    # S,I,R amounts for each population are initially symmetric and the new
    # infections are distributed proportionally across the populations, the
    # number infected should stay symmetric over time. That is, for each value
    # of t, y[t,0], y[t,1], and y[t,2] should be nearly identical.
    y = s.get_metric('I', aggregate=False, normalize=True)
    assert(np.isclose(y[:,0],y[:,1]).all())
    assert(np.isclose(y[:,0],y[:,2]).all())


def test_noninfectious_group():
    """Test zero cases in a non-infectious group."""
    K = 3
    T = 20

    # symmetric populations
    pop = 100 * np.ones(K)/K
    R0 = 100 * np.array([.1, .1, .1])
    I0 = 100 * np.array([0, .01, .01])
    S0 = pop - R0 - I0

    infections_per_contact = 1
    marginal_contacts = np.array([0,1,2])
    infection_rate = well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)

    s = sim(T, S0, I0, R0, infection_rate=infection_rate)
    s.step(T-1)

    # The group with 0 contacts should not have any infections
    assert(np.isclose(s.get_metric_for_group('I', 0),np.zeros(T)).all())
    assert is_constant_population_size(s)


def test_sim_zero_prob_discvoered():
    """Test zero discovered cases when probability of discovery is zero."""
    K = 3
    T = 20

    # symmetric populations
    pop = 100 * np.ones(K)/K
    R0 = 100 * np.array([.1, .1, .1])
    I0 = 100 * np.array([0, .01, .01])
    S0 = pop - R0 - I0

    infections_per_contact = 1
    marginal_contacts = np.array([0,1,2])
    infection_rate = well_mixed_infection_rate(pop, marginal_contacts, infections_per_contact)

    generation_time = 4/7 # in units of weeks

    s = sim(T, S0, I0, R0, infection_rate=infection_rate,
            generation_time=generation_time, infection_discovery_frac=0,
            recovered_discovery_frac=0)
    s.step(T-1)

    assert is_constant_population_size(s)
    # Since no one is discovered, the number discovered should be 0
    discovered = s.get_discovered(cumulative=True)
    assert(np.isclose(discovered, np.zeros(T)).all())


# TODO (hwr): Leaving a note here that test_sim3 containted unused analysis
# but was not a test. Hence, it was removed on 2021-12-28. Notes given below.

# Calculate the probability of infection.
# A ballpark estimate is that R0 with 1x / wk testing and a 2-day delay from
# sampling to isolation is 5, with a 4 day generation time.  We think that
# this outbreak was driven by the people with 6 contacts per period above
# because the sum of their proportions of the UG population is about 1900
# people. To achieve this, we set infections_per_contact = 5/6 so that the
# number of secondary infections from someone with 6 contacts is 5. We then
# adjust this by multiplying by the number of days infectious under our
# testing strategy, divided by the number under our December Omicron outbreak.


def test_days_infectious():

    # If there is infinite time between tests, the time infectious in the
    # presence of testing should be whatever the length of the maximum infectious
    # period is. This should be true regardless of the sensitivity.
    assert __days_infectious_perfect_sensitivity__(np.inf,1,max_infectious_days = 5) == 5
    assert days_infectious(np.inf,1,max_infectious_days = 5) == 5

    # Time infectious with testing should be less than the maximum
    assert days_infectious(2,1,max_infectious_days = 5) < 5

    # Testing with more delay should result in a longer time infectious
    assert days_infectious(7,1) > days_infectious(2, 1)

    # A shorter max infectioun time should result in fewer days infectious
    assert days_infectious(5, 1, max_infectious_days=5) < days_infectious(5, 1, max_infectious_days=7)

    # Codes should agree when sensitivity is perfect
    assert days_infectious(7,1,sensitivity=1) == __days_infectious_perfect_sensitivity__(7,1)
    assert days_infectious(4,2,sensitivity=1) == __days_infectious_perfect_sensitivity__(4,2)

    # A lower sensitivity should result in more days infectious
    assert days_infectious(5, 1, sensitivity=.5) > days_infectious(5, 1, sensitivity=.7)

    # A sensitivity of 0 should be like not testing
    assert days_infectious(2,1,max_infectious_days = 5, sensitivity = 0) == 5

    # days infectious should be larger when the days between tests is larger
    for d in range(20):
        assert days_infectious(d+2, 1) > days_infectious(d+1, 1)

    # TODO (hwr): commenting out plots created by this test to reduce clutter
    # but leaving commented out in case someone wants to view them later.

    # days_between_tests = np.arange(1, 20, .1)
    # y70 = [days_infectious(d, 1, sensitivity = 0.7) for d in days_between_tests]
    # y90 = [days_infectious(d, 1, sensitivity = 0.9) for d in days_between_tests]
    # y100 = [days_infectious(d, 1, sensitivity = 1) for d in days_between_tests]
    # plt.plot(days_between_tests, y70, label = 'Sensitivity = 70%')
    # plt.plot(days_between_tests, y90, label = 'Sensitivity = 90%')
    # plt.plot(days_between_tests, y100, label = 'Sensitivity = 100%')
    # plt.legend()
    # plt.xlabel('Number of days between tests')
    # plt.ylabel('Expected Days infectious')
    # plt.savefig('test_days_infectious1.png', facecolor='w')
    # plt.close()

    # # Compare what we are using here to what we used before
    # days_between_tests = np.arange(1, 20, .1)
    # y_old = [days_infectious(d, 1, sensitivity=1, max_infectious_days=7) for d in days_between_tests]
    # y_new = [days_infectious(d, 1) for d in days_between_tests]
    # plt.plot(days_between_tests, y_old, 'k--', label='Sensitivity = 100%, max_infectious = 7 days, 1 day delay')
    # plt.plot(days_between_tests, y_new, 'k-', label='Sensitivity = 60%, max_infectious = default, 1 day delay')

    # y_old = [days_infectious(d, 2, sensitivity=1, max_infectious_days=7) for d in days_between_tests]
    # y_new = [days_infectious(d, 2) for d in days_between_tests]
    # plt.plot(days_between_tests, y_old, 'b--', label='Sensitivity = 100%, max_infectious = 7 days, 2 day delay')
    # plt.plot(days_between_tests, y_new, 'b-', label='Sensitivity = 60%, max_infectious = default, 2 day delay')

    # plt.legend()
    # plt.xlabel('Number of days between tests')
    # plt.ylabel('Expected Days infectious')
    # plt.savefig('test_days_infectious2.png', facecolor='w')
    # plt.close()
