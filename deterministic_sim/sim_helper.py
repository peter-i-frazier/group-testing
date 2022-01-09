import numpy as np
from testing_regime import TestingRegime
from plotting import Trajectory
from strategy import Strategy
from sim import sim
from groups import population
from typing import Dict


# TODO (hwr26): popul is a redundant parameter. Try to reorganize to prevent.
def sim_test_strategy(scenario: Dict, popul: population, strategy: Strategy,
    color: str) -> Trajectory:
    """Return a trajectory in [color] representing a [strategy] on [scenario]
    with a given [popul].

    Args:
        scenario (Dict): Scenario on which the simulation is run.
        popul (population): The population of the simulation.
        strategy (Strategy): The strategy to be used with the simulation.
        color (str): The color of the trajectory.
    """
    T = scenario['T']
    GENERATION_TIME = scenario['generation_time']
    BOOSTER_EFFECTIVENESS = scenario['booster_effectiveness']
    INFECTIONS_PER_DAY_PER_CONTACT_UNIT = \
        np.array(list(scenario['infections_per_day_per_contact_unit'].values()))

    for i in range(strategy.periods):
        regime = strategy.testing_regimes[i]
        infections_per_contact_unit = BOOSTER_EFFECTIVENESS * \
                                        np.multiply(strategy.transmission_multipliers[i], \
                                                    np.multiply(
                                                        INFECTIONS_PER_DAY_PER_CONTACT_UNIT, \
                                                        regime.get_days_infectious()))
        infection_rate = popul.infection_matrix(infections_per_contact_unit)
        infection_discovery_frac = popul.metagroup2group(regime.get_infection_discovery_frac())
        recovered_discovery_frac = popul.metagroup2group(regime.get_recovered_discovery_frac())

        initial_infections = strategy.get_initial_infections(scenario)
        past_infections = strategy.get_past_infections(scenario)
        S0, I0, R0 = popul.get_init_SIR_vec(initial_infections, past_infections,
                                            weight="population x contacts")
        outside_rates = list(scenario['outside_rates'].values())
        outside_rate = popul.get_outside_rate(outside_rates)

        if i == 0: # instantiate simulation object
            s = sim(T, S0, I0, R0, infection_rate,
                    infection_discovery_frac=infection_discovery_frac ,
                    recovered_discovery_frac=recovered_discovery_frac,
                    generation_time=GENERATION_TIME, outside_rate=outside_rate)
        s.step(strategy.period_lengths[i], infection_rate=infection_rate,
            infection_discovery_frac = infection_discovery_frac,
            recovered_discovery_frac = recovered_discovery_frac)

    return Trajectory(strategy, s, color)


def sim_test_regime(scenario: Dict, popul: population, tests_per_week: int,
    delay: float, color: str):
    """Simulate a testing regime with no pre-departure or arrival testing."""

    CLASSWORK_TRANSMISSION_MULTIPLIER = \
        list(scenario['classwork_transmission_multiplier'].values())
    regime = TestingRegime(popul,tests_per_week,delay,scenario)

    strategy = \
        Strategy(name=regime.name,
            # TODO (hwr26): This matches what we had before but we should
            # choose a scenario to run with
            pct_discovered_in_pre_departure=0.25,   # Using Scenario 2
            pct_discovered_in_arrival_test=0.38,  # Using Scenario 2
            testing_regimes=[regime, regime],
            transmission_multipliers=[1, CLASSWORK_TRANSMISSION_MULTIPLIER],
            period_lengths=[3,scenario["T"]-3-1])

    return sim_test_strategy(scenario, popul, strategy, color)
