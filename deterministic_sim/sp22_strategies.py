from strategy import Strategy
from testing_regime import TestingRegime
from typing import Dict

# ========================================================
# Define testing regimes used by potential sp22 strategies
# ========================================================

def no_testing_testing_regime(scenario: Dict):
    return TestingRegime(scenario=scenario, tests_per_week=0, test_delay=1)

def ug_prof_2x_week_testing_regime(scenario: Dict):
    return TestingRegime(scenario=scenario,
                         tests_per_week={ 'UG':2, 'GR':0, 'PR':2, 'FS':0},
                         test_delay=1.5)

# =========================
# Potential sp22 strategies
# =========================

def no_testing_strategy(scenario: Dict):
    T = scenario['T']
    CLASSWORK_TRANSMISSION_MULTIPLIER = \
        list(scenario['classwork_transmission_multiplier'].values())
    return \
        Strategy(name="No Testing",
            pct_discovered_in_pre_departure=0,
            pct_discovered_in_arrival_test=0,
            testing_regimes=[no_testing_testing_regime(scenario),
                             no_testing_testing_regime(scenario)],
            transmission_multipliers=[1, CLASSWORK_TRANSMISSION_MULTIPLIER],
            period_lengths=[3,T-3-1])

# SCENARIO 1:
# all students do a pre-departure antigen test with sensitivity 50%
# all students again an arrival PCR test with sensitivity 67% and no isolation delay
#     0.5 caught in pre-departure
#     (1 - 0.5) * 0.67 = 0.335 caught in arrival testing
#     1 - 0.5 - 0.335 = 0.165 not caught

# SCENARIO 2:
# half of students do a pre-departure test with sensitivity 50%
# 75% of students (independently chosen) do an arrival PCR 67% sensitivity and no test delay
#     0.5 * 0.5 = 0.25 caught in pre-departure
#     (1 - 0.25) * (0.75 * 0.67) = 0.38 caught in arrival testing
#     1 - 0.25 - 0.38 = 0.37 not caught

def arrival_testing_strategy(scenario: Dict):
    """Pre-departure + arrival testing. No surveillance at any point"""
    T = scenario['T']
    CLASSWORK_TRANSMISSION_MULTIPLIER = \
        list(scenario['classwork_transmission_multiplier'].values())
    return \
        Strategy(name="Only Pre-Departure + Arrival Testing",
            # pct_discovered_in_pre_departure=0.5,   # Using Scenario 1
            # pct_discovered_in_arrival_test=0.335,  # Using Scenario 1
            pct_discovered_in_pre_departure=0.25,   # Using Scenario 2
            pct_discovered_in_arrival_test=0.38,  # Using Scenario 2
            testing_regimes=[no_testing_testing_regime(scenario),
                             no_testing_testing_regime(scenario)],
            transmission_multipliers=[1, CLASSWORK_TRANSMISSION_MULTIPLIER],
            period_lengths=[3,T-3-1])


def surge_testing_strategy(scenario: Dict):
    """ Pre-departure + arrival testing. Surveillance of UG and professional
    students before classes and during virtual instruction at 2x/wk. It does
    not surveil GR or FS."""
    T = scenario['T']
    CLASSWORK_TRANSMISSION_MULTIPLIER = \
        list(scenario['classwork_transmission_multiplier'].values())
    return \
        Strategy(name="UG+Prof. 2x/wk in Virtual Instr. Only",
                # pct_discovered_in_pre_departure=0.5,   # Using Scenario 1
                # pct_discovered_in_arrival_test=0.335,  # Using Scenario 1
                pct_discovered_in_pre_departure=0.25,  # Using Scenario 2
                pct_discovered_in_arrival_test=0.38,  # Using Scenario 2
                testing_regimes=[ug_prof_2x_week_testing_regime(scenario),
                                 ug_prof_2x_week_testing_regime(scenario),
                                 no_testing_testing_regime(scenario)],
                transmission_multipliers=[1,
                                          CLASSWORK_TRANSMISSION_MULTIPLIER,
                                          CLASSWORK_TRANSMISSION_MULTIPLIER],
                period_lengths=[3,3,T-6-1])
