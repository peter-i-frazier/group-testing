import numpy as np

def frequency(QFNR, R0, d0, beta=0.5**(1/7)):
    # INPUTS:
    # QFNR: Quarantine false negative rate = Prob(a person is released by group test | he/she is positive)
    # R0: basic reproduction number
    # d0: duration of infectious period (days)
    # beta: safety-factor, the daily decrease rate of the number of infections

    # RETURNS: number of days between two consecutive screenings

    assert(QFNR > 0 and QFNR < 1), "invalid QFNR; valid domain is (0,1)"
    assert(R0 > 0), "invalid R0; valid domain is >0"
    assert(d0 > 0), "invalid d0; valid domain is >0"
    assert(beta > 0 and beta <= 1), "invalid beta; valid domain is (0,1]"

    return -np.log(QFNR) / (np.log(R0)/d0 - np.log(beta))
