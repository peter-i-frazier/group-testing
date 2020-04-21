import numpy as np

def frequency(QFNR, R0, d0, beta=0.5**(1/7)):
    # INPUTS:
    # QFNR: Quarantine false negative rate = Prob(a person is released | he/she is positive)
    # R0: basic reproduction number
    # d0: duration of infectious period (days)
    # beta: safety-factor, the daily decrease rate of the number of infections

    # RETURNS: number of days between two consecutive screenings

    return -np.log(QFNR) / (np.log(R0)/d0 - np.log(beta))
