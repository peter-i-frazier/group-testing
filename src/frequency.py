import numpy as np

def frequency(QFNR, alpha, beta=0.5**(1/7)):
    # INPUTS:
    # QFNR: Quarantine false negative rate = Prob(a person is released by group test | he/she is positive)
    # alpha: 1 infected case grows to alpha^d infected cases after d days (including secondary, etc. infections)
    # beta: safety-factor, the daily decrease rate of the number of infections

    # RETURNS: number of days between two consecutive screenings, i.e., the *period* between screenings

    assert(QFNR >= 0 and QFNR < 1), "invalid QFNR; valid domain is [0,1)"
    assert(alpha > 0), "invalid R0; valid domain is >0"
    assert(beta > 0 and beta <= 1), "invalid beta; valid domain is (0,1]"

    if np.isclose(QFNR,0.):
        #print('WARNING: QFNR is nearly 0')
        # If the QFNR is 0, then the log-run frequency at which we need to test is 0 since we just quarantine all the
        # positive cases and make sure they don't infect anyone, and then the epidemic is over
        return 0

    return -np.log(QFNR) / (np.log(alpha) - np.log(beta))
