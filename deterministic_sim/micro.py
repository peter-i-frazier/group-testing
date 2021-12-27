import numpy as np
import matplotlib.pyplot as plt
import pytest

'''
This code performs "microscopic" calculations that predict how long an individual is infectious and free.

Note:
- This code does not currently support contact tracing. Later it could.
- It would be nice to predict the average generation time, the fraction of infections that will be
discovered, and when they will be discovered.
'''

'''
days_infectious

Output: The average number of days that a person is infectious before being isolated
Inputs:
  days_between_tests: the number of days between tests
  max_infectious_days: the length of the infectious period according to biology
  isolation_delay: the number of days from sampling to isolation (including delay in the lab and in contact tracing)

Pass days_between_tests = np.inf to model not testing.

We assume the following notation:
-   The start of an infection occurs between two surveillance testing days.  Let T be the distance in time
    between these two sampling days.  
-   We suppose that the beginning of the period when a person is PCR-detectable is uniformly distributed between 
    0 and T.  Call this time X so that T-X is the time until the first surveillance test.
-   We assume that the person is infectious for the whole time they are PCR-detectable, and that their infectivity
    is constant over this time period. Let R be the length of this infectious period.
-   Let D be the delay from sampling to acting on the result by isolating the person.
-   Let N be the index of the the first surveillance test to test positive, where N=0 is the first test.

With this notation, the number of days when the person is infectious is min(D+NT+T-X, R).
We compute the expected value of this quantity.
'''


def __conditional_days_infectious__(n, days_between_tests, isolation_delay, max_infectious_days):
    '''
    Computes the conditional expectation of the number of days infectious assuming that surveillance test n will be
    the first one to test positive, where n=0 is the first test.

    Conditioned on surveillance test n being the first to test positive,
    the number of days when the person is infectious is min(D+nT+T-X, R).

    We can rewrite this as D + nT + min(T-X,R-D-nT) = D + nT + T * min((T-X)/T, (R-D-nT)/T)

    U = (T-X)/T is uniformly distributed between 0 and 1.
    Let b = (R-D-NT)/T.  Let's compute y = E[min(U,b)].

    If b <= 0 then y = b
    If b > 1, this is 1/2.
    If b in (0,1), this is b*(1-b/2), according to a pencil & paper calculation

    So the conditional expected time is we'll return D + nT + T * y
    '''

    T = days_between_tests
    D = isolation_delay
    R = max_infectious_days

    assert T > 0

    if T == np.inf:
        return max_infectious_days

    b = (R - D - n*T) / T
    if b < 0:
        y = 0
    elif b > 1:
        y = 0.5
    else:
        y = b * (1 - 0.5 * b)

    return D + n * T + T * y



def days_infectious(days_between_tests, isolation_delay, sensitivity = .6, max_infectious_days = 10):
    '''
    The number of surveillance tests N that are required for a person to test positive is a geometric random variable
    with probability give by the sensitivity parameter. So the person tests positive on test n (where the first test
    is n=0) with probability P(N=n) = sensitivity*np.pow(1-sensitivity,n)

    Conditioned on surveillance test n being the first to test positive (where the first is n=0),
    the number of days when the person is infectious is min(D+nT+T-X, R).  The conditional expectation of this is
    computed by __conditional_days_infectious__
    '''

    T = days_between_tests
    D = isolation_delay
    R = max_infectious_days

    n = 0
    prob = 1 # This contains Prob(N>=n)
    y = 0 # This is the sum of Prob(N=n') * E[days_infectious | N=n'] over 0 <= n' < n
    while D+n*T < R:
        pn = sensitivity * np.power(1-sensitivity,n) # This is Prob(N=n)
        y = y + pn * __conditional_days_infectious__(n,days_between_tests,isolation_delay,max_infectious_days)
        prob = prob - pn
        n = n+1
    # Since X <= T, once D + nT >= R, we have that
    # days_infectious = min(D + nT + T - X, R) = R. This will remain true if n increases.
    # Thus, E[days_infectious | N=n] = R for this n and larger
    # Thus, the sum of Prob(N=n') * E[days_infectious | N=n'] for n' >= n is Prob(N>=n) * R
    y = y + prob*R

    return y





'''
__days_infectious_perfect_sensitivity__

Output: The average number of days that a person is infectious before being isolated
Inputs:
  days_between_tests: the number of days between tests
  max_infectious_days: the length of the infectious period
  isolation_delay: the number of days from sampling to isolation (including delay in the lab and in contact tracing)

Pass days_between_tests = np.inf to model not testing.

Assumes that tests have 100% sensitivity.

Here is the mathematics behind 
Consider the start of an infection.  This occurs between two surveillance testing days.  Let T be the distance in time 
between these two sampling days.  We will suppose that the beginning of the period when a person is PCR-detectable 
is uniformly distributed between 0 and T.  Call this time X.  We assume that the person is infectious for the whole time 
that they are PCR-detectable, and that their infectivity is constant over this time period. Let R be the length of this 
infectious period.  Let D be the delay from sampling to acting on the result by isolating the person. 

Then the number of days when the person is infectious is min(T+D-X, R).

We can rewrite this as D + min(T-X,R-D) = D + T * min((T-X)/T, (R-D)/T)

U = (T-X)/T is uniformly distributed between 0 and 1.
Let b = (R-D)/T.  Let's compute y = E[min(U,b)].

If b <= 0 then y = b
If b > 1, this is 1/2.
If b in (0,1), this is b*(1-b/2), according to a pencil & paper calculation

So we'll return D + T * y
'''
def __days_infectious_perfect_sensitivity__(days_between_tests, isolation_delay, sensitivity=1, max_infectious_days=10):
    assert(sensitivity == 1)

    T = days_between_tests
    D = isolation_delay
    R = max_infectious_days

    assert T > 0

    if T == np.inf:
        return max_infectious_days

    b = (R - D) / T
    if b < 0:
        y = 0
    elif b > 1:
        y = 0.5
    else:
        y = b * (1 - 0.5 * b)

    return D + T * y


def test_days_infectious():

    # If there is infinite time between tests, the time infectious in the presence of testing should be whatever
    # the length of the maximum infectious period is. This should be true regardless of the sensitivity.
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

    days_between_tests = np.arange(1, 20, .1)
    y70 = [days_infectious(d, 1, sensitivity = 0.7) for d in days_between_tests]
    y90 = [days_infectious(d, 1, sensitivity = 0.9) for d in days_between_tests]
    y100 = [days_infectious(d, 1, sensitivity = 1) for d in days_between_tests]
    plt.plot(days_between_tests, y70, label = 'Sensitivity = 70%')
    plt.plot(days_between_tests, y90, label = 'Sensitivity = 90%')
    plt.plot(days_between_tests, y100, label = 'Sensitivity = 100%')
    plt.legend()
    plt.xlabel('Number of days between tests')
    plt.ylabel('Expected Days infectious')
    plt.savefig('test_days_infectious1.png', facecolor='w')
    plt.close()

    # Compare what we are using here to what we used before
    days_between_tests = np.arange(1, 20, .1)
    y_old = [days_infectious(d, 1, sensitivity=1, max_infectious_days=7) for d in days_between_tests]
    y_new = [days_infectious(d, 1) for d in days_between_tests]
    plt.plot(days_between_tests, y_old, 'k--', label='Sensitivity = 100%, max_infectious = 7 days, 1 day delay')
    plt.plot(days_between_tests, y_new, 'k-', label='Sensitivity = 60%, max_infectious = 10 days, 1 day delay')

    y_old = [days_infectious(d, 2, sensitivity=1, max_infectious_days=7) for d in days_between_tests]
    y_new = [days_infectious(d, 2) for d in days_between_tests]
    plt.plot(days_between_tests, y_old, 'b--', label='Sensitivity = 100%, max_infectious = 7 days, 2 day delay')
    plt.plot(days_between_tests, y_new, 'b-', label='Sensitivity = 60%, max_infectious = 10 days, 2 day delay')

    plt.legend()
    plt.xlabel('Number of days between tests')
    plt.ylabel('Expected Days infectious')
    plt.savefig('test_days_infectious2.png', facecolor='w')
    plt.close()

