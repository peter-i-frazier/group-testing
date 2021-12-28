import numpy as np

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
  max_infectious_days: the length of the infectious period if a person is never isolated
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

The max_infectious_days parameter has an important influence on how much testing can help --- if the period is short,
then it is hard for a test to intercept a significant part of this period. We are currently setting this to 6 days as
a default. The logic here is that the pre-infectious-period + max_infectious_period / 2 should be equal to the
generation time. This assumes that the infectivity is uniform and symmetric.  So this would imply
max_infectious_period = 2 * (generation_time - pre_infectious_period). If the generation time is 5 days and the
pre-infectious period is 2 days, then this is 2 * (5-2) = 6 days.

Our simulator is pessimistic in that it assumes PCR cannot detect an infection before it becomes infectious and that
once a person becomes infectious they immediately have the maximum amount of infectivity. In reality, infectivity
builds and a person is likely PCR-detectable when their infectiousness is low.

Here is an article arguing that Omicron may become infectious and generate symptoms more quickly than Delta:
https://www.theatlantic.com/science/archive/2021/12/omicron-incubation-period-testing/621066/
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



def days_infectious(days_between_tests, isolation_delay, sensitivity = .6, max_infectious_days = 6):
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
def __days_infectious_perfect_sensitivity__(days_between_tests, isolation_delay, sensitivity=1, max_infectious_days=6):
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
