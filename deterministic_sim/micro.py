import yaml
import numpy as np

'''
This module performs "microscopic" calculations that predict how many people an infectious individual will infect in 
a fully susceptible population. This is provided by the secondary_infections function.

Note:
- Later it would be nice to put this into a class
- This code does not currently support contact tracing. Later it could, by supposing that 

It also predicts the average generation time, the fraction of infections that will be
discovered, and when they will be discovered.

.  We simulate forward in a detailed way an individual being infected, then being
discovered in a test

the path forward from the moment that an individual
becomes infecte
'''

def secondary_infections(infection_rate, days_between_tests, infectious_days, isolation_delay):
    '''
    Predicts how many people would be infected by a single infectious source in a fully susceptible population.
    :param infection_rate: The number of people exposed by a source per day of being infectious & free
    :param days_between_tests:
    :param infectious_days:
    :param isolation_delay:
    :return:
    '''
    return days_infectious(days_between_tests,infectious_days,isolation_delay)*infection_rate



'''
days_infectious

Output: The average number of days that a person is infectious before being isolated
Inputs:
  days_between_tests: the number of days between tests
  infectious_days: the length of the infectious period
  isolation_delay: the number of days from sampling to isolation (including delay in the lab and in contact tracing)

Pass days_between_tests = np.inf to model not testing.

Here is the mathematics behind 
Consider the start of an infection.  This occurs between two surveillance testing days.  Let T be the distance in time 
between these two sampling days.  We will suppose that the beginning of the period when a person is PCR-detectable 
is uniformly distributed between 0 and T.  Call this time X.  We assume that the person is infectious for the whole time 
that they are PCR-detectable, and that their infectivity is constant over this time period. Let R be the length of this 
infectious period.  Let D be the delay from sampling to acting on the result by isolating the person. We assume that 
tests have 100% sensitivity.

Then the number of days when the person is infectious is min(T+D-X, R).

We can rewrite this as D + min(T-X,R-D) = D + T * min((T-X)/T, (R-D)/T)

U = (T-X)/T is uniformly distributed between 0 and 1.
Let b = (R-D)/T.  Let's compute y = E[min(U,b)].

If b <= 0 then y = b
If b > 1, this is 1/2.
If b in (0,1), this is b*(1-b/2), according to a pencil & paper calculation

So we'll return D + T * y
'''
def days_infectious(days_between_tests, infectious_days, isolation_delay):
    T = days_between_tests
    D = isolation_delay
    R = infectious_days

    assert T > 0

    if T == np.inf:
        return infectious_days

    b = (R - D) / T
    if b < 0:
        y = 0
    elif b > 1:
        y = 0.5
    else:
        y = b * (1 - 0.5 * b)

    return D + T * y


assert days_infectious(np.inf,5,1) == 5