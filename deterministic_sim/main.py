import micro
import macro_old
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

pop = 16000 #total UG population
K = 12 #number of distinct contact groups
lamb = 0.3 #probability of infection of contact
T = 5 #num generations

# Based on Xiangyu's adjusted moment match code, this is the fraction of the population (UG-only?  or all students)
# broken out by the amount of contact that they have, starting from low contact to high contact.
q = np.array(
[0.46666204134859246,
 0.21110377424326393,
 0.13427835918082992,
 0.040993148549700376,
 0.05561449714374835,
 0.03772148571886847,
 0.020557396664413034,
 0.01829685176380435,
 0.003308325172067398,
 0.006056046991853723,
 0.0027991704708900224,
 0.002608902751968122])

# average the omicron outbreak numbers over 3 generations
init = np.array([140, 143, 106, 82, 54, 49, 34, 24, 21, 11, 8, 10])/(13/5.7)

if __name__ == '__main__':
    # sim(K, lamb, T, I0, R0)
    infected = pop * macro.sim(K, 1, 10, )


def old_main():
    T = 20
    res = np.zeros((7, T))
    lambs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    for i in range(7):
        # macro.sim(pop, K, lamb, T, q, x0)
        res[i] = pop * macro.old_sim(pop, K, lambs[i], T, q, init)

    for i in range(7):
        plt.plot(np.arange(T), res[i], label=lambs[i])
    plt.legend(title='Prob. Infect. per Contact')
    plt.xlabel('Generation')
    plt.ylabel('Cumulative Number of Infected')
    plt.title('Infections under Contact Heterogeneity 2')
    plt.savefig('conctact_get.png', facecolor='w')