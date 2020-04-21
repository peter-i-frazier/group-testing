from scipy.optimize import fsolve
import numpy as np


def eval_r(r, target_prevalence, household_dist, SAR):
    # computes probability of a primary case given population level prevalence, household size distribution,
    # and household secondary attack rate

    # INPUT:
    # r = probability of a primary case in the household
    # target_prevalence = population level prevalence
    # household_dist = probability distribution of household sizes 1,2,3,...
    # SAR = household secondary attack rate

    assert(np.absolute(np.sum(household_dist) -1 ) < 1e-6)

    exp_household_size = 0
    for i in range(len(household_dist)):
        exp_household_size += (i + 1) * household_dist[i]

    frac_tot_infected = 0
    for i in range(len(household_dist)):
        frac_tot_infected += (i + 1) * (r + SAR * (1 - r) - SAR * (1 - r) ** (i + 1)) * household_dist[
            i] / exp_household_size

    return frac_tot_infected - target_prevalence


US_dist = [0.2837, 0.3451, 0.1507, 0.1276, 0.0578, 0.0226, 0.0125]

print("find r for target prevalence = 0.005")
print(fsolve(eval_r,0.005,args=(0.005, US_dist, 0.3741)))

print("find r for target prevalence = 0.01")
print(fsolve(eval_r,0.005,args=(0.01, US_dist, 0.3741)))