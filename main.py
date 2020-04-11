import numpy as np
from two_stage_hierarchical import sim_2stage
from three_stage_hierarchical import sim_3stage
from dyadic import sim_dyadic

def simulate(sim=sim_2stage, reps=100, verbose=True):
    n_tests = [0] * reps
    for m in range(reps):
        n_tests[m] = sim()
        if verbose:
            print('Simulation {}, Result {}'.format(m,n_tests[m]))

    m = np.mean(n_tests)
    stderr = np.std(n_tests) / np.sqrt(reps)
    print("Using {} Monte Carlo Replications".format(reps))
    print("Expected number of tests per person is {:g} +/- {:g}".format(m,stderr))
    print('Point estimate for people per test is {:g}'.format(1/m))
    print("Approx 95% CI for people per test is [{:g},{:g}]".format(1/(m+1.96*stderr),1/(m-1.96*stderr)))
    print()

print('Simulating 2-stage hierarchical')
simulate(sim_2stage, reps=100000, verbose=False)

print('Simulating 3-stage hierarchical')
simulate(sim_3stage, reps=100000, verbose=False)

print('Simulating Dyadic')
print('Using only 5 replications to run quickly.  Use more to get a smaller CI.')
simulate(sim_dyadic, reps=5, verbose=False)
