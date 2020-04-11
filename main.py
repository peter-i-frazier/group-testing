import numpy as np
from dyadic import sim_dyadic
from two_stage_hierarchical import sim_2stage

def simulate(sim=sim_2stage, reps=100, verbose=True):
    n_tests = [0] * reps
    for m in range(reps):
        n_tests[m] = sim()
        if verbose:
            print('Simulation {}, Result {}'.format(m,n_tests[m]))

    m = np.mean(n_tests)
    stderr = np.std(n_tests) / np.sqrt(reps)
    print("Expected number of tests per person is {:g} +/- {:g}".format(m,stderr))
    print('Point estimate for people per test is {:g}'.format(1/m))
    print("Approx 95% CI for people per test is [{:g},{:g}]".format(1/(m+2*stderr),1/(m-2*stderr)))

print('Simulating 2-stage hierarchical')
simulate(sim_2stage, reps=1000, verbose=False)

print()
print('Simulating Dyadic')
simulate(sim_dyadic, reps=10, verbose=False)
