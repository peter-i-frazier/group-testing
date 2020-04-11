import numpy as np

def sim_2stage():
    pool_size = 32
    p = 0.001
    # x is a vector that
    x = np.random.rand(pool_size) < p

    n_tests = 1
    positive = x[1:pool_size].max()
    if positive:
        n_tests += pool_size
    return n_tests


m = int(1e4) # MC reps
y = [0]*m
for i in range(m):
    y[i] = sim_2stage()

exp_n_tests = np.mean(y)
print(exp_n_tests)
print(np.std(y)/np.sqrt(m))


#print(n_tests)

#print(sum(x)/(1.0*n))