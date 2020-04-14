def gollier(n=10000, starting_prevalence=0.001, alpha=np.exp(0.13), p_decline_rate=0.8, max_pool_size=1000, FNR=0,
            household_SAR=0.5, T=30, screening_freq = 7):
    # INPUTS:
    # n = number of individuals
    # prevalence
    # alpha = daily growth rate of the infected population size
    # p_decline = desired rate at which you want prevalence to decrease
    # max_pool_size
    # FNR = false negative rate at max pool size
    # household_SAR = household secondary attack rate = probability that an individual in an exposed household is infected (need to confirm)
    # T = length of simulation, in days

    # OUTPUTS:
    # desired frequency of screening to have prevalence decline at desired p_decline_rate
    # fraction of infected population over time
    # fraction of deconfined population over time
    # resource use: n_tests,

    prevalence = np.zeros(T + 1)
    prevalence[0] = starting_prevalence
    deconfined = np.zeros(T + 1)

    screening_freq = 7  # do a screening every 7 days; a tunable parameter

    infection_status = np.random.rand(n) < prevalence[0] # 0 if uninfected; 1 if infected
    quarantine_status = np.zeros(n)  # 0 if quarantined; 1 if released

    # NEXT TO DO: instead of individuals, keep a vector of households
    # [4/14 YZ] This will be done in the population class

    for t in range(T):


        # every day, each uninfected individual has chance ln(alpha)*(current prevalence) of becoming infected
        # [4/14 YZ] This is not exact. Only the unquarantined individuals could become infected
        #           The infection dynamics in households need to be determined
        for i in range(n):
            if infection_status[i] == 0:
                if np.random.binomial(1, np.log(alpha) * prevalence[t]) == 1:
                    infection_status[i] = 1

        # NEXT TO DO: consider intra-/inter-household correlations

        # screenings are carried out every (screening_freq) days
        if np.mod(t, screening_freq) == 0:

            # pool size as recommended in Gollier 2020
            pool_size = np.min([-1 / np.log(1 - prevalence[t]), max_pool_size])
            n_pools = int(np.ceil(n / pool_size))
            pool_counts = np.zeros(n_pools)
            pool_assignment = np.zeros([n_pools, n])

            # assign individuals to pools
            for i in range(n):
                # print(i)
                eligible = np.argwhere(pool_counts < max_pool_size)
                eligible = eligible[:, 0]
                g = np.random.choice(eligible, 1)
                pool_assignment[g, i] = 1
                pool_counts[g] += 1

            # The number of infected individuals in each pool is a g-dim vector given by pool_assignment * infection_status
            pool_infected_counts = np.matmul(pool_assignment, infection_status)

            # Pools without any infected individuals test negative (assuming 100% test accuracy)
            neg_pools = np.sum(pool_infected_counts == np.zeros(n_pools))

            # NEXT TO DO: consider false-negatives

            # Release people in negative groups from quarantine
            for g in range(len(pool_infected_counts)):
                if pool_infected_counts[g] == 0:
                    for i in np.nonzero(pool_assignment[g, :]):
                        quarantine_status[i] = 1

        # update prevalence and number of deconfined individuals
        # (4/14 YZ) This is not exact. Prevalence calculation is based on unquarantined people only, i.e.
        # prevalence at time t = (number of unquarantined infected individuals) / (number of unquarantined individuals)
        prevalence[t + 1] = np.sum(infection_status) / n
        print('prevalence on day ' + str(t + 1) + ': ' + str(prevalence[t + 1] / n))
        deconfined[t + 1] = np.sum(quarantine_status)
        print('deconfined on day ' + str(t + 1) + ': ' + str(deconfined[t + 1]))

    return prevalence, deconfined, infection_status, quarantine_status
