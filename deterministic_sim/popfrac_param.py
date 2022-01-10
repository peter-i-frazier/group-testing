import numpy as np
from scipy.stats import nbinom, pareto
import matplotlib.pyplot as plt


names = ['UG', 'Grad-Research', 'Professional', 'Employees']

egs = np.array([
        [
            0.3771077769,
            0.1744692313,
            0.2130738259,
            0.07227578963,
            0.05147885598,
            0.04364552029,
            0.02718384893,
            0.01451678525,
            0.02624836593
        ],
        [
            0.5790258783,
            0.1683859957,
            0.05712348893,
            0.04359733917,
            0.03548852695,
            0.0451325654,
            0.03935404056,
            0.01167552257,
            0.01055548401,
            0.009661158416
        ],
        [
            0.472724139,
            0.1443461065,
            0.03264544331,
            0.1245769902,
            0.04056260275,
            0.03439035051,
            0.04498079947,
            0.01334486452,
            0.04825874066,
            0.04416996305
        ],
        [
            0.5559842653,
            0.1914068446,
            0.09034172148,
            0.07325926546,
            0.03157069457,
            0.02081854892,
            0.01296646907,
            0.01615689928,
            0,
            0.005729729349,
            0.00176556195
        ]
    ])

def getBestFit(s, fidelity):
    n = len(egs[s])-1

    ps = np.arange(0,1,10**(-fidelity))
    p_best = -1
    p_best_norm = np.inf
    p_best_hist = np.zeros(n+1)

    for p in ps:
        tmp = np.array([nbinom.pmf(k,n,p) for k in range(n+1)])
        norm = np.linalg.norm(tmp-egs[s])
        if norm<p_best_norm:
            p_best = p
            p_best_norm = norm
            p_best_hist = tmp

    lambs = np.arange(0, 1, 10**(-fidelity))
    lamb_best = -1
    lamb_best_norm = np.inf
    lamb_best_hist = np.zeros(n+1)

    for lamb in lambs:
        tmp = np.array([np.exp(-lamb*k) for k in range(1,n+2)])
        tmp = tmp/np.sum(tmp)
        norm = np.linalg.norm(tmp-egs[s])
        if norm<lamb_best_norm:
            lamb_best = lamb
            lamb_best_norm = norm
            lamb_best_hist = tmp

    bs = np.arange(0.01, 10, 10**(-fidelity))
    b_best = 0
    b_best_norm = np.inf
    b_best_hist = np.zeros(n+1)

    for b in bs:
        tmp = np.array([pareto.pdf(k,b) for k in range(1,n+2)])
        tmp = tmp/np.sum(tmp)
        norm = np.linalg.norm(tmp-egs[s])
        if norm<b_best_norm:
            b_best = b
            b_best_norm = norm
            b_best_hist = tmp


    p_best = np.round(p_best,fidelity)
    lamb_best = np.round(lamb_best,fidelity)
    b_best = np.round(b_best,fidelity)

    plt.plot(np.arange(1, n+2),egs[s], label = 'Current estimate of contact distribution')
    plt.plot(np.arange(1, n+2),p_best_hist, label = 'Best Negative Binomial Fit, p = ' + str(p_best))
    plt.plot(np.arange(1, n+2),lamb_best_hist, label = 'Best Exponential Fit, lambda = ' + str(lamb_best))
    plt.plot(np.arange(1, n+2),b_best_hist, label = 'Best Pareto Fit, b = ' + str(b_best))
    plt.title(names[s] + ' pop_frac Parameterization')
    plt.xlabel('Number of contacts')
    plt.ylabel('Proportion')
    plt.legend()
    # plt.savefig(names[s]+'.png', facecolor = 'w')
    # plt.close()
    plt.show()


