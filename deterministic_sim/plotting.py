import numpy as np
from groups import population
from sim import sim
import matplotlib
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
from typing import List
from functools import reduce
from operator import iconcat, add

def plot_sm_test_regime_comparison(test_regime_names: List[str],
                                   test_regime_sims: List[sim], test_regime_colors: List[str],
                                   params):

    """Plot a comparison of various test regimes (including no surveillance).
    Args:
        test_regime_names (List[str]): Names of test regimes.
        test_regime_sims (List[sim]): List of test regime simulations.
        test_regime_colors (List[str]): List of colors for test regime trajectories.
        params: Parameters used to run the simulation.
    """
    plt.subplot(211)
    plot_infected_discovered(test_regime_names, test_regime_sims, test_regime_colors, params)
    plt.subplot(212)
    plot_oncampus_isolated(test_regime_names, test_regime_sims, test_regime_colors, params)
    plt.savefig('sp22_sim_simple.png', facecolor='w')
    plt.close()



def plot_infected_discovered(test_regime_names: List[str],
                             test_regime_sims: List[sim], test_regime_colors: List[str],
                             params, popul = None, metagroup_names : List[str] = None, legend = True):
    """Plot infected and discovered under several test regimes (including no surveillance).

    Args:
        test_regime_names (List[str]): Names of test regimes.
        test_regime_sims (List[sim]): List of test regime simulations.
        test_regime_colors (List[str]): List of colors for test regime trajectories.
        params: Parameters used to run the simulation.
        metagroup_names: list of names of meta-group(s) to plot, None to plot the sum across groups.
    """
    for i in range(len(test_regime_names)):
        label = test_regime_names[i]
        s = test_regime_sims[i]
        color = test_regime_colors[i]

        X = np.arange(s.max_T) * s.generation_time # Days in the semester, to plot on the x-axis
        if metagroup_names == None:
            discovered = s.get_discovered(aggregate=True, cumulative=True)
            infected = s.get_infected(aggregate=True, cumulative=True)
        else:
            group_idx = popul.metagroup_indices(metagroup_names)
            # Since metagroup_names is a list of metagroup names, group_idix will be a list of lists.
            # We want to flatten it. See https://stackabuse.com/python-how-to-flatten-list-of-lists/
            group_idx = reduce(iconcat, group_idx, [])
            discovered = s.get_total_discovered_for_different_groups(group_idx, cumulative=True)
            infected = s.get_total_infected_for_different_groups(group_idx, cumulative=True)
        if np.isclose(discovered,infected).all():
            # Discovered and infected are the same, or almost the same.  This occurs when we do surveillance.
            # Only plot one line
            plt.plot(X, discovered, label=label, color=color, linestyle = 'solid')
        else:
            plt.plot(X, discovered, label=label + '(Discovered)', color=color, linestyle = 'solid')
            plt.plot(X, infected, label=label + '(Infected)', color=color, linestyle = 'dashed')

    if metagroup_names == None:
        plt.title("Spring Semester Infections, Students+Employees")
    else:
        # Plots all of the metagroup names together
        plt.title("Infections " + reduce(add, metagroup_names))
    plt.rcParams.update({'font.size': 8})
    if legend:
        plt.legend()
    plt.ylabel('Cumulative Infected')

def plot_oncampus_isolated(test_regime_names: List[str],
                           test_regime_sims: List[sim], test_regime_colors: List[str],
                           params):
    """Plot the number of rooms of isolation required to isolate on-campus students under
    the passed set of test regimes.

    Args:
        test_regime_names (List[str]): Names of test regimes.
        test_regime_sims (List[sim]): List of test regime simulations.
        test_regime_colors (List[str]): List of colors for test regime trajectories.
        params: Parameters used to run the simulation.
    """
    for i in range(len(test_regime_names)):
        label = test_regime_names[i]
        s = test_regime_sims[i]
        color = test_regime_colors[i]

        X = np.arange(s.max_T) * s.generation_time # Days in the semester, to plot on the x-axis
        isolated = s.get_isolated(iso_lengths=params["isolation_durations"],
                                  iso_props=params["isolation_fracs"],
                                  on_campus_frac=params["on_campus_frac"])
        plt.plot(X, isolated, label=label, color=color)

    plt.title("On-campus Isolation")
    plt.rcParams.update({'font.size': 8})
    plt.legend()
    plt.xlabel('Days')
    plt.ylabel('Isolation (on-campus 5 day)')

def plot_comprehensive_summary(test_regime_names: List[str],
                                   test_regime_sims: List[sim], test_regime_colors: List[str],
                                   params, popul):
    """Plot a comprehensive summary of the simulation run."""
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(8.5, 11)

    window = 421

    plt.subplot(window)
    plot_infected_discovered(test_regime_names, test_regime_sims, test_regime_colors, params, legend = False)
    window += 1

    plt.subplot(window)
    plot_oncampus_isolated(test_regime_names, test_regime_sims, test_regime_colors, params)
    window += 1

    metagroups = popul.metagroup_names()

    # Plot infected and discovered for each meta-group
    if len(metagroups) > 1:
        for i in range(len(metagroups)):
            plt.subplot(window)
            window += 1
            plot_infected_discovered(test_regime_names, test_regime_sims, test_regime_colors, \
                                     params, popul, [metagroups[i]], legend = False)

    plt.subplot(window)
    plt.axis('off')
    window += 1

    param_txt = ''
    for param_name in params:
        param = params[param_name]
        if param_name == 'meta_matrix':
            np.set_printoptions(precision = 2)
            param_txt = param_txt + '\nmeta_matrix:\n' + str(np.matrix(param))
        elif param_name == 'pop_fracs':
            # skip
            1 == 1
        else:
            param_txt = param_txt + '\n' + param_name + ':' + str(param)

    plt.rcParams.update({'font.size': 8})
    plt.text(0,-0.5,param_txt)

    plt.savefig('sp22_sim_comprehensive.png', facecolor='w')
    plt.close()