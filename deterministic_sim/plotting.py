import numpy as np
from strategy import Strategy
from groups import population
from sim import sim
import matplotlib
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
from typing import List
from functools import reduce
from operator import iconcat, add
from datetime import datetime
from textwrap import fill




# TODO pf98 this really belongs somewhere else, e.g. in the yaml
def long_metagroup_name(x):
    if x == 'PR':
        return 'Professional'
    elif x == 'FS':
        return 'Employees'
    elif x == 'GR':
        return 'Grad-Research'
    else:
        return x


def plot_sm_test_regime_comparison(outfile : str, test_regime_names: List[str],
                                   test_regime_sims: List[sim], test_regime_colors: List[str],
                                   params):

    """Plot a comparison of various test regimes (including no surveillance).
    Args:
        test_regime_names (List[str]): Names of test regimes.
        test_regime_sims (List[sim]): List of test regime simulations.
        test_regime_colors (List[str]): List of colors for test regime trajectories.
        params: Parameters used to run the simulation.
    """

    plt.rcParams["figure.figsize"] = (18,16)
    plt.rcParams['font.size'] = 30
    plt.rcParams['lines.linewidth'] = 6
    plt.rcParams['legend.fontsize'] = 22
    plt.subplots_adjust(hspace = 0.8)
    plt.subplot(211)
    plot_infected_discovered(test_regime_names, test_regime_sims, test_regime_colors, params)
    plt.subplot(212)
    plot_oncampus_isolated(test_regime_names, test_regime_sims, test_regime_colors, params)
    plt.savefig(outfile, facecolor='w')
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
        plt.title("Infections " + reduce(add, [long_metagroup_name(x) for x in metagroup_names]))

    #plt.rcParams.update({'font.size': 8})
    if legend:
        # Shrink current axis by 40%
        ax = plt.gca()
        #box = ax.get_position()
        #ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])

        # Put a legend to the right of the current axis
        # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.20),
                  fancybox=True, shadow=True, ncol=2)
        # ax.legend(loc = 'lower right', prop={'size': 16}, bbox_to_anchor = (1,0.45))
    plt.ylabel('Cumulative Infected')

def plot_oncampus_isolated(strategies: List[Strategy],
                           test_regime_names: List[str],
                           test_regime_sims: List[sim], test_regime_colors: List[str],
                           params, legend = True):
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
        isolated = s.get_isolated(arrival_discovered=sum(strategies[i].get_active_discovered(params)),
                                  iso_lengths=params["isolation_durations"],
                                  iso_props=params["isolation_fracs"])
        on_campus_isolated = params["on_campus_frac"] * isolated
        plt.plot(X, on_campus_isolated, label=label, color=color)

    plt.title("On-campus Isolation")
    #plt.rcParams.update({'font.size': 8})
    if legend:
        plt.legend()
    plt.xlabel('Days')
    plt.ylabel('Isolation (on-campus 5 day)')

def plot_comprehensive_summary(outfile: str,
                                strategies: List[Strategy],
                                   test_regime_names: List[str],
                                   test_regime_sims: List[sim], test_regime_colors: List[str],
                                   params, popul, simple_param_summary = None):
    """Plot a comprehensive summary of the simulation run."""
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(8.5, 11)
    plt.rcParams.update({'font.size': 8})

    plt.subplot(411) # Take up the whole top row
    plot_infected_discovered(test_regime_names, test_regime_sims, test_regime_colors, params, legend = True)
    window = 423 # Start in the second row

    plt.subplot(window)
    plot_oncampus_isolated(strategies, test_regime_names, test_regime_sims, test_regime_colors, params, legend = False)
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


    #plt.rcParams.update({'font.size': 8})
    if simple_param_summary is None:
        plt.text(0,-0.5,param2txt(params))
    else:
        now = datetime.now()
        plt.text(0,0.5,'{}\nSimulation run {}'.format(fill(simple_param_summary, 60),now.strftime('%Y/%m/%d %H:%M')))

    plt.tight_layout(pad=1)

    plt.savefig(outfile, facecolor='w')
    plt.close()

def plot_hospitalization(outfile, test_regime_names: List[str],
                             test_regime_sims: List[sim], test_regime_colors: List[str],
                             params, popul, legend = True):
    """Plot infected and discovered under several test regimes (including no surveillance).

    Args:
        outfile: name of output .png
        test_regime_names (List[str]): Names of test regimes.
        test_regime_sims (List[sim]): List of test regime simulations.
        test_regime_colors (List[str]): List of colors for test regime trajectories.
        params: Parameters used to run the simulation.
        popul: Population object, used only for getting meta-group/group names
        metagroup_names: list of names of meta-group(s) to plot, None to plot the sum across groups.
    """
    plt.rcParams["figure.figsize"] = (8,6)
    plt.rcParams['font.size'] = 15
    plt.rcParams['lines.linewidth'] = 6
    plt.rcParams['legend.fontsize'] = 12
    for i in range(len(test_regime_names)):
        label = test_regime_names[i]
        s = test_regime_sims[i]
        color = test_regime_colors[i]
        X = np.arange(s.max_T) * s.generation_time # Days in the semester, to plot on the x-axis

        hospitalized = np.zeros(s.max_T)
        group_idxs = popul.metagroup_indices(['UG', 'GR', 'PR', 'FS'])
        for i in range(4):
            hospitalized += \
                s.get_total_infected_for_different_groups(group_idxs[i], cumulative=True) * \
                params["hospitalization_rates"][i]

        plt.plot(X, hospitalized, label=label, color=color, linestyle = 'solid')
        plt.title("Spring Semester Hospitalizations, Students+Employees")

    #plt.rcParams.update({'font.size': 8})
    if legend:
        plt.legend()
    plt.ylabel('Cumulative Hospitalized')
    plt.savefig(outfile, facecolor='w')
    plt.close()

def param2txt(params):
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
    return param_txt
