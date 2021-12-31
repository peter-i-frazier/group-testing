import numpy as np
from groups import population
from sim import sim
import matplotlib
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
from typing import List

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
    plt.savefig('sp22_sim.png', facecolor='w')



def plot_infected_discovered(test_regime_names: List[str],
                             test_regime_sims: List[sim], test_regime_colors: List[str],
                             params):
    """Plot infected and discovered under several test regimes (including no surveillance).

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
        discovered = s.get_discovered(aggregate=True, cumulative=True)
        infected = s.get_infected(aggregate=True, cumulative=True)
        if np.isclose(discovered,infected).all():
            # Discovered and infected are the same, or almost the same.  This occurs when we do surveillance.
            # Only plot one line
            plt.plot(X, discovered, label=label, color=color, linestyle = 'solid')
        else:
            plt.plot(X, discovered, label=label + '(Discovered)', color=color, linestyle = 'solid')
            plt.plot(X, infected, label=label + '(Infected)', color=color, linestyle = 'dashed')

    plt.title("Projections for Spring 2022 Semester")
    plt.rcParams.update({'font.size': 8})
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

    plt.rcParams.update({'font.size': 8})
    plt.legend()
    plt.xlabel('Days')
    plt.ylabel('Isolation (on-campus 5 day)')

def plot_comprehensive_summary(s: sim, pop: population, params):
    """Plot a comprehensive summary of the simulation run.

    Args:
        s (sim): Simulation.
        pop (population): Population used in the simulation.
        params ([type]): Parameters used in the simulation.
    """
    plt.figure(figsize=(18, 12), dpi=300)
    X = np.arange(s.max_T)*s.generation_time

    # plot infectious / discovered for each meta-group
    groups = pop.metagroup_indices(params["population_names"])
    for i in range(4):
        infectious = s.get_total_infected_for_different_groups(groups[i], cumulative=True)
        discovered = s.get_total_discovered_for_different_groups(groups[i], cumulative=True)
        infectious_lbl = "Infected: " + params["population_names"][i]
        discovered_lbl = "Discovered: " + params["population_names"][i]
        plt.subplot(int("24" + str(i + 1)))
        plt.plot(X, infectious, 'k--', label=infectious_lbl, color='r')
        plt.plot(X, discovered, 'k--', label=discovered_lbl, color='blue')
        plt.title("No surveillance")
        plt.legend()

    # TODO: Weighted sum of infections across meta-groups that is intended to
    # be hospitalizations, where the employee infections result in more
    # hospitalizations per infection because we are older
    plt.subplot(245)


    # TODO: On-campus UG + grad-research + grad-professional in isolation
    # (since this determines our housing needs)
    # NOTE: this is placeholder code
    plt.subplot(246)
    isolated = s.get_isolated(iso_lengths=params["isolation_durations"],
                              iso_props=params["isolation_fracs"],
                              on_campus_frac=params["on_campus_frac"])
    plt.plot(X, isolated, 'k', label='No surveillance')

    # TODO: All UG + grad-professional in isolation
    # NOTE: this is placeholder code
    plt.subplot(247)
    isolated = s.get_isolated(iso_lengths=params["isolation_durations"],
                              iso_props=params["isolation_fracs"],
                              on_campus_frac=params["on_campus_frac"])
    plt.plot(X, isolated, 'k', label='No surveillance')

    # TODO: Some text that includes all of the parameters that were used,
    # the github commit, and a timestamp
    plt.subplot(248)
    # plt.text()

    plt.savefig('sp22_sim_summary.png', facecolor='w')
