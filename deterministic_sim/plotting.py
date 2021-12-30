import numpy as np
from sim import sim
import matplotlib
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
from typing import List


def plot_sm_test_regime_comparison(test_regime_names: List[str],
    test_regime_sims: List[sim], test_regime_colors: List[str],
    no_surveillance_sim: sim, params):
    """Plot a comparison of various test regimes (including no surveillance).

    Args:
        test_regime_names (List[str]): Names of test regimes.
        test_regime_sims (List[sim]): List of test regime simulations.
        test_regime_colors (List[str]): List of colors for test regime trajectories.
        no_surveillance_sim (sim): No surveillance simulation.
        params: Parameters used to run the simulation.
    """

    # create test regime trajectories
    for i in range(len(test_regime_names)):
        label = test_regime_names[i]
        s = test_regime_sims[i]
        color = test_regime_colors[i]

        plt.subplot(211)
        plt.plot(np.arange(s.max_T)*s.generation_time, s.get_discovered(aggregate=True,cumulative=True), label=label, color=color)
        plt.subplot(212)
        isolated = s.get_isolated(iso_lengths=params["isolation_durations"],
                                  iso_props=params["isolation_fracs"],
                                  on_campus_frac=params["on_campus_frac"])
        plt.plot(np.arange(s.max_T)*s.generation_time, isolated, label=label, color=color)

    # create no surveillance trajectories
    s = no_surveillance_sim
    plt.subplot(211)
    plt.plot(np.arange(s.max_T)*s.generation_time, s.get_discovered(aggregate=True,cumulative=True), 'k-', label='No surveillance, Discovered')
    plt.plot(np.arange(s.max_T)*s.generation_time, s.get_infected(aggregate=True,cumulative=True), 'k--', label='No surveillance, Infected')
    plt.subplot(212)
    isolated = s.get_isolated(iso_lengths=params["isolation_durations"],
                              iso_props=params["isolation_fracs"],
                              on_campus_frac=params["on_campus_frac"])
    plt.plot(np.arange(s.max_T)*s.generation_time, isolated, 'k', label='No surveillance')

    # create and export plot
    plt.subplot(211)
    plt.title(f'dec_ug_infected_per_day_unit={params["dec_ug_infected_per_day_unit"]}, Symptomatic Rate = {params["symptomatic_rate"]}')
    plt.rcParams.update({'font.size': 8})
    plt.legend()
    plt.ylabel('UG Infected')

    plt.subplot(212)
    plt.rcParams.update({'font.size': 8})
    plt.legend()
    plt.xlabel('Weeks')
    plt.ylabel('UG in Isolation (on-campus 5day)')

    plt.savefig('sp22_sim.png', facecolor='w')
