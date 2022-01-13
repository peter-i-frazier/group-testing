import numpy as np
import pandas as pd
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


class Trajectory:

    def __init__(self, strategy: Strategy, sim: sim, color: str):
        """Manage all of the objects associated with a trajectory on a graph.

        Args:
            strategy (Strategy): Strategy that was used to run the simulation.
            sim (sim): Simulation which used the provided strategy.
            color (str): Color of the trajectory when plotting.
        """
        self.strategy = strategy
        self.sim = sim
        self.color = color


def plot_small_summary(outfile : str,
                       trajectories: List[Trajectory],
                       params):
    """Plot a small summary of the simulation run."""
    plt.rcParams["figure.figsize"] = (18,16)
    plt.rcParams['font.size'] = 30
    plt.rcParams['lines.linewidth'] = 6
    plt.rcParams['legend.fontsize'] = 22
    plt.subplots_adjust(hspace = 0.8)
    plt.subplot(211)
    plot_infected_discovered(trajectories, params)
    plt.subplot(212)
    plot_isolated(trajectories, params, oncampus=True)
    plt.savefig(outfile, facecolor='w')
    plt.close()


def plot_infected_discovered(trajectories: List[Trajectory],
                             params,
                             popul = None,
                             metagroup_names : List[str] = None,
                             legend = True):
    """Plot infected and discovered for several trajectories.

    Args:
        params: Parameters used to run the simulation.
        metagroup_names: list of names of meta-group(s) to plot, \
            None to plot the sum across groups.
    """
    # plot each trajectory
    for trajectory in trajectories:
        label = trajectory.strategy.name
        s = trajectory.sim
        color = trajectory.color

        X = np.arange(s.max_T) * s.generation_time  # days in the semester
        if metagroup_names == None:
            discovered = s.get_discovered(aggregate=True, cumulative=True)
            infected = s.get_infected(aggregate=True, cumulative=True)
        else:
            group_idx = popul.metagroup_indices(metagroup_names)
            group_idx = reduce(iconcat, group_idx, [])  # flatten
            discovered = s.get_total_discovered_for_different_groups(group_idx, cumulative=True)
            infected = s.get_total_infected_for_different_groups(group_idx, cumulative=True)
        if np.isclose(discovered, infected).all():
            # Discovered and infected are the same, or almost the same.
            # This occurs when we do surveillance.
            # Only plot one line.
            plt.plot(X, discovered, label=label, color=color, linestyle = 'solid')
        else:
            plt.plot(X, discovered, label=label + '(Discovered)', color=color, linestyle = 'solid')
            # plt.plot(X, infected, label=label + '(Infected)', color=color, linestyle = 'dashed')
            plt.axvline(7,linestyle='--',color='grey')
            plt.axvline(21,linestyle='--',color='grey')


    if metagroup_names == None:
        plt.title("Spring Semester Infections, Students+Employees")
    else:
        plt.title("Infections " + reduce(add, [params["population_names"][x] for x in metagroup_names]))

    if legend:
        ax = plt.gca()
        # Put legend below the current axis because it's too big
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.20),
                  fancybox=True, shadow=True, ncol=2)
    plt.ylabel('Cumulative Infected')


def _get_isolated(s : sim, params,
                  popul = None,
                  metagroup_names = None,
                  metagroup_idx = None,
                  active_discovered = None):
    """
    Helper function that gets the number isolated from a simulation object (which does not include the additional
    arrival positives) for some or all metagroups.  If you are getting all metagroups, popul, metagroup_names,
    and metagroup_idx should be None.
    active_discovered is either None or a numpy array (one for each metagroup)
    """

    if metagroup_names is None:
        if active_discovered is not None:
            active_discovered = sum(active_discovered) # Need to aggregate over all the metagroups
        isolated = s.get_isolated(arrival_discovered=active_discovered,
                                  iso_lengths=params["isolation_durations"],
                                  iso_props=params["isolation_fracs"])
        return isolated

    # Here, metagroup_names is not none
    # TODO: assert that metagroup_idx has the right length, population is not none

    # Get the list of group indices
    group_idx = popul.metagroup_indices(metagroup_names)  # these indices are at the group level, but in a list of lists
    idx = reduce(iconcat, group_idx, [])  # flatten to just a list

    # Get the active discovered for the desired metagroup indices
    if active_discovered is not None:
        metagroup_active_discovered = sum(active_discovered[metagroup_idx])
    else:
        metagroup_active_discovered = None

    isolated = s.get_isolated(group=idx, arrival_discovered=metagroup_active_discovered,
                              iso_lengths=params["isolation_durations"],
                              iso_props=params["isolation_fracs"])

    return isolated

# TODO pf98: Instead of defaulting metagroup_names and metagroup_idx to the explicit list of metagroups,
# set them to None and have the code do the aggregation
def plot_isolated(trajectories: List[Trajectory],
                            params,
                            legend = True,
                            popul = None, metagroup_names = None, metagroup_idx = None,
                            oncampus = False):
    """Plot the number of rooms of isolation required to isolate on-campus
    students under the passed set of test regimes.
    popul, metagroups_names and metagroup_idx are only needed if we getting specific metagroups.
    If so,  metagroups_names and metagroup_idx indicate the metagroups that we wish to include
    Turn on the oncampus flag to apply the oncampus_frac
    """
    for trajectory in trajectories:
        label = trajectory.strategy.name
        s = trajectory.sim
        color = trajectory.color

        X = np.arange(s.max_T) * s.generation_time  # days in the semester

        isolated = _get_isolated(s,params,popul=popul,
                                 metagroup_names=metagroup_names,
                                 metagroup_idx=metagroup_idx,
                                 active_discovered=trajectory.strategy.get_active_discovered(params))

        if oncampus:
            on_campus_isolated = params["on_campus_frac"] * isolated
            plt.plot(X, on_campus_isolated, label=label, color=color)
            if metagroup_names is None:
                plt.title("On-campus Isolation (Students+Employees)")
            else:
                plt.title("On-campus Isolation (" + str(metagroup_names) + ")")
        else:
            plt.plot(X, isolated, label=label, color=color)
            if metagroup_names is None:
                plt.title("Isolation (Students+Employees)")
            else:
                plt.title("Isolation (" + str(metagroup_names) + ")")

    plt.axvline(7,linestyle='--',color='grey')
    plt.axvline(21,linestyle='--',color='grey')

    if legend:
        plt.legend()
    plt.xlabel('Days')
    plt.ylabel('Isolation (5 day)')


def plot_comprehensive_summary(outfile: str,
                               trajectories: List[Trajectory],
                               params, popul, simple_param_summary = None):
    """Plot a comprehensive summary of the simulation run."""
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(8.5, 11)
    plt.rcParams.update({'font.size': 8})

    plt.subplot(421)
    plot_infected_discovered(trajectories, params, legend = True)
    window = 422

    plt.subplot(window)
    plot_isolated(trajectories, params, popul=popul, metagroup_names = ['UG'], metagroup_idx=[0],
                  legend = False, oncampus = True)
    window += 1

    plt.subplot(window)
    plot_isolated(trajectories, params, popul=popul, metagroup_names = ['UG', 'PR'], metagroup_idx=[0],
                  legend = False, oncampus = False)
    window += 1

    metagroups = popul.metagroup_names()

    # Plot infected and discovered for each meta-group
    if len(metagroups) > 1:
        for i in range(len(metagroups)):
            plt.subplot(window)
            window += 1
            plot_infected_discovered(trajectories, params, popul, [metagroups[i]], legend = False)


    def print_params():
        if simple_param_summary is None:
            plt.text(0,-0.5,param2txt(params))
        else:
            now = datetime.now()
            plt.text(0,0.5,'{}\nSimulation run {}'.format(fill(simple_param_summary, 60),now.strftime('%Y/%m/%d %H:%M')))

    #plt.subplot(window)
    #plt.axis('off')
    #window += 1
    # print_params()

    plt.tight_layout(pad=1)

    plt.savefig(outfile, facecolor='w')
    plt.close()


def plot_hospitalization(outfile,
                         trajectories: List[Trajectory],
                         params, popul, legend = True):
    """Plot total hospitalizations for multiple trajectories."""
    plt.rcParams["figure.figsize"] = (8,6)
    plt.rcParams['font.size'] = 15
    plt.rcParams['lines.linewidth'] = 6
    plt.rcParams['legend.fontsize'] = 12
    for trajectory in trajectories:
        label = trajectory.strategy.name
        s = trajectory.sim
        color = trajectory.color
        X = np.arange(s.max_T) * s.generation_time  # days in the semester

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
        elif param_name == 'population_names':
            param_txt +=  f"\n{param_name}: {str(list(param.keys()))}"
        else:
            param_txt = param_txt + '\n' + param_name + ':' + str(param)
    return param_txt


def summary_statistics(outfile: str,
                       trajectories: List[Trajectory],
                       params, popul):
    """Output a CSV file with summary statistics"""

    df = {}

    def get_peak_hotel_rooms(trajectory):
        s = trajectory.sim
        strat = trajectory.strategy
        isolated = _get_isolated(s, params,
                                 popul = popul, metagroup_names = ['UG'], metagroup_idx = [0],
                                 active_discovered = strat.get_active_discovered(params))  # Get UG isolation only
        on_campus_isolated = params["on_campus_frac"] * isolated
        return int(np.ceil(np.max(on_campus_isolated)))

    df["Hotel Room Peaks"] = \
        {t.strategy.name : get_peak_hotel_rooms(t) for t in trajectories}

    def get_total_hotel_rooms(trajectory):
        s = trajectory.sim
        strat = trajectory.strategy
        isolated = _get_isolated(s, params,
                                 popul = popul, metagroup_names = ['UG'], metagroup_idx = [0],
                                 active_discovered = strat.get_active_discovered(params))  # Get UG isolation only
        on_campus_isolated = params["on_campus_frac"] * isolated
        print(on_campus_isolated)
        return int(np.ceil(np.sum(on_campus_isolated) * params["generation_time"]))

    df["Total Hotel Rooms"] = \
        {t.strategy.name : get_total_hotel_rooms(t) for t in trajectories}


    def get_ug_prof_days_in_isolation_in_person(trajectory):
        s = trajectory.sim
        strat = trajectory.strategy
        isolated = _get_isolated(s, params, popul = popul, metagroup_names = ['UG', 'PR'], metagroup_idx = [0,2],
                                 active_discovered = strat.get_active_discovered(params))
        START_OF_IN_PERSON = 5 # generation when we start in-person instruction
        return int(np.sum(isolated[START_OF_IN_PERSON:])*params["generation_time"])

    def get_ug_prof_days_in_isolation(trajectory):
        s = trajectory.sim
        strat = trajectory.strategy
        isolated = _get_isolated(s, params, popul = popul, metagroup_names = ['UG', 'PR'], metagroup_idx = [0,2],
                                 active_discovered = strat.get_active_discovered(params))
        return int(np.sum(isolated)*params["generation_time"])

    df["UG+PR Days In Isolation In Person"] = \
        {t.strategy.name: get_ug_prof_days_in_isolation_in_person(t) for t in trajectories}

    df["UG+PR Days In Isolation (All Time)"] = \
        {t.strategy.name: get_ug_prof_days_in_isolation(t) for t in trajectories}

    def get_total_hospitalizations(trajectory):
        s = trajectory.sim
        hospitalized = np.zeros(s.max_T)
        group_idxs = popul.metagroup_indices(['UG', 'GR', 'PR', 'FS'])
        for i in range(4):
            hospitalized += \
                s.get_total_infected_for_different_groups(group_idxs[i], cumulative=True) * \
                params["hospitalization_rates"][i]
        return hospitalized[-1]

    df["Hospitalizations"] = \
        {t.strategy.name : get_total_hospitalizations(t) for t in trajectories}

    def get_cumulative_infections(trajectory):
        infected = trajectory.sim.get_infected(aggregate=True, cumulative=True)
        return int(np.ceil(infected[-1]))

    df["Cumulative Infections"] = \
        {t.strategy.name : get_cumulative_infections(t) for t in trajectories}

    pd.DataFrame(df).T.to_csv(outfile)
