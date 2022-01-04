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
    plot_on_campus_isolated(trajectories, params)
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
            plt.plot(X, infected, label=label + '(Infected)', color=color, linestyle = 'dashed')

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


def plot_on_campus_isolated(trajectories: List[Trajectory],
                            params,
                            legend = True):
    """Plot the number of rooms of isolation required to isolate on-campus
    students under the passed set of test regimes."""
    for trajectory in trajectories:
        label = trajectory.strategy.name
        s = trajectory.sim
        color = trajectory.color

        X = np.arange(s.max_T) * s.generation_time  # days in the semester
        isolated = s.get_isolated(arrival_discovered=sum(trajectory.strategy.get_active_discovered(params)),
                                  iso_lengths=params["isolation_durations"],
                                  iso_props=params["isolation_fracs"])
        on_campus_isolated = params["on_campus_frac"] * isolated
        plt.plot(X, on_campus_isolated, label=label, color=color)

    plt.title("On-campus Isolation")
    if legend:
        plt.legend()
    plt.xlabel('Days')
    plt.ylabel('Isolation (on-campus 5 day)')


def plot_comprehensive_summary(outfile: str,
                               trajectories: List[Trajectory],
                               params, popul, simple_param_summary = None):
    """Plot a comprehensive summary of the simulation run."""
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(8.5, 11)
    plt.rcParams.update({'font.size': 8})

    plt.subplot(411) # Take up the whole top row
    plot_infected_discovered(trajectories, params, legend = True)
    window = 423 # Start in the second row

    plt.subplot(window)
    plot_on_campus_isolated(trajectories, params, legend = False)
    window += 1

    metagroups = popul.metagroup_names()

    # Plot infected and discovered for each meta-group
    if len(metagroups) > 1:
        for i in range(len(metagroups)):
            plt.subplot(window)
            window += 1
            plot_infected_discovered(trajectories, params, popul, [metagroups[i]], legend = False)

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
        isolated = s.get_isolated(arrival_discovered=sum(strat.get_active_discovered(params)),
                                  iso_lengths=params["isolation_durations"],
                                  iso_props=params["isolation_fracs"])
        on_campus_isolated = params["on_campus_frac"] * isolated
        return int(np.ceil(np.max(on_campus_isolated)))

    df["Hotel Room Peaks"] = \
        {t.strategy.name : get_peak_hotel_rooms(t) for t in trajectories}

    def get_total_hotel_rooms(trajectory):
        s = trajectory.sim
        strat = trajectory.strategy
        isolated = s.get_isolated(arrival_discovered=sum(strat.get_active_discovered(params)),
                                  iso_lengths=params["isolation_durations"],
                                  iso_props=params["isolation_fracs"])
        on_campus_isolated = params["on_campus_frac"] * isolated
        return int(np.ceil(np.sum(on_campus_isolated)))

    df["Total Hotel Rooms"] = \
        {t.strategy.name : get_total_hotel_rooms(t) for t in trajectories}

    def get_ug_prof_days_in_isolation_in_person(trajectory):
        # TODO (hwr26): Waiting on unresolved TODO in get_isolated
        raise Exception("Unfinished")

    def get_ug_prof_days_in_isolation(trajectory):
        # TODO (hwr26): Waiting on unresolved TODO in get_isolated
        raise Exception("Unfinished")

    def get_total_hospitalizations(trajectory):
        s = trajectory.sim
        hospitalized = np.zeros(s.max_T)
        group_idxs = popul.metagroup_indices(['UG', 'GR', 'PR', 'FS'])
        for i in range(4):
            hospitalized += \
                s.get_total_infected_for_different_groups(group_idxs[i], cumulative=True) * \
                params["hospitalization_rates"][i]
        return int(np.ceil(hospitalized[-1]))

    df["Hospitalizations"] = \
        {t.strategy.name : get_total_hospitalizations(t) for t in trajectories}

    def get_cumulative_infections(trajectory):
        infected = trajectory.sim.get_infected(aggregate=True, cumulative=True)
        return int(np.ceil(infected[-1]))

    df["Cumulative Infections"] = \
        {t.strategy.name : get_cumulative_infections(t) for t in trajectories}

    pd.DataFrame(df).T.to_csv(outfile)
