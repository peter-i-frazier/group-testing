from functools import reduce
from operator import iconcat
import numpy as np
from trajectory import Trajectory
from typing import List


def get_isolated(trajectory: Trajectory, metagroup_names: List[str] = None):
    """Return the number of isolated individuals at each generation.

    Args:
        trajectory (Trajectory): Trajectory object.
        metagroup_names (List[str]): Limit to isolated in these metagroups.

    Returns:
        np.ndarray: Number of isolated individuals at each generation.
    """
    sim = trajectory.sim
    scenario = trajectory.scenario

    if metagroup_names is None:
        active_discovered = \
            sum(trajectory.strategy.get_active_discovered(scenario))
        isolated = sim.get_isolated(arrival_discovered=active_discovered,
                                    iso_lengths=scenario["isolation_durations"],
                                    iso_props=scenario["isolation_fracs"])
        return isolated
    else:
        group_idx = trajectory.pop.metagroup_indices(metagroup_names)
        idx = reduce(iconcat, group_idx, [])
        all_metagroup_names = trajectory.pop.metagroup_names()
        metagroup_idx = [all_metagroup_names.index(i) for i in metagroup_names]
        active_discovered = \
            sum(trajectory.strategy.get_active_discovered(scenario)[metagroup_idx])
        isolated = sim.get_isolated(group=idx, arrival_discovered=active_discovered,
                                iso_lengths=scenario["isolation_durations"],
                                iso_props=scenario["isolation_fracs"])

    return isolated


def get_peak_hotel_rooms(trajectory: Trajectory):
    """Return the peak number of hotel room used over the semester."""
    isolated = get_isolated(trajectory=trajectory,
                            metagroup_names = ['UG'])
    on_campus_isolated = trajectory.scenario["on_campus_frac"] * isolated
    return int(np.ceil(np.max(on_campus_isolated)))


def get_total_hotel_rooms(trajectory: Trajectory):
    """Return the total number of hotel rooms used over the semester."""
    isolated = get_isolated(trajectory=trajectory,
                            metagroup_names = ['UG'])
    scenario = trajectory.scenario
    on_campus_isolated = scenario["on_campus_frac"] * isolated
    return int(np.ceil(np.sum(on_campus_isolated) * scenario["generation_time"]))


def get_ug_prof_days_in_isolation_in_person(trajectory: Trajectory):
    """Return the number of undergrad and professional days in isolation
    during the portion of the semester that is in-person."""
    isolated = get_isolated(trajectory=trajectory,
                            metagroup_names = ['UG', 'PR'])
    # TODO (hwr26): Is there some way to compute this from trajectory.strategy?
    START_OF_IN_PERSON = 5 # generation when we start in-person instruction
    scenario = trajectory.scenario
    return int(np.sum(isolated[START_OF_IN_PERSON:]) * scenario["generation_time"])


def get_ug_prof_days_in_isolation(trajectory: Trajectory):
    """Return the number of undergrad and professional days in isolation."""
    isolated = get_isolated(trajectory=trajectory,
                            metagroup_names = ['UG', 'PR'])
    return int(np.sum(isolated) * trajectory.scenario["generation_time"])


def _get_cumulative_hospitalizations(trajectory: Trajectory, metagroups: List[str]):
    """Return the cumulative number of hospitalizations at each generation."""
    s = trajectory.sim
    pop = trajectory.pop
    cumulative_hospitalizations = np.zeros(s.max_T)
    group_idx = pop.metagroup_indices(metagroups)
    for i in range(len(metagroups)):
        cumulative_hospitalizations += \
            s.get_total_infected_for_different_groups(group_idx[i], cumulative=True) * \
            trajectory.scenario["hospitalization_rates"][metagroups[i]]
    return cumulative_hospitalizations


def get_cumulative_all_hospitalizations(trajectory: Trajectory):
    """Return the cumulative number of all hospitalizations at each generation."""
    pop = trajectory.pop
    metagroups = pop.metagroup_names()
    cumulative_hospitalizations = \
        _get_cumulative_hospitalizations(trajectory, metagroups)
    return cumulative_hospitalizations


def get_total_hospitalizations(trajectory: Trajectory):
    """Return the total number of hospitalizations over the semester."""
    pop = trajectory.pop
    metagroups = pop.metagroup_names()
    cumulative_hospitalizations = \
        _get_cumulative_hospitalizations(trajectory, metagroups)
    return cumulative_hospitalizations[-1]


def get_total_employee_hospitalizations(trajectory: Trajectory):
    """Return the total number of employee hospitalizations over the semester."""
    cumulative_hospitalizations = \
        _get_cumulative_hospitalizations(trajectory, ["FS"])
    return cumulative_hospitalizations[-1]


def get_cumulative_infections(trajectory: Trajectory):
    infected = trajectory.sim.get_infected(aggregate=True, cumulative=True)
    return int(np.ceil(infected[-1]))
