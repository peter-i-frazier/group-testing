from strategy import Strategy
from groups import population
from sim import sim
from typing import Dict


class Trajectory:

    def __init__(self, scenario: Dict, strategy: Strategy, sim: sim,
        color: str, name: str=None):
        """Manage all of the objects associated with a trajectory on a graph.

        Args:
            scenario (Dict): Scenario that the simulation was run under.
            strategy (Strategy): Strategy that was used to run the simulation.
            sim (sim): Simulation which used the provided strategy.
            color (str): Color of the trajectory when plotting.
            name (str): Name of the trajectory.
        """
        self.scenario = scenario
        # TODO (hwr26): patch--come back and think about the right thing here
        self.pop = population.from_scenario(scenario)
        self.strategy = strategy
        self.sim = sim
        self.color = color
        self.name = strategy.name if name is None else name
