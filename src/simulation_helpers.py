import matplotlib.pyplot as plt
from population import Population
from group_testing import HouseholdGroupTest
from simulation import Simulation

def initiate_simulation():
    group_test = HouseholdGroupTest(group_test_size,
                                     group_test_participation_rate,
                                     false_negative_rate)

    population = Population(n_households,
                            household_size,
                            initial_prevalence,
                            disease_length,
                            non_quarantine_alpha,
                            daily_secondary_attack_rate,
                            fatality_pct)


    simulation = Simulation(population, group_test, test_day_frequency)
    return simulation

def summarize(simulation):
    print("Total number of tests performed over {} days: {}".format(simulation.current_day,
                                                                   simulation.cumulative_tests_to_date))
    days = range(simulation.current_day)
    cumulative_infected_pct = [simulation.recorded_data[day]['cumulative_infected_fraction'] for day in days]
    quarantine_pct = [simulation.recorded_data[day]['in_quarantine_fraction'] for day in days]
    plt.figure(figsize=(10,6))
    plt.ylim((-0.1,1.1))
    plt.plot(days, cumulative_infected_pct, label="Cumulative Fraction of Infected Population")
    plt.plot(days, quarantine_pct, label="Fraction of Population in Quarantine")
    plt.legend(loc='best')
    plt.show()
    
def run(simulation, number_of_days):
    for _ in range(number_of_days):
        simulation.step()
