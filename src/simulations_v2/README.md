# Stochastic Population-Level Simulation

This folder implements the population-level stochastic simulation that is used in our analysis.
The main simulation logic is contained in the file `stochastic\_simulation.py`.

For a complete description of the simulation model see the associated writeup.

## Running a Sensitivity Analysis

A sensitivity analysis fixes a collection of parameter scenarios and plots the distribution
of hospitalizations and infections for each scenario as a single parameter is varied. 

We include the scripts `run\_sensitivity.py` and `run\_sensitivity.sh` to streamline the process
of running a sensitivity analysis. Here is a summary of the options for `run\_sensitivity.py` and what they do:
* `-o` `--outputdir` Specify directory to store simulation output.  Defaults to `/nfs01/covid\_sims/`.
* `-V` `--verbose` Include verbose output.  Defaults to `False`.
* `-s` `--scenarios` Specify a list of parameter scenarios to use in the sensitivity analysis.
	Each parameter scenario should be specified as a yaml config file that can be parsed by the `load\_params()`
	function in the file `load\_params.py`.  For a pre-defined set of parameter scenarios, see the `params/` folder.
* `-p` `--param-to-vary` Specify which parameter should be varied in the sensitivity analysis.  Right now the script
	only supports the following parameters:
** `contact\_tracing\_isolations`
** `expected\_contacts\_per\_day`
** `mild\_symptoms\_daily\_self\_report\_p`
** `severe\_symptoms\_daily\_self\_report\_p`
** `initial\_ID\_prevalence`
** `exposed\_infection\_p`
** `asymptomatic\_p`
** `contact\_tracing\_delay`
** `test\_protocol\_QFNR`
** `test\_population\_fraction`
* `-v` `--values` Specify the list of values that the varying parameter should take
* `-n` `--ntrajectories` Specify the total number of trajectories to simulate for each (scenario, parameter value) pair. Defaults to 500.
* `-t` `--time-horizon` Specify the time horizon (in days) that each trajectory should use. Defaults to 112 days.
* `-f` `--fig-dir` Specify a directory where the sensitivity plots should be saved.  Defaults to `outputdir`.

As an example, consider the following command:
```
python3 run_sensitivity.py -s params/fall_nominal.yaml params/fall_nominal_testing.yaml -p initial_ID_prevalence -v 0.0001 0.0005 0.001 0.0025 -n 15 -t 100
```

This uses the `-s` flag to specify the `params/fall_nominal_testing.yaml` and `params/fall_nominal.yaml` parameter scenarios,
it uses the `-p` flag to specify that `initial_ID_prevalence` is the parameter which is varying, and it uses the `-v` flag 
to specify that we should simulate outcomes when the variable takes on the following values: `0.0001 0.0005 0.001 0.0025`.
It also uses the `-n` flag to specify only 15 trajectories, which is too small to get meaningful observations, but it is small
enough so that the sims complete quickly.   

Here is the output when I run that command:
```
Running simulations for 2 scenarios and 4 parameter values across 8 separate processes.
Results being saved in output directory /nfs01/covid_sims//1591042932-initial_ID_prevalence.
Waiting for simulations to finish...
Simulations done. Generating plots now...
Saved plots to directory /nfs01/covid_sims//1591042932-initial_ID_prevalence
```
