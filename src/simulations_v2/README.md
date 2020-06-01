# Stochastic Population-Level Simulation

This folder implements the population-level stochastic simulation that is used in our analysis.
The main simulation logic is contained in the file `stochastic_simulation.py`.

For a complete description of the simulation model see the associated writeup.

## Running a Sensitivity Analysis

A sensitivity analysis fixes a collection of parameter scenarios and plots the distribution
of hospitalizations and infections for each scenario as a single parameter is varied. 

We include the scripts `run_sensitivity.py` and `run_sensitivity.sh` to streamline the process
of running a sensitivity analysis. Here is a summary of the options for `run_sensitivity.py` and what they do:
* `-o` `--outputdir` Specify directory to store simulation output.  Defaults to `/nfs01/covid_sims/`.
* `-V` `--verbose` Include verbose output.  Defaults to `False`.
* `-s` `--scenarios` Specify a list of parameter scenarios to use in the sensitivity analysis.
	Each parameter scenario should be specified as a yaml config file that can be parsed by the `load_params()`
	function in the file `load_params.py`.  For a pre-defined set of parameter scenarios, see the `params/` folder.
* `-p` `--param-to-vary` Specify which parameter should be varied in the sensitivity analysis.  Right now the script
	only supports the following parameters:
	* `contact_tracing_isolations`
	* `expected_contacts_per_day`
	* `mild_symptoms_daily_self_report_p`
	* `severe_symptoms_daily_self_report_p`
	* `initial_ID_prevalence`
	* `exposed_infection_p`
	* `asymptomatic_p`
	* `contact_tracing_delay`
	* `test_protocol_QFNR`
	* `test_population_fraction`
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

### Specifying Script Arguments in a File

Since the number of arguments can get quite large, and since we might want to run the same sensitivity analysis many times
and thus would like a method to recall what parameters were used for a particular simulation, we include a bash script
`run_sensitivity.sh` which calls the `run_sensitivity.py` script using parameters that are specified in a particular file.

Suppose you had the following text inside a file `initial_ID_prevalence.txt`:
```
-s params/fall_nominal.yaml params/fall_nominal_testing.yaml
-p initial_ID_prevalence
-v 0.0001 0.0005 0.001 0.0025
```

Then running `bash run_sensitivity.sh initial_ID_prevalence.txt` would be equivalent to running the python script with the 
parameters specified above.  

You can also specify additional parameters to the bash script that are not included in the file.  For example,
`bash run_sensitivity.sh initial_ID_prevalence.txt -n 15 -t 100` is equivalent to the first example we gave.

For more parameter configuration scripts, see the `sensitivity_config/` folder.
