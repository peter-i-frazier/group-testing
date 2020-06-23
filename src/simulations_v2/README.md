# Stochastic Population-Level Simulation

This folder implements the population-level stochastic simulation to quantify the spread
of a virus (e.g. COVID-19) under contact-tracing and frequent testing of the entire population. 

While we have tried hard to make the code as configurable as possible, all existing parameter
files (see `params` folder and discussion below) are tailored to modeling outcomes for reopening
Cornell's Ithaca campus in the fall.  If you'd like to use this code to run simulations for a
different institution or under a different set of parameters: please feel free to contact us,
we'd be happy to help get that set up.

For a complete description of the simulation model and assumptions see the associated writeup 
located in the `reports/` folder.

**Disclaimer:** This codebase is rapidly evolving and some of the information below might soon become
out-of-date.  We will try to keep this README accurate, but if you notice anything is amiss please let
us know.

## Important files in this directory
* `stochastic_simulation.py` 
    * This file implements the main simulation logic via the class
     `StochasticSimulation`.  This class takes as input a large parameters dictionary; to simplify
     these parameters, we include a script to instantiate the parameters dictionary from a simple YAML file
     (see below).   
    * The main external-facing function in this class is `run_new_trajectory(T)`, which resets all the
     state variables to their initial values and then simulates a single infection trajectory over `T` days.
     This function returns a pandas dictionary which contains a row for each day and includes data
     summarizing the simulated trajectory.
* `load_params.py`
    * Implements functions to load simulation model parameters from configurable YAML files.
    * For example parameter config files see the files: `params/fall_nominal.yaml` and 
        `params/fall_nominal_testing.yaml`.
    * The main external facing function is `load_params(param_file)` which takes a YAML-file path
     argument `param_file` and returns a corresponding parameters dictionary that can be used to instantiate
     a `StochasticSimulation` object
* `run_sensitivity.py` and `run_sensitivity.sh`
    * Helper scripts to run multiple simulations over multiple parameter configurations using multiprocessing.
    * Supports single-parameter and multi-parameter sensitivity analysis.  In the case of a single-parameter
      sensitivity-analysis, automatically generates infection and hospitalization percentage plots (i.e. automatically
      generates the plots which appear in the Results section of our report).
    * The folder `sensitivity_config` contains many `.txt` files which serve as parameters for the bash
      version of this script.  E.g. the following command `bash run_sensitivity.sh sensitivity_config/asymptomatic_p.txt`
      kicks off a single-parameter sensitivity analysis which varies the asymptomatic-infection-fraction parameter.
    * Saves all the resulting trajectory data-frames to a configurable output directory location.
    * Uses `dill` to save a copy of the parameters dictionary that was used to run each trajectory.
    * See below section for more detail about running sensitivity analyses.
* `multiparam_output_loader.py`
    * Implements a class to load all the trajectory data saved to the output folder by the `run_sensivity` scripts.
    * Main external facing class is `MultiParamOutputLoader`; see the comments at the top of this file for how to
     use it.
    * (The name is a bit of a misnomer -- this script should be used to load parameters regardless of whether you did
        a single-parameter or multi-parameter sensitivity analysis).


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
* `-p` `--param-to-vary` Add a parameter to vary in the sensitivity analysis.  
       One or more parameters can be added.  At the time of writing
       the script supports the following parameters:
	* `contact_tracing_isolations`
	* `expected_contacts_per_day`
	* `symptomatic_daily_self_report_p`
	* `asymptomatic_daily_self_report_p`
	* `initial_ID_prevalence`
	* `exposed_infection_p`
	* `asymptomatic_p`
	* `contact_tracing_delay`
	* `test_protocol_QFNR`
	* `test_population_fraction`
    * `daily_outside_infection_p`
    * `initial_E_count`
    * `initial_pre_ID_count`
    * `initial_ID_count`
    * `initial_SyID_mild_count`
    * `initial_SyID_severe_count`
    * `population_size`
* `-v` `--values` Specify the list of values that the varying parameter should take.  
        The number of times `-v` is specified should be the same number of times that `-p` is specified and the order of
        the parameter values should match the order of the varying parameters.
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
