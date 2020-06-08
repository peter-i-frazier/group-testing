"""
Implements a class to automate loading/processing the output directories
created by the run_sensitivity script

The MultiParamOutputLoader class takes a top-level directory created by the run_sensitivity
script as input, and performs the following tasks:
    * Each parameter scenario is identified by looking for the subfolders of the top-level folder
        that is passed as input
        * the name of each subfolder is taken to be the name of the corresponding parameter scenario
    * For each parameter scenario subfolder:
        * the scenario-parameter dill file is loaded from the scn_params.dill file
            * if the file is not found then an exception is thrown
        * then, each subfolder of the scenario-parameter folder is iterated over, and
          assumed to contain multiple CSVs each corresponding to one simulation trajectory
          for the particular parameter configuration considered in that subfolder
        * the subfolder-specific parameter configuration is assumed to be specified in the param_specifier.yaml
          file

The resulting object, sim_output = MultiParamOutputLoader(sim_dir), has the following variables:
    * sim_output.param_scenarios: a list of names of the different parameter scenarios considered in the sim
    * sim_output.scn_params: a dictionary mapping param scenario names to the corresponding parameters dictionary object
    * sim_output.sim_results: a dictionary mapping pram scenario names to a dictionary of sim results
    * sim_output.sim_results[scenario_name]: a dictionary mapping varied-parameter-value tuples to a list of dataframes
                                             in one-to-one correspondence with the resulting trajectories
                                             The order of valalues in the tuple-key corresponnds to the order of variables
                                             in the sim_output.varying_params list
    * sim_output.varying_params: a list of parameters that were varied in the multiparameter simulation
"""

import yaml
import os
import pandas as pd
import dill


class MultiParamOutputLoader:
    def __init__(self, sim_dir):
        self.param_scenarios = []
        self.sim_results = {}
        self.scn_params = {}
        self.varying_params = None
        subfolders = [f.path for f in os.scandir(sim_dir) if f.is_dir()]
        for folder in subfolders:
            scn_name = folder.split('/')[-1]
            dill_path = '{}/scn_params.dill'.format(folder)
            if os.path.exists(dill_path) and os.path.isfile(dill_path):
                with open(dill_path, 'rb') as params_file:
                    params = dill.load(params_file)
            else:
                raise(Exception("Could not find params dill file in the subfolder {}".format(folder)))
            self.scn_params[scn_name] = params
            self.param_scenarios.append(scn_name)
            self.sim_results[scn_name] = {}
            self.load_param_scenario(folder, scn_name)

    def load_param_scenario(self, folder, scn_name):
        subfolders = [f.path for f in os.scandir(folder) if f.is_dir()]
        for subfolder in subfolders:
            with open('{}/param_specifier.yaml'.format(subfolder)) as f:
                param_specifier = yaml.load(f)
            if self.varying_params == None:
                self.varying_params = list(param_specifier.keys())
            assert(set(self.varying_params) == set(param_specifier.keys()))
            param_specifier_key = [param_specifier[param] for param in self.varying_params]
            param_specifier_key = tuple(param_specifier_key)

            self.sim_results[scn_name][param_specifier_key] = []
            for f in os.scandir(subfolder):
                if f.path.split('.')[-1] == 'csv':
                    self.sim_results[scn_name][param_specifier_key].append(
                                                        pd.read_csv(f.path, index_col=0))



