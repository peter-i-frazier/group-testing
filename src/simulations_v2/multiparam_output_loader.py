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
                params = None
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



