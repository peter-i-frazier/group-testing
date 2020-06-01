import yaml
from subdivide_severity import subdivide_severity

def load_age_severity_params(param_file):


# reads stochastic-simulation parameters from a yaml config file
# supports depence between config files, so that one param file
# can point to another file and params from the pointed-to-file
# are loaded first
def load_params(param_file, param_file_stack=[]):
    with open(param_file) as f:
        params = yaml.load(f)

    if '_age_severity_config' in params:
        age_sev_file = params['_age_severity_config']
        age_sev_params = load_age_sev_params(age_sev_file)

    if '_inherit_config' in params:


