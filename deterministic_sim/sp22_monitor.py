import sys
import numpy as np
import json
import yaml
from groups import population
from sim_helper import sim_test_regime, sim_test_strategy
from sp22_strategies import (no_testing_strategy, arrival_testing_strategy,
                             surge_testing_strategy)
import plotting


def main(yaml_file='nominal.yaml', out_file='sp22_sim.png', **kwargs):

    # =======================
    # [Initialize Parameters]
    # =======================

    params = yaml.safe_load(open(yaml_file, "r"))
    params.update(kwargs)

    # convert to numpy matrix
    params["meta_matrix"] = \
        np.array([list(row.values()) for row in params["meta_matrix"].values()])

    # If set, this replaces the detailed description of parameters in the plot
    # with a simple summary
    if 'simple_param_summary' in params:
        SIMPLE_PARAM_SUMMARY = params['simple_param_summary']
    else:
        SIMPLE_PARAM_SUMMARY = None

    # ==================================
    # [Run] Compare a list of strategies
    # ==================================

    trajectories = [sim_test_strategy(params, surge_testing_strategy(params), 'purple')]

    # =================
    # [Plot] Make plots
    # ==================

    plotting.plot_arrival_linear(out_file, trajectories, 7)


def usage():
    ''' Print usage message '''
    print('Usage:')
    print('python sp22_monitor.py [-h] [--help] [--yaml=file.yaml] [--out=file.png] [yaml-overrides]')
    print('--yaml=file.yaml uses file.yaml instead of nominal.yaml')
    print('--outfile=file.png saves the plot to file.png instead of sp22_sim.png')
    print('--help and -h print this message and exit')
    print('yaml-overrides replace parameters in the yaml file and are of the form parameter_name=value')
    print('Example: python sp22_sim.py --yaml=nominal_ug.yaml --out=nominal_ug.png T=10')


if __name__ == "__main__":

    # Default values for arguments not specified in the yaml
    yaml_file = 'nominal.yaml'
    out_file = 'sp22_monitor.png' # The filename for the plots

    # Parameters from the yaml file to override
    override_params = {}

    for arg in sys.argv[1:]:
        # Handle arguments without an equal sign
        if arg == '-h' or arg == '--help':
            usage()
            exit(0)
            continue

        # If we get to this point, the argument should have the form key=value
        # Arguments that are not part of the things specified in the YAML file
        assert(len(arg.split('='))==2)
        k,v = arg.split('=') # key, value

        # Arguments that are not YAML overrides
        if k == "--yaml":
            yaml_file = str(v)
        elif k == '--out':
            out_file = str(v)

        # YAML overrides
        elif k == "T":
            # This yaml over
            override_params[k] = int(v)
        elif k == 'simple_param_summary':
            override_params[k] = v # string
        else:
            override_params[k] = float(v)

    main(yaml_file, out_file, **override_params)
