import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
from matplotlib import ticker
from analysis_helpers import load_sim_dir
from textwrap import wrap
from multiparam_output_loader import MultiParamOutputLoader

plot_labels = {
    'daily_contacts': 'Average Contacts per Person per Day',
    'severe_self_reporting': 'Daily Symptomatic Self-Reporting Likelihood(%)',
    'prevalence': 'Initial Prevalence (% of Total Population)',
    'mild_self_reporting': 'Daily Asymptomatic Self-Reporting Likelihood(%)',
    'exposed_infection_p': 'Transmission Likelihood per Contact (%)',
    'contact_isolations': 'Isolations per Contact Trace',
    'contact_delay': 'Contact Trace Delay (Days)',
    'asymptomatic_p': 'Percentage of Cases Asymptomatic (%)',
    'test_protocol_QFNR': 'Testing False-Negative Rate',
    'test_population_fraction': 'Percentage of Population Tested Daily',
    'daily_outside_infection_p': 'Daily Individual Probability of Outside Infection (%)'
}

# normalize x-axis to be a percentage 
normalize_params = {
    'daily_contacts': False,
    'severe_self_reporting': True,
    'prevalence': True,
    'mild_self_reporting': True,
    'exposed_infection_p': True,
    'contact_isolations': False,
    'contact_delay': False,
    'asymptomatic_p': True,
    'default': False,
    'test_protocol_QFNR': True,
    'test_population_fraction': True,
    'daily_outside_infection_p': True,
}

# put the x-axis on a log-scale
plot_log_scale = {
    'daily_contacts': True,
    'severe_self_reporting': False,
    'prevalence': True,
    'mild_self_reporting': False,
    'exposed_infection_p': True,
    'contact_isolations': False,
    'contact_delay': False,
    'asymptomatic_p': True,
    'test_protocol_QFNR': False,
    'test_population_fraction': True,
    'daily_outside_infection_p': True
}

# rount x-axis labels to int
use_x_int_labels = {
    'daily_contacts': True,
    'severe_self_reporting': True,
    'prevalence': False,
    'mild_self_reporting': True,
    'exposed_infection_p': False,
    'contact_isolations': False,
    'contact_delay': True,
    'asymptomatic_p': False,
    'test_protocol_QFNR': False,
    'test_population_fraction': False,
    'daily_outside_infection_p': False
}

key_mapping = {
    'contact_tracing_isolations': 'contact_isolations',
    'expected_contacts_per_day': 'daily_contacts',
    'asymptomatic_daily_self_report_p': 'mild_self_reporting',
    'symptomatic_daily_self_report_p': 'severe_self_reporting',
    'initial_ID_prevalence': 'prevalence',
    'exposed_infection_p': 'exposed_infection_p',
    'asymptomatic_p': 'asymptomatic_p',
    'contact_tracing_delay': 'contact_delay',
    'test_protocol_QFNR': 'test_protocol_QFNR',
    'test_population_fraction': 'test_population_fraction',
    'daily_outside_infection_p': 'daily_outside_infection_p'
}

def plot_from_folder(folder, savefig_dir):
    sim_output = MultiParamOutputLoader(folder)
    if len(sim_output.varying_params) != 1:
        raise(Exception("Can only make sensitivity plots for single-parameter varying simulations"))

    param_varying = sim_output.varying_params[0]
    if param_varying in key_mapping:
        sublabel = key_mapping[param_varying]
        plot_label = plot_labels[sublabel]
    else:
        raise(Exception("could not find parameter {} in the list of acceptable keys".format(param_varying)))

    plot_many_dfs_quantiles(sim_output.sim_results, 
                            cum_severe_quantiles,
                            normalize_params[sublabel], x_log_scale=plot_log_scale[sublabel],
                     xlabel=plot_labels[sublabel], ylabel="Percentage of Population Hospitalized",
                     title="Nominal Parameters: Hospitalization Percentage vs. {}".format(plot_label), 
                               q_low=0.1, q_high=0.9, alpha=0.1, y_min = 0, y_max=5, y_log_scale=True,
                               savefig_path="{}/hospitalization_{}.pdf".format(savefig_dir, param_varying), 
                             use_x_int_labels=use_x_int_labels[sublabel])
    
    plot_many_dfs_quantiles(sim_output.sim_results, 
                            cum_infection_quantiles,
                            normalize_params[sublabel], x_log_scale=plot_log_scale[sublabel],
                     xlabel=plot_labels[sublabel], ylabel="Percentage of Population Infected",
                     title="Nominal Parameters: Infection Percentage vs. {}".format(plot_label), 
                               q_low=0.1, q_high=0.9, alpha=0.1, y_min = 0, y_max=5, y_log_scale=True,
                               savefig_path="{}/infection_{}.pdf".format(savefig_dir, param_varying), 
                             use_x_int_labels=use_x_int_labels[sublabel])
    plot_many_dfs_quantiles(sim_output.sim_results, 
                            cum_non_outside_infection_quantiles,
                            normalize_params[sublabel], x_log_scale=plot_log_scale[sublabel],
                     xlabel=plot_labels[sublabel], ylabel="Percentage of Population Infected (from Non-Outside Interaction)",
                     title="Nominal Parameters: Non-Outside Infection Percentage vs. {}".format(plot_label), 
                               q_low=0.1, q_high=0.9, alpha=0.1, y_min = 0, y_max=5, y_log_scale=True,
                               savefig_path="{}/infection_non_outside_{}.pdf".format(savefig_dir, param_varying), 
                             use_x_int_labels=use_x_int_labels[sublabel])

def plot_from_folders_legacy(folder_map, param_varying, savefig_dir):
    sim_map = {scn_name: load_sim_dir(scn_folder, verbose=False) for scn_name, scn_folder in folder_map.items()}

    if param_varying in key_mapping:
        sublabel = key_mapping[param_varying]
        plot_label = plot_labels[sublabel]
    else:
        raise(Exception("could not find parameter {} in the list of acceptable keys".format(param_varying)))

    plot_many_dfs_quantiles(sim_map, 
                            cum_severe_quantiles,
                            normalize_params[sublabel], x_log_scale=plot_log_scale[sublabel],
                     xlabel=plot_labels[sublabel], ylabel="Percentage of Population Hospitalized",
                     title="Nominal Parameters: Hospitalization Percentage vs. {}".format(plot_label), 
                               q_low=0.1, q_high=0.9, alpha=0.1, y_min = 0, y_max=5, y_log_scale=True,
                               savefig_path="{}/hospitalization_{}.pdf".format(savefig_dir, param_varying), 
                             use_x_int_labels=use_x_int_labels[sublabel])
    
    plot_many_dfs_quantiles(sim_map, 
                            cum_infection_quantiles,
                            normalize_params[sublabel], x_log_scale=plot_log_scale[sublabel],
                     xlabel=plot_labels[sublabel], ylabel="Percentage of Population Infected",
                     title="Nominal Parameters: Infection Percentage vs. {}".format(plot_label), 
                               q_low=0.1, q_high=0.9, alpha=0.1, y_min = 0, y_max=5, y_log_scale=True,
                               savefig_path="{}/infection_{}.pdf".format(savefig_dir, param_varying), 
                             use_x_int_labels=use_x_int_labels[sublabel])
    plot_many_dfs_quantiles(sim_map, 
                            cum_non_outside_infection_quantiles,
                            normalize_params[sublabel], x_log_scale=plot_log_scale[sublabel],
                     xlabel=plot_labels[sublabel], ylabel="Percentage of Population Infected (from Non-Outside Interaction)",
                     title="Nominal Parameters: Non-Outside Infection Percentage vs. {}".format(plot_label), 
                               q_low=0.1, q_high=0.9, alpha=0.1, y_min = 0, y_max=5, y_log_scale=True,
                               savefig_path="{}/infection_non_outside_{}.pdf".format(savefig_dir, param_varying), 
                             use_x_int_labels=use_x_int_labels[sublabel])

def modify_tick_labels(labels, axis, use_int_labels, max_str_length=7):
    assert(axis in ['x', 'y'])
    new_labels = []
    for label in labels:
        txt = label.get_text()
        if txt == '':
            new_labels.append('')
        else: 
            if axis == 'x':
                val = label._x
            else:
                val = label._y
            
            if use_int_labels:
                newtxt = str(int(val))
            else:
                newtxt = str(float(val))

            if len(newtxt) >= max_str_length:
                newtxt = newtxt[0:max_str_length]
            
            new_labels.append(newtxt)
    return new_labels
            

def extract_severities(dfs, hospitalizations_only=True, subtract_outside_infections=False):
    severities = []
    for df in dfs:
        all_cols = set(df.columns)
        new_cols = set(['cumulative_mild', 
                        'cumulative_severe', 
                        'cumulative_outside_infections',
                        'severity_0', 
                        'severity_1', 
                        'severity_2', 
                        'severity_3'])
        main_cols = all_cols - new_cols
        subdf = df[list(main_cols)]
        popsize = sum(subdf.iloc[0])
        
        if hospitalizations_only:
            subdf = df[['severity_2', 'severity_3']]
        else:
            subdf = df[['severity_0', 'severity_1', 'severity_2', 'severity_3']]

        severe = sum(subdf.iloc[subdf.shape[0]-1])

        if subtract_outside_infections:
            subdf = df[['cumulative_outside_infections']]
            outside_infections = sum(subdf.iloc[subdf.shape[0] - 1])
            severe = severe - outside_infections
        
        severities.append(100 * severe / popsize)
    return severities
def average_cumulative_severe(dfs):
    severities = []
    for df in dfs:
        all_cols = set(df.columns)
        new_cols = set(['cumulative_outside_infections',
                        'cumulative_mild', 
                        'cumulative_severe', 
                        'severity_0', 
                        'severity_1', 
                        'severity_2', 
                        'severity_3'])
        main_cols = all_cols - new_cols
        subdf = df[list(main_cols)]
        popsize = sum(subdf.iloc[0])
        
        subdf = df[['severity_2', 'severity_3']]
        severe = sum(subdf.iloc[subdf.shape[0]-1])
        
        severities.append(100 * severe / popsize)
    return np.mean(severities)

def avg_cum_severe_quantile(dfs, q_low, q_high):
    severities = extract_severities(dfs)
    return np.quantile(severities, q_low), np.mean(severities), np.quantile(severities, q_high)

def cum_severe_quantiles(dfs, q_low, q_high):
    severities = extract_severities(dfs)
    return np.quantile(severities, q_low), np.quantile(severities, 0.5), np.quantile(severities, q_high)

def cum_infection_quantiles(dfs, q_low, q_high):
    severities = extract_severities(dfs, hospitalizations_only=False)
    return np.quantile(severities, q_low), np.quantile(severities, 0.5), np.quantile(severities, q_high)

def cum_non_outside_infection_quantiles(dfs, q_low, q_high):
    severities = extract_severities(dfs, hospitalizations_only=False, subtract_outside_infections=True)
    return np.quantile(severities, q_low), np.quantile(severities, 0.5), np.quantile(severities, q_high)

def plot_many_dfs(sim_output_dict, yaxisfn, ylabel="", xlabel="", title="", figsize=(10,6)):
    plt.figure(figsize=figsize)
    for sim_label, sim_output in sim_output_dict.items():
        xs = []
        ys = []
        for sim_parameter_name, dfs in sim_output.items():
            # compute x-value assuming that sim_param_name is of form 'varied_param_name.value'
            param_val = float('.'.join(sim_parameter_name.split('.')[1:]))
            xs.append(param_val)
            # yaxisfn is a function that takes in a list of trajectory dataframes and
            # produces an output metric
            ys.append(yaxisfn(dfs))
        plt.plot(xs, ys, marker='o', label=sim_label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(loc='best')
    plt.show()


def truncate(val, y_min, y_max):
    return max(min(val, y_max), y_min)
    
def plot_many_dfs_quantiles(sim_output_dict, yaxisfn, normalize_x_axis, x_log_scale=False, 
                            y_log_scale=False,
                            q_low=0.05, q_high=0.95, 
                            y_min = 0, y_max = 5,
                            ylabel="", xlabel="", title="", 
                            figsize=(10,6), alpha=0.1, color=None, savefig_path = None,
                            use_x_int_labels=False, use_y_int_labels=False):
    # assn: yaxisfn(dfs) returns a tuple (q_low_val, avg, q_high_val)
    epsilon=1e-6
    plt.figure(figsize=figsize)
    if x_log_scale:
        plt.xscale("log")
    if y_log_scale:
        plt.yscale("log")
        #_, _, ymin, _ = plt.axis()
        #plt.ylim(bottom=max(1e-2, ymin))
    else:
        plt.ylim(y_min, y_max)
    
    for sim_label, sim_output in sim_output_dict.items():
        xs = []
        ys = []
        q_low_vals = []
        q_high_vals = []
        for sim_parameter_name, dfs in sim_output.items():
            # compute x-value assuming that sim_param_name is of form 'varied_param_name-value'
            if type(sim_parameter_name) == str:
                param_val = float('-'.join(sim_parameter_name.split('-')[1:]))
            elif type(sim_parameter_name) == tuple:
                assert(len(sim_parameter_name)) == 1
                param_val = sim_parameter_name[0]
            else:
                raise(Exception("unrecognized type for sim_parameter_name"))
            if normalize_x_axis:
                xs.append(param_val * 100)
            else:
                xs.append(param_val)
            # yaxisfn is a function that takes in a list of trajectory dataframes and
            # produces an output metric
            q_low_val, avg, q_high_val = yaxisfn(dfs, q_low, q_high)
            if not y_log_scale:
                ys.append(truncate(avg, y_min, y_max))
                q_low_vals.append(truncate(q_low_val, y_min, y_max))
                q_high_vals.append(truncate(q_high_val, y_min, y_max))
            else:
                ys.append(avg + epsilon)
                q_low_vals.append(q_low_val + epsilon)
                q_high_vals.append(q_high_val + epsilon)
        
        sorted_points = sorted(zip(xs, ys, q_low_vals, q_high_vals), key=lambda x: x[0])
        xs = [p[0] for p in sorted_points]
        ys = [p[1] for p in sorted_points]
        q_low_vals = [p[2] for p in sorted_points]
        q_high_vals = [p[3] for p in sorted_points]
        
        if color == None:
            plt.plot(xs, ys, marker='o', label=sim_label)
            plt.fill_between(xs, q_low_vals, q_high_vals, alpha=alpha)
        else:
            plt.plot(xs, ys, marker='o', label=sim_label, color=color)
            plt.fill_between(xs, q_low_vals, q_high_vals, alpha=alpha, color=color)
    
    plt.xlabel("\n".join(wrap(xlabel, 60)))
    plt.ylabel("\n".join(wrap(ylabel, 60)))
    plt.title("\n".join(wrap(title, 60)))
    plt.legend(loc='best')
    
    plt.draw()
    
    ax = plt.gca()
    
    
    labels = ax.get_xmajorticklabels()
    new_labels = modify_tick_labels(labels, "x", use_x_int_labels)
    ax.set_xticklabels(new_labels)
    
    labels = ax.get_xminorticklabels()
    new_labels = modify_tick_labels(labels, "x", use_x_int_labels)
    ax.set_xticklabels(new_labels, minor=True)
    
    labels = ax.get_ymajorticklabels()
    new_labels = modify_tick_labels(labels, "y", use_y_int_labels)
    ax.set_yticklabels(new_labels)
    
    labels = ax.get_yminorticklabels()
    new_labels = modify_tick_labels(labels, "y", use_y_int_labels)
    ax.set_yticklabels(new_labels, minor=True)
    
    plt.tight_layout()
    if savefig_path:
        plt.savefig(savefig_path)
