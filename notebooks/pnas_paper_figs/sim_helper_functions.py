import * from lhs_space_transformations

def run_parallel_sims_lhs_space(uncertainty_points, output_folder):
    run_parallel_sims_params_dict(
        [uncertainty_point_to_params_dict(uncertainty_point) for uncertainty_point in uncertainty_points],
        output_folder)

def run_parallel_sims_params_dict(params_dicts, output_folder):
