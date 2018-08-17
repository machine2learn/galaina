import configparser
import os
from typing import Dict

info = 'info'
copula_factor_algorithm = 'copula_factor_algorithm'
edge_weight_algorithm = 'edge_weight_algorithm'
pc_algorithm = 'pc_algorithm'
plot_and_display = 'plot_and_display'
input_paths_and_related_parameters = 'input_paths_and_related_parameters'
output_paths = 'output_paths'


def abs_path_of(rel_path):
    return os.path.join(os.path.dirname(__file__), rel_path)


class CustomConfigParser(configparser.ConfigParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_as_slice(self, *args, **kwargs):
        raw_get = self.get(*args, **kwargs)
        if ":" in raw_get:
            return slice(*map(int, raw_get.split(":")))
        else:
            return int(raw_get)

    def get_rel_path(self, *args, **kwargs):
        raw_get = self.get(*args, **kwargs)
        if not raw_get:
            return ""
        if raw_get.startswith('/'):
            return raw_get

        return abs_path_of(raw_get)


    def _from_info(self, param):
        return self.get(info, param)

    def _from_copula_factor(self, param):
        return self.get(copula_factor_algorithm, param)

    def _from_edge_weight(self, param):
        return self.get(edge_weight_algorithm, param)

    def _from_pc_algorithm(self, param):
        return self.get(pc_algorithm, param)

    def _from_plot_and_display(self, param):
        return self.get(plot_and_display, param)

    def _from_input_paths_and_related_parameters(self, param):
        return self.get(input_paths_and_related_parameters, param)

    def _from_output_paths(self, param):
        return self.get(output_paths, param)

    def get_path(self):
        return self._from_info('config_path')

    def get_name(self):
        return self._from_info('config_name')

    def get_gibbs_sampling_n(self):
        return self._from_copula_factor('gibbs_sampling_n')

    def get_gibbs_burn_in_n(self):
        return self._from_copula_factor('gibbs_burn_in_n')

    def get_gibbs_first_random_seed_n(self):
        return self._from_copula_factor('gibbs_first_random_seed_n')

    def get_gibbs_random_seed_update_parameter_n(self):
        return self._from_copula_factor('gibbs_random_seed_update_parameter_n')

    def get_bootstrap_n(self):
        return self._from_edge_weight('bootstrap_n')

    def get_bootstrap_first_random_seed_n(self):
        return self._from_edge_weight('bootstrap_first_random_seed_n')

    def get_bootstrap_random_seed_update_parameter_n(self):
        return self._from_edge_weight('bootstrap_random_seed_update_parameter_n')

    def get_causal_discovery_observation_n(self):
        return self._from_pc_algorithm('causal_discovery_observation_n')

    def get_indeptest(self):
        return self._from_pc_algorithm('indeptest')

    def get_alpha(self):
        return self._from_pc_algorithm('alpha')

    def get_numcores(self):
        return self._from_pc_algorithm('numcores')

    def get_verbose(self):
        return self._from_pc_algorithm('verbose')

    def get_nadelete(self):
        return self._from_pc_algorithm('nadelete')

    def get_m_max(self):
        return self._from_pc_algorithm('m_max')

    def get_u2pd(self):
        return self._from_pc_algorithm('u2pd')

    def get_skel_method(self):
        return self._from_pc_algorithm('skel_method')

    def get_conservative(self):
        return self._from_pc_algorithm('conservative')

    def get_maj_rule(self):
        return self._from_pc_algorithm('maj_rule')

    def get_solve_confl(self):
        return self._from_pc_algorithm('solve_confl')

    def get_fixedgaps(self):
        return self._from_pc_algorithm('fixedgaps')

    def get_fixededges(self):
        return self._from_pc_algorithm('fixededges')

    def get_core_plot_title_str(self):
        return self._from_plot_and_display('core_plot_title_str')

    def get_input_data_ls(self):
        return self._from_input_paths_and_related_parameters('input_data_ls')

    def get_input_factor_ls(self):
        return self._from_input_paths_and_related_parameters('input_factor_ls')

    def get_ls_separator_str(self):
        return self._from_input_paths_and_related_parameters('ls_separator_str')

    def get_column_separator_str(self):
        return self._from_input_paths_and_related_parameters('column_separator_str')

    def get_input_path_directed_edges_blacklist(self):
        return self._from_input_paths_and_related_parameters('input_path_directed_edges_blacklist')

    def get_output_path_merged_data(self):
        return self._from_output_paths('output_path_merged_data')

    def get_output_path_merged_factor_model_table(self):
        return self._from_output_paths('output_path_merged_factor_model_table')

    def get_output_path_merged_factor_model_loading(self):
        return self._from_output_paths('output_path_merged_factor_model_loading')

    def get_output_path_fig(self):
        return self._from_output_paths('output_path_fig')

    def get_output_path_log(self):
        return self._from_output_paths('output_path_log')

    def get_output_path_suffstat(self):
        return self._from_output_paths('output_path_suffstat')

    def get_output_path_pc_algo_obj(self):
        return self._from_output_paths('output_path_pc_algo_obj')

    def get_output_path_bn_strength_obj(self):
        return self._from_output_paths('output_path_bn_strength_obj')

    def get_output_path_avg_bn_obj(self):
        return self._from_output_paths('output_path_avg_bn_obj')


def read_config(path):
    config = CustomConfigParser(inline_comment_prefixes=['#'], interpolation=configparser.ExtendedInterpolation())
    config.read(path)
    return config
