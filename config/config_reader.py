import configparser
import os
from typing import Dict

INFO = 'INFO'
COPULA_FACTOR = 'COPULA_FACTOR'
EDGE_WEIGHT = 'EDGE_WEIGHT'
PC_ALGORITHM = 'PC_ALGORITHM'
PLOT_AND_DISPLAY = 'PLOT_AND_DISPLAY'
input_paths_and_related_parameters = 'input_paths_and_related_parameters'


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
        return self.get(INFO, param)

    def _from_copula_factor(self, param):
        return self.get(COPULA_FACTOR, param)

    def _from_edge_weight(self, param):
        return self.get(EDGE_WEIGHT, param)

    def _from_pc_algorithm(self, param):
        return self.get(PC_ALGORITHM, param)

    def _from_plot_and_display(self, param):
        return self.get(PLOT_AND_DISPLAY, param)

    def _from_input_paths_and_related_parameters(self, param):
        return self.get(input_paths_and_related_parameters, param)

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

    # def train_batch_size(self) -> int:
    #     return int(self._from_training('batch_size'))
    #
    # def learning_rate(self) -> int:
    #     return int(self._from_training('learning_rate'))
    #
    # def validation_batch_size(self) -> int:
    #     return int(self._from_training('validation_batch_size'))
    #
    # def optimizer(self) -> str:
    #     return self._from_training('optimizer')
    #
    # def l1_reqularization(self) -> float:
    #     return float(self._from_training('l1_regularization'))
    #
    # def l2_reqularization(self) -> float:
    #     return float(self._from_training('l2_regularization'))
    #
    # def num_epochs(self) -> int:
    #     return int(self._from_training('num_epochs'))
    #
    # def hidden_layers(self):
    #     return [int(x) for x in self.get('NETWORK', 'hidden_layers').split(',')]
    #
    #
    # def training_path(self):
    #     return abs_path_of(self._from_paths('train_file'))
    #
    # def validation_path(self):
    #     return abs_path_of(self._from_paths('validation_file'))
    #
    # def targets(self):
    #     return [x for x in self.get('TARGETS', 'targets').split(',')]


def read_config(path):
    config = CustomConfigParser(inline_comment_prefixes=['#'], interpolation=configparser.ExtendedInterpolation())
    config.read(path)
    return config
