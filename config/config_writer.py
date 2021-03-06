from configparser import ConfigParser
from pprint import pprint

from config.config_reader import CustomConfigParser
from config import config_reader
import os
import configparser

sections = ['info', 'copula_factor_algorithm', 'edge_weight_algorithm', 'pc_algorithm', 'plot_and_display',
            'input_paths_and_related_parameters', 'output_paths']
info = 'info'
copula_factor_algorithm = 'copula_factor_algorithm'
EDGE_WEIGHT = 'edge_weight_algorithm'
PC_ALGORITHM = 'pc_algorithm'
PLOT_AND_DISPLAY = 'plot_and_display'
input_paths_and_related_parameters = 'input_paths_and_related_parameters'
output_paths = 'output_paths'


class ConfigWriter:
    def __init__(self, path, name, dataset_name):
        self.path = path
        self.config = CustomConfigParser()
        self.config.optionxform = str
        self.create_config(name, dataset_name)

    def itemize(self, form):
        result = []
        for k, value in form.items():
            print(k, value)
            if 'csrf_token' not in k and not k.endswith('_check'):
                section, key = k.split('-', 1)
                section = section.lower()
                result.append((section, key, value))
        return result

    def keys(self):
        return self.config._sections.keys()

    def get_section(self, section):
        return self.config._sections[section]

    def get(self, section, param):
        return self.get_section(section)[param]

    def populate_config(self, form):
        for section, key, value in self.itemize(form):
            self.add_item(section, key, value)

    def add_item(self, section, key, value):
        if section not in self.config.sections():
            self.config.add_section(section)
        if section == 'pc_algorithm' or section == 'PC_ALGORITHM':
            if key == 'indeptest':
                key = 'indepTest'
            if key == 'numcores':
                key = 'numCores'
            if key == 'fixedgaps':
                key = 'fixedGaps'
            if key == 'fixededges':
                key = 'fixedEdges'
            if key == 'nadelete':
                key = 'NAdelete'
            if key == 'm_max':
                key = 'm.max'
            if key == 'm_max':
                key = 'm.max'
            if key == 'skel_method':
                key = 'skel.method'
            if key == 'maj_rule':
                key = 'maj.rule'
            if key == 'solve_confl':
                key = 'solve.confl'
        self.config.set(section, key, value)

    def write_config(self):
        with open(self.path, 'w') as f:
            self.config.write(f)

    def append_config(self):
        with open(self.path, 'a') as f:
            self.config.write(f)

    def replace_config(self):
        with open(self.path, 'wb') as f:
            self.config.write(f)

    def create_config(self, name, dataset_name):
        self.add_item('info', 'dataset_name', dataset_name)
        self.add_item('info', 'config_name', name)
        self.add_item('info', 'config_path', self.path)
        if 'info' not in config_reader.read_config(self.path).sections():
            self.write_config()

    def add_input_paths(self, dict_files):
        self.add_item('input_paths_and_related_parameters', 'input_data_ls', ','.join(dict_files['input_data_ls']))
        self.add_item('input_paths_and_related_parameters', 'input_factor_ls', ','.join(dict_files['input_factor_ls']))
        self.add_item('input_paths_and_related_parameters', 'ls_separator_str', ',')
        self.add_item('input_paths_and_related_parameters', 'column_separator_str', ',')
        self.add_item('input_paths_and_related_parameters', 'input_path_directed_edges_blacklist', '')
        self.write_config()

    def add_output_paths(self, path):
        self.add_item('output_paths', 'output_path_merged_data',
                      os.path.join(path, 'output', 'output_path_merged_data.csv'))
        self.add_item('output_paths', 'output_path_merged_factor_model_table',
                      os.path.join(path, 'output', 'output_path_merged_factor_model_table.csv'))
        self.add_item('output_paths', 'output_path_merged_factor_model_loading',
                      os.path.join(path, 'output', 'output_path_merged_factor_model_loading.csv'))
        self.add_item('output_paths', 'output_path_fig', os.path.join(path, 'output', 'output_path_fig.pdf'))
        self.add_item('output_paths', 'output_path_log', os.path.join(path, 'output', 'output_path_log.txt'))
        self.add_item('output_paths', 'output_path_suffstat', os.path.join(path, 'output', 'output_suff_stat.rds'))
        self.add_item('output_paths', 'output_path_pc_algo_obj',
                      os.path.join(path, 'output', 'output_path_pc_algo_obj.rds'))
        self.add_item('output_paths', 'output_path_bn_strength_obj',
                      os.path.join(path, 'output', 'output_path_bn_strength_obj.rds'))
        self.add_item('output_paths', 'output_path_avg_bn_obj',
                      os.path.join(path, 'output', 'output_path_avg_bn_obj.rds'))
        self.write_config()

    def load_section(self, section):
        reader = config_reader.read_config(self.path)
        if section in reader.sections():
            for key, val in reader[section].items():
                self.add_item(section, key, val)

    def load_all_sections(self):
        for sec in sections:
            self.load_section(sec)

    def add_r_front_end(self, APP_ROOT, path_r_binary_command):
        if 'r_front_end' not in self.keys():
            # self.add_item('r_front_end', 'path_r_binary_command', '/usr/local/bin/Rscript')
            self.add_item('r_front_end', 'path_r_binary_command', path_r_binary_command)
            self.add_item('r_front_end', 'r_binary_options', '--vanilla')
            self.add_item(
                'r_front_end', 'path_r_last_part_program',
                os.path.join(APP_ROOT, 'R_code', 'backend_galaina_final_part.R')
                # os.path.join(APP_ROOT, 'R_code', '20180725_use_config_ini_final_part.R')
            )
            self.add_item('r_front_end', 'path_r_infer_copula_factor_script',
                          os.path.join(APP_ROOT, 'R_code', 'my_inferCopulaFactorModel.R'))
            self.write_config()

    def get_info(self):
        path = self.get('info', 'config_path')
        dataset_name = self.get('info', 'dataset_name')
        return path, dataset_name


    def get_output_path_fig(self):
        return self.get('output_paths', 'output_path_fig')