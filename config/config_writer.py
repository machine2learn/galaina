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
    def __init__(self, path, name):
        self.path = path
        self.config = CustomConfigParser()
        self.config.optionxform = str
        self.create_config(name)

    def itemize(self, form):
        result = []
        for k, value in form.items():
            print(k, value)
            if 'csrf_token' not in k:
                section, key = k.split('-', 1)
                section  = section.lower()
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
            if key == 'addelete':
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

    def create_config(self, name):
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
        self.add_item('output_paths', 'output_path_suffstat', '')
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
