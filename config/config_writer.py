from configparser import ConfigParser
from pprint import pprint

from config.config_reader import CustomConfigParser
from config import config_reader
import os


class ConfigWriter:
    def __init__(self, path, name):
        self.path = path
        self.config = CustomConfigParser()
        self.create_config(name)

    def itemize(self, form):
        result = []
        for k, value in form.items():
            print(k, value)
            if 'csrf_token' not in k:
                section, key = k.split('-', 1)
                result.append((section.upper(), key, value))
        return result

    def populate_config(self, form):
        for section, key, value in self.itemize(form):
            self.add_item(section, key, value)

    def add_item(self, section, key, value):
        if section not in self.config.sections():
            self.config.add_section(section)
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
        self.add_item('INFO', 'config_name', name)
        self.add_item('INFO', 'config_path', self.path)
        if 'INFO' not in config_reader.read_config(self.path).sections():
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

    def load_input_output(self):
        self.load_section('input_paths_and_related_parameters')
        self.load_section('output_paths')
