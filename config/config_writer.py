from configparser import ConfigParser
from pprint import pprint

from config.config_reader import CustomConfigParser

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

    def create_config(self, name):
        self.add_item('INFO', 'config_name', name)
        self.add_item('INFO', 'config_path', self.path)
        self.write_config()

    def append_config(self):
        with open(self.path, 'a') as f:
            self.config.write(f)
