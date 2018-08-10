import configparser
import os
from typing import Dict


INFO = 'INFO'


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


    def get_path(self):
        return self._from_info('config_path')

    def get_name(self):
        return self._from_info('config_name')

    def train_batch_size(self) -> int:
        return int(self._from_training('batch_size'))

    def learning_rate(self) -> int:
        return int(self._from_training('learning_rate'))

    def validation_batch_size(self) -> int:
        return int(self._from_training('validation_batch_size'))

    def optimizer(self) -> str:
        return self._from_training('optimizer')

    def l1_reqularization(self) -> float:
        return float(self._from_training('l1_regularization'))

    def l2_reqularization(self) -> float:
        return float(self._from_training('l2_regularization'))

    def num_epochs(self) -> int:
        return int(self._from_training('num_epochs'))

    def hidden_layers(self):
        return [int(x) for x in self.get('NETWORK', 'hidden_layers').split(',')]


    def training_path(self):
        return abs_path_of(self._from_paths('train_file'))

    def validation_path(self):
        return abs_path_of(self._from_paths('validation_file'))

    def targets(self):
        return [x for x in self.get('TARGETS', 'targets').split(',')]


def read_config(path):
    config = CustomConfigParser(inline_comment_prefixes=['#'], interpolation=configparser.ExtendedInterpolation())
    config.read(path)
    return config



