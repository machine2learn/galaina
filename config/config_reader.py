import configparser
import os
from typing import Dict


INFO = 'info'


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
        return self._from_info('path')

    def get_name(self):
        return self._from_info('config_name')


    def path(self):
        return dict(self.items(PATHS))

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

    def features(self):
        return dict(self.items(FEATURES))

    def feature_slice(self):
        return self.get_as_slice(FEATURES, 'columns')

    def checkpoint_dir(self):
        return self.get_rel_path(PATHS, 'checkpoint_dir')

    def export_dir(self):
        return self.get_rel_path(PATHS, 'export_dir')

    def training_path(self):
        return abs_path_of(self._from_paths('train_file'))

    def validation_path(self):
        return abs_path_of(self._from_paths('validation_file'))

    def targets(self):
        return [x for x in self.get('TARGETS', 'targets').split(',')]

    # TODO TASK0?
    def label_slice(self):
        return self.get_as_slice(TASK0, 'ground_truth_column')

    def all(self):
        result = dict(self.items(TRAINING))
        result.update(self.items(EXPERIMENT))

        result.update(self.items(NETWORK))
        result.update(self.items(PATHS))
        result.update(self.items(CUSTOM_MODEL))

        int_columns = ["num_epochs", "batch_size", "validation_batch_size", "save_summary_steps",
                       "keep_checkpoint_max", "throttle", "validation_interval", "save_checkpoints_steps"]

        float_columns = ["learning_rate", "l1_regularization", "l2_regularization", "dropout"]


        for key in int_columns:
            if key in result:
                result[key] = int(result[key])

        for key in float_columns:
            if key in result:
                result[key] = float(result[key])

        result.update({'hidden_layers': self.hidden_layers()})
        result.update({'targets': self.targets()})

        return result


def read_config(path):
    config = CustomConfigParser(inline_comment_prefixes=['#'], interpolation=configparser.ExtendedInterpolation())
    config.read(path)
    return config


def get_task_sections(config):
    return {section_name: config[section_name] for section_name in config.sections() if
            section_name.startswith("TASK")}

# config = read_config("config/default.ini")
# print(config.get_slice("FEATURES","columns"))
# print ([1,2,3][config.get_slice("FEATURES","columns")])
