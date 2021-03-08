import os

from flask import session, redirect, url_for

from config.config_writer import ConfigWriter
# from config import config_reader
import os

CONFIG_NAME = 'config.ini'


class Session:
    def __init__(self, app):
        self._config_writer = {}
        # self._config_reader = {}
        self._config = {}
        self._app = app
        self.log_path = 'output.txt'  # Later this is copied to the log path specified in the INI

    def read_log(self):
        return self.file_pointer.read()

    def open_log(self):
        if os.path.exists(self.log_path):
            os.remove(self.log_path)
        file = open(self.log_path, 'a')
        file.close()
        self.file_pointer = open(self.log_path)

    def add_user(self, user):
        self._config[user] = {}

    def reset_user(self):
        user = self.get_session('user')
        self._config[user] = {}

    def get_session(self, user_id):
        with self._app.app_context():
            if user_id not in session:
                return redirect(url_for('login'))
            return session[user_id]

    def create_config(self, path, name, dataset_name):
        user = self.get_session('user')
        config_path = os.path.join(path, CONFIG_NAME)
        self._config_writer[user] = ConfigWriter(config_path, name, dataset_name)
        # self._config_reader[user] = config_reader.read_config(config_path)

    def set_config(self, path, name, dataset_name):
        user = self.get_session('user')
        config_path = os.path.join(path, CONFIG_NAME)
        self._config_writer[user] = ConfigWriter(config_path, name, dataset_name)
        # self._config_reader[user] = config_reader.read_config(config_path)

    def get_writer(self):
        user = self.get_session('user')
        return self._config_writer[user]

    def get_reader(self):
        user = self.get_session('user')
        # return self._config_reader[user]