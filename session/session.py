import os

from flask import session, redirect, url_for

from config.config_writer import ConfigWriter
from config import config_reader

CONFIG_NAME = 'config.ini'


class Session:
    def __init__(self, app):
        self._config_writer = {}
        self._config_reader = {}
        self._config = {}
        self._app = app

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

    def create_config(self, path, name):
        user = self.get_session('user')
        config_path = os.path.join(path, CONFIG_NAME)
        self._config_writer[user] = ConfigWriter(config_path, name)
        self._config_reader[user] = config_reader.read_config(config_path)

    def set_config(self, path, name):
        user = self.get_session('user')
        config_path = os.path.join(path, CONFIG_NAME)
        self._config_writer[user] = ConfigWriter(config_path, name)
        self._config_reader[user] = config_reader.read_config(config_path)

    def get_writer(self):
        user = self.get_session('user')
        return self._config_writer[user]

    def get_reader(self):
        user = self.get_session('user')
        return self._config_reader[user]