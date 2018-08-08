from flask import session, redirect, url_for

from config.config_writer import ConfigWriter


class Session:
    def __init__(self, app):
        self._config_writer = {}
        self._config = {}
        self._app = app

    def add_user(self, user):
        self._config[user] = {}
        self._config_writer[user] = ConfigWriter()

    def reset_user(self):
        user = self.get_session('user')
        self._config[user] = {}
        self._config_writer[user] = ConfigWriter()

    def get_session(self, user_id):
        with self._app.app_context():
            if user_id not in session:
                return redirect(url_for('login'))
            return session[user_id]
