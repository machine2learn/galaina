import os


def create_user_path(username):
    if not os.path.exists('user_data'):
        os.mkdir('user_data')
    if not os.path.exists(os.path.join('user_data/', username)):
        os.mkdir(os.path.join('user_data/', username))