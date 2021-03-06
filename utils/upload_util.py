import os
from config import config_reader
from shutil import copyfile

from utils.profile import define_new_config_file


def existing_data(form, user_configs, username, sess, APP_ROOT):
    dataset_name = form['exisiting_files-train_file_exist']
    path = os.path.join(APP_ROOT, 'user_data', username, dataset_name)
    if form['exisiting_files-configuration'] != 'new_config':
        config_name = form['exisiting_files-configuration']
        sess.set('config_file', os.path.join(path, config_name, 'config.ini'))
        sess.load_config()
        return 'parameters'
    else:
        config_name = define_new_config_file(dataset_name, APP_ROOT, username)
        sess.set('config_file', os.path.join(path, config_name, 'config.ini'))
        if user_configs[dataset_name] and os.path.isfile(
                os.path.join(path, user_configs[dataset_name][0], 'config.ini')):
            reader = config_reader.read_config(os.path.join(path, user_configs[dataset_name][0], 'config.ini'))
            copyfile(os.path.join(path, user_configs[dataset_name][0], 'config.ini'),
                     os.path.join(path, config_name, 'config.ini'))
            filename = reader['PATHS']['file']
        elif os.path.isfile(os.path.join(path, dataset_name + '.csv')):
            filename = os.path.join(path, dataset_name + '.csv')
        else:
            filename = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and '.csv' in f][0]
        sess.set('file', os.path.join(path, filename))
        sess.get_writer().add_item('PATHS', 'file', os.path.join(path, filename))
        sess.get_writer().write_config(sess.get('config_file'))
        return 'slider'


