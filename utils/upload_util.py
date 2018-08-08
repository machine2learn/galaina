import os
from config import config_reader
from utils import config_ops
from shutil import copyfile


def existing_data(form, user_configs, username, sess, APP_ROOT):
    dataset_name = form['exisiting_files-train_file_exist']
    path = os.path.join(APP_ROOT, 'user_data', username, dataset_name)
    if form['exisiting_files-configuration'] != 'new_config':
        config_name = form['exisiting_files-configuration']
        sess.set('config_file', os.path.join(path, config_name, 'config.ini'))
        sess.load_config()
        return 'parameters'
    else:
        config_name = config_ops.define_new_config_file(dataset_name, APP_ROOT, username, sess.get_writer())
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


def generate_dataset_name(app_root, username, dataset_name):
    user_datasets = []
    if os.path.isdir(os.path.join(app_root, 'user_data', username)):
        user_datasets = [a for a in os.listdir(os.path.join(app_root, 'user_data', username))
                         if os.path.isdir(os.path.join(app_root, 'user_data', username, a))]
    cont = 1
    while dataset_name + '_' + str(cont) in user_datasets:
        cont += 1
    new_dataset_name = dataset_name + '_' + str(cont)
    return new_dataset_name

