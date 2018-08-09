import configparser
import os

from werkzeug.utils import secure_filename

from config.config_writer import ConfigWriter

USER_DATA = 'user_data'
ALLOWED_EXTENSIONS = ['csv']


def generate_dataset_name(app_root, username, datasetname):
    datasetname = datasetname or 'dataset'
    user_datasets = []
    if os.path.isdir(os.path.join(app_root, 'user_data', username)):
        user_datasets = [a for a in os.listdir(os.path.join(app_root, 'user_data', username))
                         if os.path.isdir(os.path.join(app_root, 'user_data', username, a))]

    latest = 1
    if datasetname in user_datasets:
        datasets = [conf_name.split(datasetname + '_')[1] for conf_name in user_datasets]
        if datasets:
            latest = max(datasets) + 1

    return datasetname + '_' + str(latest)


def get_configs_files(app_root, username):
    user_configs = {}
    existing_datasets = []
    configs = {}
    user_path = os.path.join(app_root, USER_DATA, username)
    user_datasets = [dataset for dataset in os.listdir(user_path) if os.path.isdir(os.path.join(user_path, dataset))]
    for user_dataset in user_datasets:
        user_configs[user_dataset] = [config for config in os.listdir(os.path.join(user_path, user_dataset)) if
                                      os.path.isdir(os.path.join(user_path, user_dataset, config))]

        for config_file in user_configs[user_dataset]:
            dataset_config = user_dataset + '_' + config_file
            config = configparser.ConfigParser()
            config.read(os.path.join(user_path, user_dataset, config_file, 'config.ini'))
            configs[dataset_config] = {'name': config.get('INFO', 'config_name')}

        existing_datasets.append((user_dataset, user_dataset))
    return existing_datasets, user_configs, configs


def copyfile(src, dst):
    from shutil import copyfile
    if os.path.exists(src): copyfile(src, dst)


def create_config(username, APP_ROOT, dataset, config_name, sess):
    path = os.path.join(APP_ROOT, USER_DATA, username, dataset, config_name)
    os.makedirs(path, exist_ok=True)
    sess.create_config(path, config_name)


def generate_config_name(app_root, username, dataset_name):
    user_configs = []
    if os.path.isdir(os.path.join(app_root, USER_DATA, username, dataset_name)):
        user_configs = [a for a in os.listdir(os.path.join(app_root, USER_DATA, username, dataset_name))
                        if os.path.isdir(os.path.join(app_root, USER_DATA, username, dataset_name, a))]
    dataset_configs = [int(conf_name.rsplit('_')[1]) for conf_name in user_configs]
    latest = 1 if not dataset_configs else max(dataset_configs) + 1
    return 'config_' + str(latest)


def define_new_config_file(dataset_name, APP_ROOT, username, sess):
    config_name = generate_config_name(APP_ROOT, username, dataset_name)
    target = os.path.join(APP_ROOT, USER_DATA, username, dataset_name, config_name)
    if not os.path.isdir(target):
        os.makedirs(target, exist_ok=True)
    create_config(username, APP_ROOT, dataset_name, config_name, sess)
    return config_name


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


def save_files(target, files):
    for filename, file in files.items():
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            destination = os.path.join(target, filename)
            file.save(destination)
    return True


def new_dataset(dataset_name, filenames, APP_ROOT, username, sess):
    dataset_name = generate_dataset_name(APP_ROOT, username, dataset_name)
    destination = os.path.join(APP_ROOT, USER_DATA, username, dataset_name)

    if not os.path.isdir(destination):
        os.makedirs(destination, exist_ok=True)

    save_files(destination, filenames)
    config_name = define_new_config_file(dataset_name, APP_ROOT, username, sess)

    return 'slider'
