import os
import shutil


def delete_configs(config, dataset, username):
    if config != 'all':
        paths = [os.path.join('user_data', username, dataset, config)]
    else:
        paths = [os.path.join('user_data', username, dataset, d) for d in
                 os.listdir(os.path.join('user_data', username, dataset)) if
                 os.path.isdir(os.path.join('user_data', username, dataset, d))]
    for path in paths:
        shutil.rmtree(path)


def delete_dataset(APP_ROOT, username, dataset):
    path = os.path.join(APP_ROOT,'user_data', username, dataset)
    print('removing ...' + str(path))
    shutil.rmtree(path)


