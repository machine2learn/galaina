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


def delete_dataset(username, dataset):
    path = os.path.join('user_data', username, dataset,'input')
    shutil.rmtree(path)
    path = os.path.join('user_data', username, dataset,'factor')
    shutil.rmtree(path)

def delete_entire_dataset(username, dataset):
    path = os.path.join('user_data', username, dataset)
    print('removing ...' + str(path))
    shutil.rmtree(path)
