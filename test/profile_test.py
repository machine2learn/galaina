
import pytest
from _pytest import monkeypatch

import utils
from utils.profile import generate_dataset_name


def test_generate_dataset_name_dataset_only_one(monkeypatch):

    def mockreturn(app_root, username):
        return ['dataset_1']

    def path_already_exists(app_root, username):
        return True

    monkeypatch.setattr(utils.profile, 'find_all_datasets', mockreturn)
    monkeypatch.setattr(utils.profile, 'path_already_exists', path_already_exists)
    result = generate_dataset_name('a', 'b', 'dataset')
    assert result == 'dataset_2'

def test_generate_dataset_name_if_empty(monkeypatch):

    def mockreturn(app_root, username):
        return []

    def path_already_exists(app_root, username):
        return True

    monkeypatch.setattr(utils.profile, 'find_all_datasets', mockreturn)
    monkeypatch.setattr(utils.profile, 'path_already_exists', path_already_exists)
    result = generate_dataset_name('a', 'b', 'dataset')
    assert result == 'dataset_1'

def test_generate_dataset_name_non_existing_name(monkeypatch):

    def mockreturn(app_root, username):
        return ['dataset_1']

    def path_already_exists(app_root, username):
        return True

    monkeypatch.setattr(utils.profile, 'find_all_datasets', mockreturn)
    monkeypatch.setattr(utils.profile, 'path_already_exists', path_already_exists)
    result = generate_dataset_name('a', 'b', 'userdataset')
    assert result == 'userdataset'

def test_generate_dataset_name_existing_dataset_name(monkeypatch):

    def mockreturn(app_root, username):
        return ['userdataset']

    def path_already_exists(app_root, username):
        return True

    monkeypatch.setattr(utils.profile, 'find_all_datasets', mockreturn)
    monkeypatch.setattr(utils.profile, 'path_already_exists', path_already_exists)
    result = generate_dataset_name('a', 'b', 'userdataset')
    assert result == 'userdataset_1'

def test_generate_dataset_name_existing_dataset_name_with_one_index(monkeypatch):

    def mockreturn(app_root, username):
        return ['userdataset', 'userdataset_1']

    def path_already_exists(app_root, username):
        return True

    monkeypatch.setattr(utils.profile, 'find_all_datasets', mockreturn)
    monkeypatch.setattr(utils.profile, 'path_already_exists', path_already_exists)
    result = generate_dataset_name('a', 'b', 'userdataset')
    assert result == 'userdataset_2'
