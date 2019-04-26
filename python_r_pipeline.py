#!/usr/bin/env python3
# coding: utf-8

import sys
import functools
import re
import subprocess
import configparser

import pandas as pd


class MyConfigReader:
    def __init__(self, configparser_obj):
        self.configparser_obj = configparser_obj

    def get_input_factor_ls(self):
        input_factor_ls_str = self.configparser_obj.get(
            section="input_paths_and_related_parameters",
            option="input_factor_ls"
        )
        return self.get_list(input_factor_ls_str)

    def get_input_data_ls(self):
        input_data_ls_str = self.configparser_obj.get(
            section="input_paths_and_related_parameters",
            option="input_data_ls"
        )
        return self.get_list(input_data_ls_str)

    def get_input_path_directed_edges_blacklist(self):
        return self.configparser_obj.get(
            section="input_paths_and_related_parameters",
            option="input_path_directed_edges_blacklist",
            # fallback=None
        )

    def get_output_path_merged_data(self):
        return self.configparser_obj.get(
            section="output_paths",
            option="output_path_merged_data"
        )

    def get_output_path_merged_factor_model_loading(self):
        return self.configparser_obj.get(
            section="output_paths",
            option="output_path_merged_factor_model_loading"
        )

    def get_output_path_merged_factor_model_table(self):
        return self.configparser_obj.get(
            section="output_paths",
            option="output_path_merged_factor_model_table"
        )

    def get_csv_separator(self):
        return self.configparser_obj.get(
            section="input_paths_and_related_parameters",
            option="column_separator_str"
        )

    def get_ls_separator(self):
        return self.configparser_obj.get(
            section="input_paths_and_related_parameters",
            option="ls_separator_str"
        )

    def get_list(self, input_str):
        return input_str.split(sep=self.get_ls_separator())

    def get_path_to_config_ini(self):
        return self.configparser_obj.get(
            section="added_section_01", option="path_to_config_ini"
        )

    def get_r_binary_command(self):
        return self.configparser_obj.get(
            section="r_front_end", option="path_r_binary_command"
        )

    def get_r_binary_options(self):
        return self.configparser_obj.get(
            section="r_front_end", option="r_binary_options"
        )

    def get_path_r_last_part_program(self):
        return self.configparser_obj.get(
            section="r_front_end", option="path_r_last_part_program"
        )


class FactorInfoFactory:
    def __init__(self, factor_df):
        self.factor_df = factor_df

    def factory(self):
        factor_type_str = self.determine_factor_type()
        # instead of if, use dictionary of function names
        type2obj = {
            'factor_model_df': FactorModelDF,
            'factor_loading_df': FactorLoadingDF,
            'factor_var_set_df': FactorVarSetDF
        }
        return type2obj[factor_type_str](self.factor_df)

    def determine_factor_type(self):
        """
        Check the index name to determine factor types.
        3 possible factor types:
        - Factor Model (used also internally to merge the multiple sources):
            index.name = 'Variable'
            columns = ['Factor', 'Loading']
        - Factor Loading Matrix (must be used to export data for R script):
            index.name = None - for compatibility with R
            columns = All factors (not checked now),
            loading_df.loc[variable_x, factor_y] = info2fd['factor_model_df'].loc[variable_x, 'Loading'] iff factor_y is latent factor of variable_x
        - Factor Variable Set
            index.name = 'Factor'
            column = 'Variable_set'
            df.loc[factor, variable_set] = formula expressing factor as a function of the variable in variable_set
        """
        column_tuple2detected_format_str = {
            ('Factor', 'Loading'): 'factor_model_df',
            ('Variable_set',): 'factor_var_set_df',
        }
        tmp_tuple = tuple(self.factor_df.columns.tolist())
        # if no compatible key is in the dictionary, return factor_loading_df
        factor_type_str = column_tuple2detected_format_str.get(tmp_tuple, 'factor_loading_df')
        # factor_type_str = 'factor_loading_df'  # assume it is a Loading Matrix if index has no name
        # print(factor_df.index.name)
        detected_format_str2index_name = {
            'factor_model_df': 'Variable',
            'factor_loading_df': None,
            'factor_var_set_df': 'Factor'
        }
        assert self.factor_df.index.name == detected_format_str2index_name[factor_type_str], \
            'Factor DataFrame index name incompatible with columns'
        return factor_type_str


class FactorInfo:
    def __init__(self, factor_df):
        # TODO maybe define factor_df as a property also later, see below. So we fix issues with modifying the index name
        # https://stackoverflow.com/questions/7374748/whats-the-difference-between-a-python-property-and-attribute
        # TODO use FactorInfoFactory
        tmp_df = factor_df.copy()
        tmp_df.columns.rename(None, inplace=True)
        # factor_type_str = FactorInfoFactory(factor_df).determine_factor_type()
        self.factor_df = tmp_df
        # self.factor_df.columns.name = None
        # self.factor_df.columns.rename(None, inplace=True)

    # Otherwise it is = info2df['factor_var_set_df'].index.name
    # TODO actually here we could put the adjustment of columns.name and index.nme
    # self.factor_type_str


class FactorModelDF(FactorInfo):
    def __init__(self, factor_df):
        super().__init__(factor_df)

    def adjust_factor_model_index_and_columns_name(
            self):  # TODO use inheritance and replace with adjust_factor_index_and_columns_names
        tmp_df = self.factor_df.copy()
        tmp_df.index.rename('Variable', inplace=True)
        self.factor_df = tmp_df

    def from_factor_model_to_factor_loading(self):
        """
        Index = variables
        Columns = all factors
        loading_df.loc[variable_x, factor_y] = info2fd['factor_model_df'].loc[variable_x, 'Loading'] iff factor_y is latent factor of variable_x
        :return: FactorLoadingDF()
        """
        # Create loading matrix
        loading_df = self.factor_df.set_index([self.factor_df.index, 'Factor']).unstack(
            level=1).fillna(0.0).astype(float)
        loading_df.columns = loading_df.columns.droplevel(0)  # Remove the column header 'Loading'
        factor_loading_obj = FactorLoadingDF(loading_df)
        factor_loading_obj.adjust_factor_loading_index_and_columns_name()  # TODO inheritance
        return factor_loading_obj

    def get_unstructured_and_structured_variables(self):
        """
        unstructured variables = variables that are also factor of its own, and the only factor
        Used for merged_info2df
        """
        unstructured_variable_index = self.factor_df.index.intersection(self.factor_df['Factor'])
        structured_variable_index = self.factor_df.index.drop(unstructured_variable_index)
        return unstructured_variable_index, structured_variable_index


class FactorLoadingDF(FactorInfo):
    def __init__(self, factor_df):
        super().__init__(factor_df)

    def adjust_factor_loading_index_and_columns_name(self):  # TODO use inheritance
        tmp_df = self.factor_df.copy()
        tmp_df.index.rename(None, inplace=True)
        self.factor_df = tmp_df

    def check_factor_loading(self):
        factor_bdf = (self.factor_df != 0)
        assert (factor_bdf.sum(axis=1) == 1).all(), \
            "Zero or more than one factor per variable, factor models are not pure"
        used_factor_bse = (factor_bdf.sum(axis=0) >= 1)
        if not used_factor_bse.all():
            # log warning
            print("Following factors are not used:")
            print(*used_factor_bse[~used_factor_bse].index.values, sep="\n")

    def from_loading_to_factor_model(self):
        """
        From a factor loading DataFrame, make the Factor model Dataframe
        """
        self.check_factor_loading()
        factor_model_df = self.factor_df.copy()
        # TODO rewrite it without loops?
        factor_model_df['Factor'] = (factor_model_df != 0).idxmax(axis=1)  # Make Factor column
        for iFactor, iDF in factor_model_df.groupby('Factor'):  # Make Loading column
            factor_model_df.loc[iDF.index, 'Loading'] = iDF[iFactor]
        factor_model_obj = FactorModelDF(factor_model_df[['Factor', 'Loading']])
        factor_model_obj.adjust_factor_model_index_and_columns_name()  # TODO use inheritance
        return factor_model_obj


class FactorVarSetDF(FactorInfo):
    def __init__(self, factor_df):
        super().__init__(factor_df)

    def from_var_set_to_factor_loading(self):
        """
        From expression "a1 * x1 + a2 * x2 +..." with index "f1" create row of loading matrix s.t.
        loading_matrix_df.loc["f1", "x1"] = a1 ** -1. Variables absent in the formula will be set to zero
        """
        pattern_str = r"([+-])?\s*(?:(\d+)\s*\*\s*)?([a-zA-Z]\w*)"
        loading_df = self.factor_df['Variable_set'].apply(
            lambda formula_str: pd.Series(
                {iVar: int(iSign + (iCoeff or '1')) ** -1
                 for (iSign, iCoeff, iVar) in re.findall(pattern_str, formula_str)
                 }
            )
        ).fillna(0.0).astype(float).T
        factor_loading_df = FactorLoadingDF(loading_df)
        factor_loading_df.adjust_factor_loading_index_and_columns_name()

        return factor_loading_df


class InfoToDF:
    def __init__(self, name, data_df, factor_df):
        self.name = name
        self.data_df = data_df
        info_factory_obj = FactorInfoFactory(factor_df)
        factor_type_str = info_factory_obj.determine_factor_type()  # actually used inside the factory() method
        factor_obj = info_factory_obj.factory()
        self.type2factor = {factor_type_str: factor_obj}

    def convert_to_factor_model_type_if_needed(self):
        if 'factor_model_df' not in self.type2factor.keys() and 'factor_loading_df' in self.type2factor.keys():
            self.type2factor['factor_model_df'] = self.type2factor['factor_loading_df'].from_loading_to_factor_model()
        assert 'factor_model_df' in self.type2factor.keys(), "Conversion to FactorModel object failed"
        return self

    def check_factor_model(self):
        # TODO only works with factorInfo of FactorModelDF type
        assert set(self.type2factor['factor_model_df'].factor_df.index) == set(self.data_df.columns), \
            "Some variables have no factor assigned: \n Variables for Factors but not in Data:{}\n Variables for Data but not in Factors:{}".format(
                set(self.type2factor['factor_model_df'].factor_df.index).difference(self.data_df.columns),
                set(self.data_df.columns).difference(self.type2factor['factor_model_df'].factor_df.index)
            )  # check match data and factor
        assert self.type2factor['factor_model_df'].factor_df['Factor'].notnull().all().all(), \
            "Some variables have null factor assigned"
        assert (self.type2factor['factor_model_df'].factor_df['Loading'] != 0).all().all(), \
            "Some variable-factor combinations have zero or null loading coefficient"

    def sort_factors_and_variables_with_unstructured_first(self):
        """
        Just doing it because it looks to me that Ruifei code indirectly requires it
        Used for merged_info2df
        """
        # TODO only works with factorInfo of FactorModelDF type
        unstructured_variable_index, structured_variable_index = \
            self.type2factor['factor_model_df'].get_unstructured_and_structured_variables()
        unstructured_variable_first_index = unstructured_variable_index.append(structured_variable_index)
        self.type2factor['factor_model_df'].factor_df = self.type2factor['factor_model_df'].factor_df.loc[
            unstructured_variable_first_index].copy()  # put factor-unstructured variables first
        # not needed - columns of data_df sorted as rows of factor_model_df
        self.data_df = self.data_df.loc[:, unstructured_variable_first_index].copy()

    def check_background_knowledge_and_merged_info(self, config):
        """
        Used for merged_info2df
        :param config:
        :return:
        """
        tmp_path_to_bgk = config.get_input_path_directed_edges_blacklist()
        if tmp_path_to_bgk:  # != "":  #is not None:
            bgk_df = pd.read_csv(tmp_path_to_bgk, sep=config.get_csv_separator())
            assert set(bgk_df.values.flatten()).issubset(self.type2factor['factor_loading_df'].factor_df.columns), \
                "Some nodes of blacklist directed edges are neither factors nor unstructured variables"

    # TODO now added to InfoToDF, but we use only after merging
    def check_merged_all_factor_model_from_list(self):  # merged_info2df, ls_of_info2df=None):
        # TODO check that variables and factors are unique per data type
        # PROBLEM requires both a FactorModelDF and a FactorLoadingDF object
        #  After concatenation, check that there are
        # 1. no duplicates among factors of different source (difficult to check)
        #     assert (merged_info2df['Factor'].duplicated() == False).all(), "Duplicated Factor"
        #     if ls_of_info2df:
        #         name2factor_ls = {}
        #         for iInfo2df in ls_of_info2df:
        #             if iInfo2df['name'] == 'unstructured':
        #                 unstructured_variable_index = iInfo2df['factor_model_df'].index
        #             name2factor_ls[iInfo2df['name']] = set(iInfo2df['factor_model_df']['Factor'].tolist())
        #         for iNameSetTuple1, iNameSetTuple2 in itertools.combinations(name2factor_ls.items(), 2):
        #             assert bool(iNameSetTuple1[1] & iNameSetTuple2[1]) == False, "Some Factors are shared among {} and {}".format(iNameSetTuple1[0], iNameSetTuple2[0])
        #     else:

        # unstructured_variable_index, structured_variable_index = get_unstructured_and_structured_variables(merged_info2df)
        unstructured_variable_index, structured_variable_index = self.type2factor[
            'factor_model_df'].get_unstructured_and_structured_variables()

        assert (self.type2factor['factor_model_df'].factor_df.loc[
                    self.type2factor['factor_model_df'].factor_df['Factor'].isin(
                        unstructured_variable_index), 'Factor'].tolist()
                == unstructured_variable_index.tolist()), 'Unstructured variable are also factors of structured variables'

        assert (self.type2factor['factor_model_df'].factor_df.loc[unstructured_variable_index, 'Loading'] == 1).all(), \
            'Some unstructured variables do not have Loading = 1'
        assert (self.type2factor['factor_model_df'].factor_df.loc[unstructured_variable_index, 'Factor'].tolist()
                == unstructured_variable_index.tolist()), "Some unstructured variables are not factors of themselves"

        #         assert merged_info2df['factor_loading_df'].factor_df.loc[unstructured_variable_index, unstructured_variable_index]

        assert self.type2factor['factor_loading_df'].factor_df.columns.is_unique, "Duplicated factors"
        # 2. no duplicates among variables
        assert self.type2factor['factor_model_df'].factor_df.index.is_unique, "Duplicated variables"
        # 3. no duplicates among factors

        #     assert (merged_info2df.index.duplicated() == False).all(), "No duplicate Variables"
        # 3. intersection between variables and factors == unstructured variables
        # TODO check that this concat is append
        #     structured_variable_index = merged_info2df['factor_model_df'].index.drop(unstructured_variable_index)
        assert not structured_variable_index.isin(self.type2factor['factor_model_df'].factor_df['Factor']).any(), \
            "Some factors are also variables"
        # TODO TEST maybe check that all proper factor have at least 2 associated variables
        assert (self.type2factor['factor_model_df'].factor_df.loc[
                    structured_variable_index, 'Factor'].value_counts() > 1).all(), \
            "Some factors of structured variables have just one structured variable associated"


class ListInfoToDF:
    def __init__(self, ls_of_info2df):
        self.ls_info2df = ls_of_info2df

    def prepare_for_merging(self):
        self.ls_info2df = [iInfo2df.convert_to_factor_model_type_if_needed() for iInfo2df in self.ls_info2df]

    def check_before_merging(self):
        for iInfo2df in self.ls_info2df:  # Check factor model before merging
            iInfo2df.check_factor_model()

    def merge_all_info2df(self):
        """
        Merge factor models, then put unstructured variables as first columns of factor_model_df.
        Now sort also rows s.t. unstructured var are on top, and then use this variable sorting for data_df columns
        :return: InfoToDF obj
        """

        tmp_info2df = InfoToDF(
            name='intersection',
            data_df=self.merge_all_data(),
            factor_df=self.merge_all_factor_model()
        )
        tmp_info2df.sort_factors_and_variables_with_unstructured_first()
        return tmp_info2df

    def merge_all_data(self):
        """
        TODO decide if set the index.name to None here or just use the saving
        :return: pd.DataFrame()
        """
        if len(self.ls_info2df) == 1:
            tmp_df = self.ls_info2df[0].data_df.copy()
        else:
            tmp_df = functools.reduce(
                lambda left, right: left.join(right, how='inner'), [iInfo2df.data_df for iInfo2df in self.ls_info2df]
            )
        return tmp_df

    def merge_all_factor_model(self):
        return pd.concat([iInfo2df.type2factor['factor_model_df'].factor_df for iInfo2df in self.ls_info2df])


def save_to_csv_for_r(input_df, path_to_data_file, config):
    """
    :param input_df:
    :param path_to_data_file:
    :param config:
    :return:
    """
    input_df.to_csv(path_to_data_file, sep=config.get_csv_separator(), index_label=False)


def load_input_csv(path_to_data_file, config):
    return pd.read_csv(path_to_data_file, sep=config.get_csv_separator(), index_col=0)


def load_data(path_to_data_file, config):
    return load_input_csv(path_to_data_file, config)


def load_factor(path_to_factor_file, config):
    return load_input_csv(path_to_factor_file, config)


def check_data_factor_number_match(config):
    input_data_ls = config.get_input_data_ls()
    input_factor_ls = config.get_input_factor_ls()
    assert len(input_data_ls) == len(
        input_factor_ls), "No one-to-one match between data and factor models"  # is a check


def create_ls_of_info2df(config):
    input_data_ls = config.get_input_data_ls()
    input_factor_ls = config.get_input_factor_ls()
    ls_of_info2df = []
    for iPos_n in range(len(input_data_ls)):
        # Maybe we
        print('Process {}'.format(input_data_ls[iPos_n]))
        tmp_input_data_df = load_data(input_data_ls[iPos_n], config)
        print('Process {}'.format(input_factor_ls[iPos_n]))
        tmp_input_factor_df = load_factor(input_factor_ls[iPos_n], config)
        tmp_info2df_obj = InfoToDF(
            name='source_{}'.format(iPos_n),
            data_df=tmp_input_data_df,
            factor_df=tmp_input_factor_df
        )
        if 'factor_var_set_df' in tmp_info2df_obj.type2factor.keys():
            tmp_info2df_obj.type2factor['factor_loading_df'] = tmp_info2df_obj.type2factor[
                'factor_var_set_df'].from_var_set_to_factor_loading()

        ls_of_info2df += [tmp_info2df_obj]
    return ls_of_info2df


def main(configparser_obj):
    config = MyConfigReader(configparser_obj)
    # , Load input DataFrames
    check_data_factor_number_match(config)
    ls_of_info2df = create_ls_of_info2df(config)

    # . Merge them all
    ls_info_to_df_obj = ListInfoToDF(ls_of_info2df)
    ls_info_to_df_obj.prepare_for_merging()
    ls_info_to_df_obj.check_before_merging()
    merged_info2df = ls_info_to_df_obj.merge_all_info2df()

    # . Create the loading matrix of the merged data
    merged_info2df.type2factor['factor_loading_df'] = \
        merged_info2df.type2factor['factor_model_df'].from_factor_model_to_factor_loading()
    # . Additional checks
    merged_info2df.check_merged_all_factor_model_from_list()
    # check_background_knowledge_and_merged_info(merged_info2df, config)

    save_to_csv_for_r(
        merged_info2df.data_df,
        config.get_output_path_merged_data(),
        config
    )

    merged_info2df.type2factor['factor_loading_df'].factor_df.to_csv(
        config.get_output_path_merged_factor_model_loading(), sep=config.get_csv_separator()
    )

    # Check if write index.name to CSV for FactorModelDF
    merged_info2df.type2factor['factor_model_df'].factor_df.to_csv(
        config.get_output_path_merged_factor_model_table(), sep=config.get_csv_separator(),
    )

    # Continue in R
    subprocess.call(
        [
            config.get_r_binary_command(),
            config.get_r_binary_options(),
            config.get_path_r_last_part_program(),
            config.get_path_to_config_ini()
        ],
        shell=False
    )


def run(cfg, path):
    # Create section with config INI file name and path to R script to be run
    if 'added_section_01' not in cfg.sections():
        cfg.add_section("added_section_01")
        cfg.set(section="added_section_01", option="path_to_config_ini", value=path)
    main(cfg)


if __name__ == "__main__":
    cfg = configparser.ConfigParser()
    cfg.read(sys.argv[1])
    run(cfg, sys.argv[1])
