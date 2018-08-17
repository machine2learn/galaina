from config.config_reader import read_config


def check_causal_discovery_ob(form):
    if form['pc_algorithm-causal_discovery_observation_n'].lower() == 'automatic':
        form['pc_algorithm-causal_discovery_observation_n'] = '0'


def set_form(form, config):
    config = read_config(config.get('info', 'config_path'))
    if 'copula_factor_algorithm' in config.keys():
        form.copula_factor_algorithm.form.gibbs_sampling_n.default = config.get('copula_factor_algorithm',
                                                                                 'gibbs_sampling_n')
        form.copula_factor_algorithm.form.gibbs_burn_in_n.default = config.get('copula_factor_algorithm',
                                                                               'gibbs_burn_in_n')
        form.copula_factor_algorithm.form.gibbs_first_random_seed_n.default = config.get('copula_factor_algorithm',
                                                                                         'gibbs_first_random_seed_n')
        form.copula_factor_algorithm.form.gibbs_random_seed_update_parameter_n.default = config.get(
            'copula_factor_algorithm',
            'gibbs_random_seed_update_parameter_n')
        form.copula_factor_algorithm.form.process()
    if 'edge_weight_algorithm' in config.keys():
        form.edge_weight_algorithm.form.bootstrap_n.default = config.get('edge_weight_algorithm', 'bootstrap_n')
        form.edge_weight_algorithm.form.bootstrap_first_random_seed_n.default = config.get('edge_weight_algorithm',
                                                                                           'bootstrap_first_random_seed_n')
        form.edge_weight_algorithm.form.bootstrap_random_seed_update_parameter_n.default = config.get(
            'edge_weight_algorithm',
            'bootstrap_random_seed_update_parameter_n')
        form.edge_weight_algorithm.form.process()
    if 'pc_algorithm' in config.keys():
        form.pc_algorithm.form.causal_discovery_observation_n.default = config.get('pc_algorithm',
                                                                                   'causal_discovery_observation_n')
        form.pc_algorithm.form.indepTest.default = config.get('pc_algorithm', 'indepTest')
        form.pc_algorithm.form.alpha.default = config.get('pc_algorithm', 'alpha')
        form.pc_algorithm.form.numCores.default = config.get('pc_algorithm', 'numCores')
        form.pc_algorithm.form.verbose.default = config.get('pc_algorithm', 'verbose')
        form.pc_algorithm.form.NAdelete.default = config.get('pc_algorithm', 'NAdelete')
        form.pc_algorithm.form.m_max.default = config.get('pc_algorithm', 'm.max')
        form.pc_algorithm.form.u2pd.default = config.get('pc_algorithm', 'u2pd')
        form.pc_algorithm.form.skel_method.default = config.get('pc_algorithm', 'skel.method')
        form.pc_algorithm.form.conservative.default = config.get('pc_algorithm', 'conservative')
        form.pc_algorithm.form.maj_rule.default = config.get('pc_algorithm', 'maj.rule')
        form.pc_algorithm.form.solve_confl.default = config.get('pc_algorithm', 'solve.confl')
        form.pc_algorithm.form.fixedGaps.default = config.get('pc_algorithm', 'fixedGaps')
        form.pc_algorithm.form.fixedEdges.default = config.get('pc_algorithm', 'fixedEdges')
        form.pc_algorithm.form.process()
    if 'plot_and_display' in config.keys():
        form.plot_and_display.form.core_plot_title_str.default = config.get('plot_and_display', 'core_plot_title_str')
        form.plot_and_display.form.process()

# def set_form(form, reader):
#     if 'COPULA_FACTOR' in reader.keys():
#         form.copula_factor.form.gibbs_sampling_n.default = reader.get_gibbs_sampling_n()
#         form.copula_factor.form.gibbs_burn_in_n.default = reader.get_gibbs_burn_in_n()
#         form.copula_factor.form.gibbs_first_random_seed_n.default = reader.get_gibbs_first_random_seed_n()
#         form.copula_factor.form.gibbs_random_seed_update_parameter_n.default = reader.get_gibbs_random_seed_update_parameter_n()
#         form.copula_factor.form.process()
#     if 'EDGE_WEIGHT' in reader.keys():
#         form.edge_weight.form.bootstrap_n.default = reader.get_bootstrap_n()
#         form.edge_weight.form.bootstrap_first_random_seed_n.default = reader.get_bootstrap_first_random_seed_n()
#         form.edge_weight.form.bootstrap_random_seed_update_parameter_n.default = reader.get_bootstrap_random_seed_update_parameter_n()
#         form.edge_weight.form.process()
#     if 'PC_ALGORITHM' in reader.keys():
#         form.pc_algorithm.form.causal_discovery_observation_n.default = reader.get_causal_discovery_observation_n()
#         form.pc_algorithm.form.indepTest.default = reader.get_indeptest()
#         form.pc_algorithm.form.alpha.default = reader.get_alpha()
#         form.pc_algorithm.form.numCores.default = reader.get_numcores()
#         form.pc_algorithm.form.verbose.default = reader.get_verbose()
#         form.pc_algorithm.form.NAdelete.default = reader.get_nadelete()
#         form.pc_algorithm.form.m_max.default = reader.get_m_max()
#         form.pc_algorithm.form.u2pd.default = reader.get_u2pd()
#         form.pc_algorithm.form.skel_method.default = reader.get_skel_method()
#         form.pc_algorithm.form.conservative.default = reader.get_conservative()
#         form.pc_algorithm.form.maj_rule.default = reader.get_maj_rule()
#         form.pc_algorithm.form.solve_confl.default = reader.get_solve_confl()
#         form.pc_algorithm.form.fixedGaps.default = reader.get_fixedgaps()
#         form.pc_algorithm.form.fixedEdges.default = reader.get_fixededges()
#         form.pc_algorithm.form.process()
#     if 'PLOT_AND_DISPLAY' in reader.keys():
#         form.plot_and_display.form.core_plot_title_str.default = reader.get_core_plot_title_str()
#         form.plot_and_display.form.process()
