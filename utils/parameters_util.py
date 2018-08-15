def check_causal_discovery_ob(form):
    if form['pc_algorithm-causal_discovery_observation_n'].lower() == 'automatic':
        form['pc_algorithm-causal_discovery_observation_n'] = '0'


def set_form(form, reader):
    if 'COPULA_FACTOR' in reader.keys():
        form.copula_factor.form.gibbs_sampling_n.default = reader.get_gibbs_sampling_n()
        form.copula_factor.form.gibbs_burn_in_n = reader.get_gibbs_burn_in_n()
        form.copula_factor.form.gibbs_first_random_seed_n = reader.get_gibbs_first_random_seed_n()
        form.copula_factor.form.gibbs_random_seed_update_parameter_n = reader.get_gibbs_random_seed_update_parameter_n()
        form.copula_factor.form.process()
    if 'EDGE_WEIGHT' in reader.keys():
        form.edge_weight.form.bootstrap_n = reader.get_bootstrap_n()
        form.edge_weight.form.bootstrap_first_random_seed_n = reader.get_bootstrap_first_random_seed_n()
        form.edge_weight.form.bootstrap_random_seed_update_parameter_n = reader.get_bootstrap_random_seed_update_parameter_n()
        form.edge_weight.form.process()
    if 'PC_ALGORITHM' in reader.keys():
        form.pc_algorithm.form.causal_discovery_observation_n = reader.get_causal_discovery_observation_n()
        form.pc_algorithm.form.indeptest = reader.get_indeptest()
        form.pc_algorithm.form.alpha = reader.get_alpha()
        form.pc_algorithm.form.numcores = reader.get_numcores()
        form.pc_algorithm.form.verbose = reader.get_verbose()
        form.pc_algorithm.form.nadelete = reader.get_nadelete()
        form.pc_algorithm.form.m_max = reader.get_m_max()
        form.pc_algorithm.form.u2pd = reader.get_u2pd()
        form.pc_algorithm.form.skel_method = reader.get_skel_method()
        form.pc_algorithm.form.conservative = reader.get_conservative()
        form.pc_algorithm.form.maj_rule = reader.get_maj_rule()
        form.pc_algorithm.form.solve_confl = reader.get_solve_confl()
        form.pc_algorithm.form.fixedgaps = reader.get_fixedgaps()
        form.pc_algorithm.form.fixededges = reader.get_fixededges()
        form.pc_algorithm.form.process()
    if 'PLOT_AND_DISPLAY' in reader.keys():
        form.plot_and_display.form.core_plot_title_str = reader.get_core_plot_title_str()
        form.plot_and_display.form.process()
