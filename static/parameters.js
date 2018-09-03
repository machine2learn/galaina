$(document).ready(function () {
    document.getElementById('edge_weight_algorithm-bootstrap_first_random_seed_n').classList.add('hidden');
    document.getElementById('edge_weight_algorithm-bootstrap_random_seed_update_parameter_n').classList.add('hidden');
    document.getElementById("copula_factor_algorithm-gibbs_burn_in_n").classList.add('hidden');
    document.getElementById("copula_factor_algorithm-gibbs_first_random_seed_n").classList.add('hidden');
    document.getElementById("copula_factor_algorithm-gibbs_random_seed_update_parameter_n").classList.add('hidden');
    document.getElementById("pc_algo").classList.add('hidden');

    $('#edge_weight_check').click(function () {
        if (document.getElementById("edge_weight_check").checked === false) {
            document.getElementById('edge_weight_algorithm-bootstrap_first_random_seed_n').classList.add('hidden');
            document.getElementById('edge_weight_algorithm-bootstrap_random_seed_update_parameter_n').classList.add('hidden');
        } else {
            document.getElementById('edge_weight_algorithm-bootstrap_first_random_seed_n').classList.remove('hidden');
            document.getElementById('edge_weight_algorithm-bootstrap_random_seed_update_parameter_n').classList.remove('hidden');
        }
    });

    $('#copula_factor_check').click(function () {
        if (document.getElementById("copula_factor_check").checked === false) {
            document.getElementById("copula_factor_algorithm-gibbs_burn_in_n").classList.add('hidden');
            document.getElementById("copula_factor_algorithm-gibbs_first_random_seed_n").classList.add('hidden');
            document.getElementById("copula_factor_algorithm-gibbs_random_seed_update_parameter_n").classList.add('hidden');
        } else {
            document.getElementById("copula_factor_algorithm-gibbs_burn_in_n").classList.remove('hidden');
            document.getElementById("copula_factor_algorithm-gibbs_first_random_seed_n").classList.remove('hidden');
            document.getElementById("copula_factor_algorithm-gibbs_random_seed_update_parameter_n").classList.remove('hidden');

        }
    });

    $('#pc_algo_check').click(function () {
        if (document.getElementById("pc_algo_check").checked === false) {
            document.getElementById("pc_algo").classList.add('hidden');
        } else {
            document.getElementById("pc_algo").classList.remove('hidden');
        }
    });

    // $("form").submit(function () {
    //     if (document.getElementsByName('copula_factor_algorithm-gibbs_burn_in_n')[0].value < document.getElementsByName('copula_factor_algorithm-gibbs_sampling_n')[0].value) {
    //         alert('ERROR :  Gibbs burn in < Gibbs sampling is not TRUE')  ;
    //         return false;
    //     }
    // });

});