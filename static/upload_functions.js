function isEmpty(obj) {
  return Object.keys(obj).length === 0;
}

function get_rows(dataset_name, user_configs, param_configs) {
    var my_rows = [];
    if (dataset_name in user_configs) {
        var len = user_configs[dataset_name].length;
        for (var i = 0; i < len; i++) {
            var config_name = user_configs[dataset_name][i];
            var bootstrap = param_configs[dataset_name + '_' + config_name]['bootstrap'];
            var gibbs_sampling_n = param_configs[dataset_name + '_' + config_name]['gibbs_sampling_n'];
            var gibbs_burn_in_n = param_configs[dataset_name + '_' + config_name]['gibbs_burn_in_n'];
            var delet_input = '<a data-id=' + config_name + ' onclick="ConfirmDelete(this, false)" ><span class="glyphicon glyphicon-remove"></span></a>';
            my_rows.push([config_name, bootstrap,gibbs_sampling_n,gibbs_burn_in_n, delet_input]);
        }
        my_rows.push(['new_config', '', '','','']);
    }
    return my_rows;
}


function validate_new_submit() {
    var datasetname = $('#datasetname').val();
    // Check dataset name ! contains '_'
    if (datasetname.indexOf('_') >= 0) {
        $('#error_text').text('Dataset name cannot contains "_".');
        return false;
    }
    // Check dataset name ! exists
    if (appConfig.handle_key.user_datasets.indexOf(datasetname) >= 0) {
        $('#error_text').text('Dataset name already exists.');
        return false;
    }
    // Check exists input data and factor model files
    var all_empty = true;
    for (var i = 0; i < document.getElementsByName('input-1/').length; i++) {
        if ((document.getElementsByName('input-1/')[i].files.length !== 0) || (document.getElementsByName('factor-1/')[i].files.length !== 0)) {
            all_empty = false;
        }
        if ((document.getElementsByName('input-1/')[i].files.length !== 0) && (document.getElementsByName('factor-1/')[i].files.length === 0)) {
            $('#error_text').text('Input data file without Factor model file');
            return false;
        }
        if ((document.getElementsByName('input-1/')[i].files.length === 0) && (document.getElementsByName('input-1/')[i].files.length !== 0)) {
            $('#error_text').text('Factor model file without Input data file');
            return false;
        }
    }
    // Check input data and factor model files not empty
    if (all_empty) {
        $('#error_text').text('Please choose files to upload');
        return false;
    }
    return true;
}


function validate_submit(submt) {
    if ($('#existingdataset').is(":checked") === false) {
        if (validate_new_submit() == false) {
            document.getElementsByName('datasetname')[0].remove();
            document.getElementsByName('selected_config')[0].remove();
            return false;
        }
    } else {
        // Check if there is a config selected
        if (submt == false) {
            $('#error_text').text('Please select a configuration.');
            document.getElementsByName('selected_config')[0].remove();
            datasetname
        }
        return submt;
    }
    return true;
}

function ConfirmDelete(elem, all) {
    var message = "Are you sure you want to delete the selected configuration?";
    if (all == true) {
        message = "Are you sure you want to delete all saved configurations?";
    }
    if (confirm(message)) {
        $.ajax({
            url: "/delete_config",
            type: 'POST',
            dataType: 'json',
            contentType: 'application/json;charset=UTF-8',
            accepts: {
                json: 'application/json',
            },
            data: JSON.stringify({
                'config': $(elem).attr('data-id'),
                'dataset': $('#existing-select').find("option:selected").text()
            }),
            success: function (data) {
                var dataset_selected = $('#existing-select').find("option:selected").text();
                handle_key.user_configs = data.user_configs;
                handle_key.param_configs = data.param_configs;
                data = get_rows(dataset_selected, handle_key.user_configs, handle_key.param_configs);
                $('#dataset-table').DataTable().clear().rows.add(data).draw();
            }
        })
    }
}

