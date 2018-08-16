$(document).ready(function () {
    if (isEmpty(appConfig.handle_key.user_configs)){
        $('#existingdataset-check').toggle();
    }

    var table = $('#upload').DataTable({
        'select': 'multiple',
        'searching': false,
        fixedHeader: false
    });

    var tableContainer = $(table.table().container());
    var counter = 1;
    var visible = true;

    var handle_key = appConfig.handle_key;
    var my_rows = get_rows(handle_key.user_datasets[0], handle_key.user_configs, handle_key.param_configs);
    var foot = 'Clear all configs <a data-id=all onclick="ConfirmDelete(this, true)"> <span class="glyphicon glyphicon-trash" style="color:#ff0000"></span></a>';


    var dataset_table = $('#dataset-table').DataTable({
        data: my_rows, 'select': 'single',
        columns: [{title: 'Config Name'}, {title: 'Bootstrap'}, {title: '', width: "5%"}],
        fixedHeader: true,
        footer: foot,
    }).draw(false);

    var dataset_tableContainer = $(dataset_table.table().container());
    dataset_tableContainer.css('display', 'none');
    dataset_table.fixedHeader.adjust();

    $('#addRow').on('click', function () {
        table.row.add([
            `<input type=\"file\" name=input-${counter}/>`,
            `<input type=\"file\"  name=factor-${counter}/>`
        ]).draw(true);

        counter++;
    });

    $('#removeRow').on('click', function () {
        table.rows({selected: true}).every(function (rowIdx, tableLoop, rowLoop) {
            table.row('.selected').remove().draw(true);
        });
    });


    $('form').submit(function () {
        let input = $("<input>")
            .attr("type", "hidden")
            .attr("name", "datasetname").val($('#datasetname').val());
        $('form').append($(input));
    });

    $.each(handle_key.user_datasets, function (i, item) {
        $('#existing-select').append($('<option>', {
            value: item,
            text: item
        }));
    });
    $('#existing-select-title').hide();
    $('#existing-select').hide();

    $('#existingdataset').on('change', function () {
        $('#existing-select').toggle();
        $('#existing-select-title').toggle();
        $('#datasetname-row').toggle();
        $('#buttons-row').toggle();

        tableContainer.css('display', visible ? 'none' : 'block');
        table.fixedHeader.adjust();

        dataset_tableContainer.css('display', visible ? 'block' : 'none');
        dataset_table.fixedHeader.adjust();

        visible = !visible;
        $('#error_text').empty();
    });

    $('#existing-select').on('change', function () {
        var dataset_selected = "";
        $("#existing-select option:selected").each(function () {
            dataset_selected = $(this).text();
        });
        // data = handle_key.user_configs[dataset_selected];
        // data = data.map(x => [x]);
        data = get_rows(dataset_selected, handle_key.user_configs, handle_key.param_configs);
        dataset_table.clear().rows.add(data).draw();
    });

    $('form').submit(function () {
        var submt = false;
        let selected_config = [];
        dataset_table.rows({selected: true}).every(function (rowIdx, tableLoop, rowLoop) {
            selected_config.push(this.data()[0]);
            submt = true;
        });
        let input = $("<input>")
            .attr("type", "hidden")
            .attr("name", "selected_config").val(JSON.stringify(selected_config));

        $('form').append($(input));
        return validate_submit(submt);
    });

});
