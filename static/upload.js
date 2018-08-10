$(document).ready(function () {
    var table = $('#upload').DataTable({
        'select': 'multiple',
        'searching': false,
        fixedHeader: true

    });

    var tableContainer = $(table.table().container());
    var counter = 1;
    var visible = true;

    var handle_key = appConfig.handle_key;
    var first_dataset_configs = handle_key.user_configs[handle_key.user_datasets[0]];
    var dataset_table = $('#dataset-table').DataTable({
        data: [first_dataset_configs],
        'select': 'single',
        columns: [{title: 'config'}],
        fixedHeader: true
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

    $('#existing-select').hide();

    $('#existingdataset').on('click', function () {

        $('#existing-select').toggle();

        $('#datasetname-row').toggle();

        tableContainer.css('display', visible ? 'none' : 'block');
        table.fixedHeader.adjust();

        dataset_tableContainer.css('display', visible ? 'block' : 'none');
        dataset_table.fixedHeader.adjust();

        visible = !visible;
    });

    $('#existing-select').on('change', function () {
        var dataset_selected = "";
        $("#existing-select option:selected").each(function () {
            dataset_selected = $(this).text();
        });
        data = handle_key.user_configs[dataset_selected];
        data = data.map(x => [x]);
        dataset_table.clear().rows.add(data).draw();

    });

    $('form').submit(function () {
        let selected_config = [];
        dataset_table.rows({selected: true}).every(function (rowIdx, tableLoop, rowLoop) {
            selected_config.push(this.data()[0]);
        });
        let input = $("<input>")
            .attr("type", "hidden")
            .attr("name", "selected_config").val(JSON.stringify(selected_config));

        $('form').append($(input));
    });


});

