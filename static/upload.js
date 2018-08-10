$(document).ready(function () {
    var table = $('#upload').DataTable({
        'select': 'multiple',
        'searching': false,
        fixedHeader: true

    });

    var dataset_table;

    var counter = 1;

    $('#addRow').on('click', function () {
        table.row.add([
            `<input type=\"file\" name=input-${counter}/>`,
            `<input type=\"file\"  name=factor-${counter}/>`
        ]).draw(true);

        counter++;
    });

    var visible = true;
    var tableContainer = $(table.table().container());
    var dataset_tableContainer;


    $('#removeRow').on('click', function () {
        table.rows({selected: true}).every(function (rowIdx, tableLoop, rowLoop) {
            table.row('.selected').remove().draw(true);
        });
    });


    $('form').submit(function () {
        let selected_rows = [];
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


    handle_key = appConfig.handle_key;
    var first_dataset_configs = handle_key.user_configs[handle_key.user_datasets[0]];
    dataset_table = $('#dataset-table').DataTable({
        data: [first_dataset_configs],
        'select': 'multiple',
        columns: [{title: 'config'}],
        fixedHeader: true
    }).draw(false);

    dataset_tableContainer = $(dataset_table.table().container());
    dataset_tableContainer.css('display', 'none');
    dataset_table.fixedHeader.adjust();

    $('#existingdataset').on('click', function () {

        $('#existing-select').toggle();

        tableContainer.css('display', visible ? 'none' : 'block');
        table.fixedHeader.adjust();

        dataset_tableContainer.css('display', visible ? 'block' : 'none');
        dataset_table.fixedHeader.adjust();

        visible = !visible;
    });


});

