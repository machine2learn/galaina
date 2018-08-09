$(document).ready(function () {
    var table = $('#upload').DataTable({
        'select': 'multiple',
        'searching': false,
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


    $('#existing-select').hide();

    var handle_key = {};

    $.ajax({
        url: "/dataset_configs",
        type: 'GET',
        dataType: 'json',
        contentType: 'application/json;charset=UTF-8',
        accepts: {
            json: 'application/json',
        },
        success: function (data) {
            items = data.user_dataset[0];
            $.each(items, function (i, item) {
                $('#existing-select').append($('<option>', {
                    value: item,
                    text: item
                }));
            });
            handle_key.user_datasets = data.user_dataset[0];
            handle_key.user_configs = data.user_configs;

        }

    });

    $('#existingdataset').on('click', function () {

        $('#existing-select').toggle();
        // $('#dataset-table-div').toggle();

        if (visible) {
            table.destroy();
            $('#upload').empty();

            var first_dataset_configs = handle_key.user_configs[handle_key.user_datasets[0]];
            dataset_table = $('#dataset-table').DataTable({
                data: [first_dataset_configs],
                columns: [{title: 'config'}],
            });
        }
        else {
            table = $('#upload').DataTable({
                select: 'multiple',
                searching: false,
                columns: [{title: 'Input data'}, {title: 'Factor Model'}],
                data: [[`<input type=\"file\" name=input-${counter}/>`, `<input type=\"file\"  name=factor-${counter}/>`]]
            });
            dataset_table.destroy();
            $('#dataset-table').empty();

        }
        visible = !visible;

    });


});

