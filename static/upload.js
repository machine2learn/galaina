$(document).ready(function () {
    const table = $('#upload').DataTable({
        'select': 'multiple',
        'searching': false,

    });
    var counter = 1;

    $('#addRow').on('click', function () {
        table.row.add([
            `<input type=\"file\" name=input-${counter}/>`,
            `<input type=\"file\"  name=factor-${counter}/>`
        ]).draw(false);

        counter++;

    });

    $('#addRow').click();


    $('#removeRow').on('click', function () {
        table.rows({selected: true}).every(function (rowIdx, tableLoop, rowLoop) {
            table.row('.selected').remove().draw(false);
        });
    });


    $('form').submit(function () {
        let selected_rows = [];
        let input = $("<input>")
            .attr("type", "hidden")
            .attr("name", "datasetname").val($('#datasetname').val());
        $('form').append($(input));
    });

});

