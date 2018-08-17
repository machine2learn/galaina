$(document).ready(function () {
    $("#run_button").click(function () {
        console.log('you click in play');
        document.getElementById("span_run").classList.remove('glyphicon-play');
        document.getElementById("span_run").classList.add('glyphicon-cd');
        $('#run_button').prop('disabled', true);

        $.ajax({
            url: "/run",
            type: 'POST',
            data: $("#run_form").serialize(),
            success: function (data) {
                console.log('server ok ');
                document.getElementById("span_run").classList.remove('glyphicon-cd');
                document.getElementById("span_run").classList.add('glyphicon-ok');
                document.getElementById("run_button").classList.remove('btn-primary');
                document.getElementById("run_button").classList.add('btn-success');
            }
        })
    })

    var output = document.getElementById('log');
    setInterval(function () {
        $.ajax({
            url: "/stream",
            type: 'GET',
            dataType: 'json',
            contentType: 'application/json;charset=UTF-8',
            accepts: {
                json: 'application/json',
            },
            data: JSON.stringify(""),
            success: function (data) {
                if (data.data != '') {
                    output.append(data.data)
                    $('#log').scrollTop($('#log')[0].scrollHeight);
                }
            }
        })
    }, 1000);

});
