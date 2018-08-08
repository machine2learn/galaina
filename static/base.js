function define_steps(page) {
    var element = document.getElementById(page);
    element.className = element.className.replace("btn btn-default", "btn btn-primary");
    var i;
    for (i = page + 1; i < 6; i++) {
        var nelement = document.getElementById(i);
        nelement.classList.add('disabled');
        nelement.disabled = true;
    }
};