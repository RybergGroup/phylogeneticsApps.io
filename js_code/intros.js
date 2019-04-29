function bio_apps_intro(field) {
    var text = "<a href=\"../index.html\">HOME</a>\n";
    text += "<p>This is a javascript application and everything is done on your computer by your web browser. It is free to use but no waranty, no guarantees, and no liability of the provider. The maximum size of blobs for your browser will set a limit to how big files can be processed. It has been tested in Firefox, but only to a limited extent in other browsers.</p>";
    field.innerHTML = text;
}

function set_style_bio (field) {
    var text = "table, th, td {\n";
    text += "border: 1px solid black;\n";
    text += "}\n";
    text += "td {\n";
    text += "text-align: center;\n";
    text += "}\n";
    field.innerHTML = text;
}
