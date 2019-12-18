//console.log("Read table parser");

function tableParser ( input, sep = "\t" ) {
    const text=input;
    this.p = 0;
    this.sep = sep;
    function isNewline ( character ) {
	if (character === '\n' || character === '\r') return true;
	else return false;
    }
    this.atEnd = function () {
	if (this.p < text.length) return false;
	else return true;
    }
    this.fromStart = function () { this.p = 0; }
    this.getNextRow = function () {
	var columns = [];
	var word = "";
	//console.log("Get row");
	for (;this.p<text.length;++this.p) {
	    if (isNewline(text[this.p])) {
		while(this.p<text.length && isNewline(text[this.p])) ++this.p;
		if (word) { columns[columns.length] = word; word = ""; }
		if (columns.length > 0) return columns;
		if (skip) { skip = false; }
	    }
	    else if (this.sep.search(text[this.p]) >= 0) {
    		columns[columns.length] = word;
		word = "";
	    }
	    else if (text[this.p] != ' ') { word += text[this.p]; }
	    else if (word) { word += text[this.p]; }
	}
	return columns;
    }
}
