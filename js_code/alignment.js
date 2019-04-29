function alignmentObject () {
    //Variables
    this.linebreak = "\n";
    this.type = "alignment";
    this.OTUs={};
    this.partition_order=[];
    this.partitions={};

    //Functions
    this.pars_fasta = function (text,file_name,prefix="",sufix="") {
	var sequences = text.split('>');
        var max_length = 0;
        for (var i=0; i < sequences.length; ++i) {
            if (sequences[i]) {
                var seq = sequences[i].split(/\r\n|\n|\r/,);
                var name = seq.shift();
		if (prefix) {
		    var q=0;
		    for (var p = 0; p < name.length; ++p) {
			if (name[p] === prefix[0]) {
			    while (q < prefix.length && name[p+q] === prefix[q]) { ++q; }
			    if (q === prefix.length && p+q < name.length) {
				name = name.substring(p+q);
				break;
			    }
			}
		    }
		}
		if (sufix) {
		    //alert (sufix);
		    var last_match=0;
		    for (var p = 0; p < name.length; ++p) {
			if (name[p] === sufix[0]) {
			    var q=0;
			    while (q < sufix.length && name[p+q] === sufix[q]) { ++q; }
			    if (q === sufix.length && p+q < name.length) {
				last_match = p;
			    }
			    //alert ( last_match );
			}
		    }
		    if (last_match > 0) { name = name.substring(0,last_match); }
		}
                while (seq[0] === undefined || seq[0] === null || seq[0] === "") { seq.shift(); }
                seq = seq.join("");
                seq.replace(/\s/g,"");
                if (seq.length > max_length) { max_length = seq.length; }
                this.add_seq(name,file_name,seq);
            }
        }
	this.add_partition(file_name, max_length);
	return this.nOTUs();
    }
    this.alignment_length = function () {
	length = 0;
	for (part in this.partitions) {
	    if (this.partitions.hasOwnProperty(part)) { length += this.partitions[part]; }
	}
	return length;
    }
    this.nOTUs = function () {
	var n=0;
	for (OTU in this.OTUs) {
	    if (this.OTUs.hasOwnProperty(OTU)) { ++n; }
	}
	return n;
    }
    this.add_partition = function (partition, length=0) {
	this.partitions[partition] = length;
    }
    this.print_all_partitions = function () {
	this.partition_order=[];
	for (part in this.partitions) {
	    if (this.partitions.hasOwnProperty(part)) {
		this.partition_order.push(part);
	    }
	}
    }
    this.get_all_partitions = function () {
	var output =[];
	for (part in this.partitions) {
	    if (this.partitions.hasOwnProperty(part)) {
		output.push(part);
	    }
	}
	return output;
    }
    this.get_partition_in_order_as_rows = function () {
	var output = "";
	var start = 1;
	for (var i=0; i < this.partition_order.length; ++i) {
	    output += "<tr><td>" + this.partition_order[i] + "</td><td>" + start + "-" + (start-1+this.partitions[this.partition_order[i]]) + "</td></tr>\n";
	    start += this.partitions[this.partition_order[i]];
	}
	return output;
    }
    this.get_all_partitions_as_table_rows = function () {
	var output = "";
	for (part in this.partitions) {
            if (this.partitions.hasOwnProperty(part)) {
                output += "<tr><td>" + part + "</td><td>" + this.partitions[part] + "</td></tr>\n";
            }
        }
        return output;
    }
    this.add_partition_to_order= function (partition) {
	if (this.partition_order.includes(partition)) { alert(partition + " is already included among the partitions."); }
	this.partition_order.push(partition);
    }
    this.add_seq = function (OTU, partition, sequence) {
	if (OTU !== undefined && partition !== undefined && sequence !== undefined) {
	    if (this.OTUs[OTU] === undefined) { this.OTUs[OTU] = {}; }
	    this.OTUs[OTU][partition] = sequence;
	}
    }
    this.split_partition = function (to_split, new_partitions) {
	for (OTU in this.OTUs) {
	    if (this.OTUs.hasOwnProperty(OTU)) {
		for (part in new_partitions) {
		    if (new_partitions.hasOwnProperty(part)) {
			if (this.OTUs[OTU][to_split]) {
			    for (i=0; i<new_partitions[part].length; i+=2) {
				var slash = new RegExp("/");
				if (slash.test(new_partitions[part][i+1])) {
				    interval = new_partitions[part][i+1].split(/\//);
				    this.OTUs[OTU][part] = "";
				    for (j=new_partitions[part][i]; j<interval[0] && j<this.OTUs[OTU][to_split].length;j+=interval[1]) {
					this.OTUs[OTU][part] += this.OTUs[OTU][to_split][j-1];
				    }
				}
				else {
				    this.OTUs[OTU][part] = this.OTUs[OTU][to_split].substring(new_partitions[part][0]-1,new_partitions[part][1]);
				}
			    }
			    if (this.partitions[part] === undefined) { this.partitions[part] = this.OTUs[OTU][part].length; }
			    else if (this.partitions[part].length < this.OTUs[OTU][part].length) { this.partitions[part] = this.OTUs[OTU][part].length; }
			}
		    }
		}
		delete this.OTUs[OTU][to_split];
	    }
	}
	delete this.partitions[to_split];
    }
    this.getGC = function () {
	var nGC=0;
	var n=0;
	for (OTU in this.OTUs) {
	    if (this.OTUs.hasOwnProperty(OTU)) {
		for (var i = 0; i < this.partition_order.length; ++i) {
	    	    if (this.OTUs[OTU][this.partition_order[i]] !== undefined) {
	    		for (var j=0; j < this.OTUs[OTU][this.partition_order[i]].length; ++j) {
			    ++n;
			    if (this.OTUs[OTU][this.partition_order[i]][j] === "G") { ++nGC; }
			    else if (this.OTUs[OTU][this.partition_order[i]][j] === "C") { ++nGC; }
			    else if (this.OTUs[OTU][this.partition_order[i]][j] === "g") { ++nGC; }
			    else if (this.OTUs[OTU][this.partition_order[i]][j] === "c") { ++nGC; }
			}
		    }
		}
	    }
	}
	return nGC/n;
    }
    this.get_sequences_as_fasta = function () {
	var fasta = "";
	var lengths = [];
	for (OTU in this.OTUs) {
	    for (var i =0; i < this.partition_order.length; ++i) {
		if (OTU[this.partition_order[i]] !== undefined && length[i] < OTU[this.partition_order[i]].length) {
		    length[i] = OTU[this.partition_order[i]].length;
		}
	    }
	}
	for (OTU in this.OTUs) {
	    if (this.OTUs.hasOwnProperty(OTU)) {
		fasta += ">" + OTU + this.linebreak;
		for (var i = 0; i < this.partition_order.length; ++i) {
		    if (this.OTUs[OTU][this.partition_order[i]] !== undefined) {
			fasta += this.OTUs[OTU][this.partition_order[i]];
			if (this.OTUs[OTU][this.partition_order[i]].length < this.partitions[this.partition_order[i]]) {
			    fasta += "-".repeat(this.partitions[this.partition_order[i]]-this.OTUs[OTU][this.partition_order[i]].length);
		       	}
		    }
		    else {
			fasta += "-".repeat(this.partitions[this.partition_order[i]]);
		    }
		}
		fasta += this.linebreak;
	    }
	}
	return(fasta);
    }
    this.bootRep = function () {
	var rep = new alignmentObject();
	rep.linebreak = this.linebreak;
	var chars = {};
	for (var OTU in this.OTUs) {
	    if (this.OTUs.hasOwnProperty(OTU)) {
		rep.OTUs[OTU] = {};
		for (var pars in this.OTUs[OTU]) {
		    if (this.OTUs[OTU].hasOwnProperty(pars)) {
			if (!chars.hasOwnProperty(pars)) {
			    chars[pars] = []
			    for (var col=0; col < this.OTUs[OTU][pars].length; ++col) {
				var selected = Math.floor(Math.random()*this.OTUs[OTU][pars].length);
				chars[pars].push(selected);
				//console.log(selected);
			    }
			}
			rep.OTUs[OTU][pars] = "";
			for (var i = 0; i < chars[pars].length; ++i) {
			    //console.log();
			    rep.OTUs[OTU][pars] += this.OTUs[OTU][pars].charAt(chars[pars][i]);
			}
		    }
		}
	    }
	}
	rep.partition_order=this.partition_order;
	rep.partitions=this.partitions;
	return(rep);
    }
}

function alphabet () {
    this.nonASuncertain = true;
    this.characters = {};
    this.uncertain;
    this.set_uncertainty = function (uncertainty) {
	this.uncertain = uncertainty;
	for (var chars in this.characters) {
	    if (this.characters.hasOwnProperty(chars)) {
		addBit(this.characters[chars],this.characters[uncertainty]);
	    }
	}
    }
    function setBit (bit) {
	var value=0;
	if (bit > 0) {
	    value = 1;
	    value = value << (bit-1);
	}
	return value;
    }
    function addBit (bit, bit_store) {
	var value = setBit (bit);
	bit_store = value | bit_store;
    }
    this.DNA = function (gapsASfifth = false) {
	this.uncertain = 'N';
	this.characters['A'] = this.characters['a'] = setBit(1);
	this.characters['C'] = this.characters['c'] = setBit(2);
	this.characters['G'] = this.characters['g'] = setBit(3);
	this.characters['T'] = this.characters['t'] = this.characters['U'] = this.characters['u']= setBit(4);
	this.characters['R'] = this.characters['r'] = this.characters['A'] | this.characters['G'];
	this.characters['Y'] = this.characters['y'] = this.characters['C'] | this.characters['T'];
	this.characters['S'] = this.characters['s'] = this.characters['G'] | this.characters['C']; 
	this.characters['W'] = this.characters['w'] = this.characters['A'] | this.characters['T'];
       	this.characters['K'] = this.characters['k'] = this.characters['G'] | this.characters['T'];
       	this.characters['M'] = this.characters['m'] = this.characters['A'] | this.characters['C'];
       	this.characters['B'] = this.characters['b'] = this.characters['C'] | this.characters['G'] | this.characters['T'];
	this.characters['D'] = this.characters['d'] = this.characters['A'] | this.characters['G'] | this.characters['T'];
	this.characters['H'] = this.characters['h'] = this.characters['A'] | this.characters['C'] | this.characters['T'];
       	this.characters['V'] = this.characters['v'] = this.characters['A'] | this.characters['C'] | this.characters['G'];
	this.characters['N'] = this.characters['n'] = this.characters['A'] | this.characters['C'] | this.characters['G'] | this.characters['T'];
	if (gapsASfifth) {
	    this.characters['.'] = this.characters['-'] = setBit(5);
	    this.characters['?'] = this.characters['N'] | this.characters['-'];
	}
	else {
	    this.characters['.'] = this.characters['-'] = 0;
	    this.characters['?'] = this.characters['N'];
	}
    }
    this.get_char = function (bit) {
	for (chars in this.characters) {
	    if (this.characters.hasOwnProperty(chars) && chars !== 'uncertain' && this.characters[chars] === bit) { return chars;}
	}
	return this.uncertain;
    }
    this.pars_compare = function (A, B) {
	if (this.nonASuncertain && this.uncertain !== undefined && this.uncertain !== null) {
	    if (this.characters[A] === 0) A = this.uncertain;
	    if (this.characters[B] === 0) B = this.uncertain;
	}
	if (this.characters[B] === undefined || this.characters[B] === null) { alert("No char " + B); return [A,0]; }
	else if (this.characters[A] === undefined || this.characters[A] === null) { alert("No char " + A); return [B,0]; }
	else if (this.characters[B] & this.characters[A]) { return [this.get_char(this.characters[B] & this.characters[A]),0]; }
	else { return [this.get_char(this.characters[B] | this.characters[A]),1]; }
    }
}
