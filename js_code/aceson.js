class aceson {
    constructor ( ) {
	this.contig = []; //JSON.parse( JSONstring );
	this.colors = { "A": "green",
			"C": "blue",
			"G": "black",
			"T": "red" }
	console.log(this.contig.length);
    }
    addContigs ( JSONstring ) {
	//console.log(JSON.parse( JSONstring ));
	let to_add = JSON.parse( JSONstring );
	if (to_add["contigs"])
	    this.contig = this.contig.concat(JSON.parse( JSONstring )["contigs"]);
	else if (to_add["reads"]) {
	    this.contig = this.contig.concat(JSON.parse( JSONstring ));
	}
	else
	    console.log("No recognized format");
	console.log(this.contig);
    }
    haveTraces( contig_no ) {
	if (this.contig[contig_no]["reads"])
	    for (let i=0; i<this.contig[contig_no]["reads"].length; ++i)
		if (this.haveTrace(contig_no, i))
		    return true;
	return false;
    }
    haveTrace ( contig_no, read_no) {
	if (this.contig[contig_no]["reads"][read_no]["chrom"] && this.contig[contig_no]["reads"][read_no]["chrom"]["traces"])
	    return true;
	return false;
    }
    haveAlignment( contig_no ) {
	if (this.contig[contig_no]["alignment"])
	    return true;
	return false;
    }
    getNcontigs () { return this.contig.length; }
    getNreads ( contig_no ) {
	if (contig_no < this.contig.length && this.contig[contig_no]["reads"])
	    return this.contig[contig_no]["reads"].length;
	else
	    return 0;
    }
    getTraceMaxval( contig_no, read_no ) {
	//console.log(this.contig[contig_no]);
	if (contig_no < this.contig.length && read_no < this.contig[contig_no]["reads"].length)
	    return (this.contig[contig_no]["reads"][read_no]["chrom"]["traces"]["max_val"]);
	else
	    return null;
    }
    getContigPos ( contig_no, pos ) {
	return pos - this.contig[contig_no]["alignment"]["contig"]["start"];
    }
    getContigSeqPos ( contig_no, pos ) {
	let contig_pos = this.getContigPos( contig_no, pos );
	let n_prev_gap = 0;
	for (let i=0; i<this.contig[contig_no]["alignment"]["contig"]["gaps"].length; ++i) {
	    if (contig_pos >= this.contig[contig_no]["alignment"]["contig"]["gaps"][i]["pos"] && contig_pos < this.contig[contig_no]["alignment"]["contig"]["gaps"][i]["pos"]+this.contig[contig_no]["alignment"]["contig"]["gaps"][i]["length"])
		return -1;
	    else if (contig_pos > this.contig[contig_no]["alignment"]["contig"]["gaps"][i]["pos"])
		n_prev_gap += this.contig[contig_no]["alignment"]["contig"]["gaps"][i]["length"];
	    else
    		break;
	}
	return pos - this.contig[contig_no]["alignment"]["contig"]["start"] - n_prev_gap;

    }
    getReadSeqPos ( contig_no, read_no, pos ) {
        let read_pos = pos - this.contig[contig_no]["alignment"]["reads"][read_no]["start"];
        let n_prev_gap = 0;
        for (let i=0; i<this.contig[contig_no]["alignment"]["reads"][read_no]["gaps"].length; ++i) {
            if (read_pos >= this.contig[contig_no]["alignment"]["reads"][read_no]["gaps"][i]["pos"] && read_pos < this.contig[contig_no]["alignment"]["reads"][read_no]["gaps"][i]["pos"]+this.contig[contig_no]["alignment"]["reads"][read_no]["gaps"][i]["length"])
                return -1;
            else if (read_pos > this.contig[contig_no]["alignment"]["reads"][read_no]["gaps"][i]["pos"])
		n_prev_gap += this.contig[contig_no]["alignment"]["reads"][read_no]["gaps"][i]["length"];
            else
                break;
        }
	if (read_pos-n_prev_gap >= 0 && read_pos-n_prev_gap < this.contig[contig_no]["reads"][read_no]["chrom"]["n_bases"]) {
	    if (this.contig[contig_no]["reads"][read_no]["reverse"])
		return this.contig[contig_no]["reads"][read_no]["chrom"]["n_bases"]-(read_pos-n_prev_gap)-1;
	    else
		return read_pos-n_prev_gap;
	}
	else
	    return -2;
    }
    getReadBase( contig_no, read_no, pos ) {
	let read_pos = this.getReadSeqPos( contig_no, read_no, pos );
	if (read_pos === -1) return {"base":'*', "qual":0};
	else if (read_pos < 0) return {"base":'-', "qual":0};
	else if (this.contig[contig_no]["reads"][read_no]["reverse"]) {
	    //read_pos = this.contig[contig_no]["reads"][read_no]["chrom"]["n_bases"]-read_pos-1;
	    return {"base":this.complementTo(this.contig[contig_no]["reads"][read_no]["chrom"]["calls"]["seq"][read_pos]), "qual":this.contig[contig_no]["reads"][read_no]["chrom"]["calls"]["base_prob"][read_pos][this.contig[contig_no]["reads"][read_no]["chrom"]["calls"]["seq"][read_pos]]};
	    
	}
	else
	    return {"base":this.contig[contig_no]["reads"][read_no]["chrom"]["calls"]["seq"][read_pos], "qual":this.contig[contig_no]["reads"][read_no]["chrom"]["calls"]["base_prob"][read_pos][this.contig[contig_no]["reads"][read_no]["chrom"]["calls"]["seq"][read_pos]]};
    }
    getContigBase( contig_no, pos ) {
	let contig_pos = this.getContigPos( contig_no, pos );
	let n_prev_gap = 0;
	if (contig_pos >= 0) {
	    for (let i=0; i<this.contig[contig_no]["alignment"]["contig"]["gaps"].length; ++i) {
		console.log("" + contig_pos + " " + this.contig[contig_no]["alignment"]["contig"]["gaps"][i]["pos"] + " " + this.contig[contig_no]["alignment"]["contig"]["gaps"][i]["pos"] + this.contig[contig_no]["alignment"]["contig"]["gaps"][i]["length"]);
		if (contig_pos >= this.contig[contig_no]["alignment"]["contig"]["gaps"][i]["pos"] && contig_pos < this.contig[contig_no]["alignment"]["contig"]["gaps"][i]["pos"] + this.contig[contig_no]["alignment"]["contig"]["gaps"][i]["length"])
		    return {"base":'*', "qual":this.contig[contig_no]["alignment"]["contig"]["gaps"][i]["qual"][contig_pos-this.contig[contig_no]["alignment"]["contig"]["gaps"][i]["pos"]]};
		else if (contig_pos > this.contig[contig_no]["alignment"]["contig"]["gaps"][i]["pos"])
		    n_prev_gap += this.contig[contig_no]["alignment"]["contig"]["gaps"][i]["length"];
		else
		    break;
	    }
	    if (contig_pos-n_prev_gap < this.contig[contig_no]["seq"].length)
		return {"base":this.contig[contig_no]["seq"][contig_pos-n_prev_gap],"qual":this.contig[contig_no]["qual"][contig_pos-n_prev_gap]};
	    else
		return {"base":'-',"qual":0};
	}
	else
	    return {"base":'-',"qual":0};
    }
    setContigBase( contig_no, pos, base ) {
	let contig_pos = this.getContigSeqPos( contig_no, pos );
	if (contig_pos >= 0 && contig_pos < this.contig[contig_no]["seq"].length) {
	    this.contig[contig_no]["seq"][contig_pos] = base;
	    this.contig[contig_no]["qual"][contig_pos] = 100;
	}
    }
    tracesSVG ( contig_no, read_no, start_x, start_y, scale_x, scale_y ) {
	//console.log("Giving svg for read " + read_no + " (" + this.contig[contig_no]["reads"].length + ") of contig " + contig_no + " (" + this.contig.length + ")");
	if (contig_no < this.contig.length && this.contig[contig_no]["reads"] && read_no < this.contig[contig_no]["reads"].length) {
	    let svg_string = "<svg id=\"chrom_graph_" + read_no + "\" width=\"" + (this.contig[contig_no]["reads"][read_no]["chrom"]["traces"]["n_points"]*scale_x) + "\" height=\"" + (this.contig[contig_no]["reads"][read_no]["chrom"]["traces"]["max_val"]*scale_y) + "\">\n";
	    svg_string += "<line id=\"line_" + read_no + "\" x1=\"0\" y1=\"0\" x2=\"0\" y2=\"0\"style=\"stroke:rgb(0,0,0);stroke-width:1\" />\n";
	    for (let base of ["A","C","G","T"]) {
		svg_string += "<path d=\" M" + start_x + " " + start_y;
		//console.log(this.contig[contig_no]["reads"][read_no]["chrom"]["traces"][base].length);
		for (let i=0; i < this.contig[contig_no]["reads"][read_no]["chrom"]["traces"][base].length; ++i) {
		    let y_coord;
		    if (this.contig[contig_no]["reads"][read_no]["reverse"])
			y_coord = this.contig[contig_no]["reads"][read_no]["chrom"]["traces"][base][this.contig[contig_no]["reads"][read_no]["chrom"]["traces"][base].length-i-1];
		    else
			y_coord = this.contig[contig_no]["reads"][read_no]["chrom"]["traces"][base][i];
		    svg_string += " L" + ((i*scale_x)+start_x) + " " + ((start_y- y_coord* scale_y));
		}
		if (this.contig[contig_no]["reads"][read_no]["reverse"])
		    svg_string += "\" stroke=\"" + this.colors[this.complementTo(base)] + "\" fill=\"none\" />\n"
		else
		    svg_string += "\" stroke=\"" + this.colors[base] + "\" fill=\"none\" />\n"
	    }
	    svg_string += "</svg>\n"
	    return svg_string;
	}
	else
	    return "<svg></svg>";

    }
    contigAlignStart ( contig_no ) {
	let start_pos = 0;
	for (let i = 0; i < this.contig[contig_no]["n_reads"]; ++i) {
            if (this.contig[contig_no]["reads"][i]["start"] < start_pos)
                start_pos = this.contig[contig_no]["reads"][i]["start"];
        }
	return start_pos*-1;
    }
    complementTo ( base ) {
	if (base === 'A') return 'T';
	else if (base === 'C') return 'G';
	else if (base === 'G') return 'C';
	else if (base === 'T') return 'A';
	else return '*';
    }
    getQualCoord(contig_no, read, pos) {
	let read_pos = this.getReadSeqPos( contig_no, read, pos );
	while (read_pos === -1)
	    read_pos = this.getReadSeqPos( contig_no, read, ++pos );
	if (this.contig[contig_no]["reads"][read]["reverse"] )
	    return this.contig[contig_no]["reads"][read]["chrom"]["traces"]["n_points"]-this.contig[contig_no]["reads"][read]["chrom"]["calls"]["trace_pos"][read_pos];
	else
	    return this.contig[contig_no]["reads"][read]["chrom"]["calls"]["trace_pos"][read_pos];
    }
    alignPosSVG(contig_no, pos, fontsize) {
	let svg_string = "";
	let x_coord = pos*fontsize;
	for (let i=0; i < this.contig[contig_no]["n_reads"]; ++i) {
	    let read_base = this.getReadBase(contig_no, i, pos);
	    if (read_base["base"] !== '-') {
		if (read_base["base"] === '*')
		    console.log(read_base);
		let read_pos = this.getReadSeqPos(contig_no, i, pos);
		read_base["qual"] *= 4;
		if (read_base["qual"] > 255)
		    read_base["qual"] = 255;
		svg_string += "<rect x=\"" + (x_coord-fontsize/2) + "\" y=\"" + ((fontsize+1)*(i)) + "\" height=\"" + (fontsize+1) + "\" width=\"" + fontsize + "\" style=\"fill:rgb(" + read_base["qual"] + "," + read_base["qual"] + "," + read_base["qual"] + ")\"/>"
		svg_string += "<text x=\"" + x_coord + "\" y=\"" + (fontsize+1)*(i+1) + "\" text-anchor=\"middle\" fill=\"";
		let white_cut = 120;
		if (read_base["qual"] > white_cut && (read_pos < this.contig[contig_no]["reads"][i]["qual_begin"]-1 || read_pos > this.contig[contig_no]["reads"][i]["qual_end"])) svg_string += "gray";
	       	else if (read_base["base"] !== '*' && read_base["qual"] > white_cut) svg_string += this.colors[read_base["base"]];
		else if (read_base["qual"] <= white_cut) svg_string += "white";
		else svg_string += "black";
	       	svg_string += "\" font-size=\"" + fontsize + "\">" + read_base["base"] + "</text>";
	    }
	}
	//let contig_pos = pos - contig_start -2;
	let contig_base = this.getContigBase(contig_no, pos);
	if (contig_base["base"] !== '-') {
	    //let base = this.contig[contig_no]["seq"][contig_pos];
	    svg_string += "<text x=\"" + x_coord + "\" y=\"" + ((fontsize+1)*(this.contig[contig_no]["n_reads"]+1)) + "\" text-anchor=\"middle\" fill=\"" + this.colors[contig_base["base"]] + "\" font-size=\"" + fontsize + "\" onclick=\"baseClick(" + contig_no + ", " + pos + ")\">" + contig_base["base"] + "</text>";
	}
	return svg_string;
    }
    alignmentSVG ( contig_no, fontsize ) {
	let SVG = [("<svg width=\"" + (this.contig[contig_no]["alignment"]["length"]*fontsize) + "\" height=\"" + ((fontsize+2)*(this.contig[contig_no]["n_reads"]+1) + "\">\n"))];
	for (let i = 0; i <= this.contig[contig_no]["alignment"]["length"]; ++i)
	    SVG.push(this.alignPosSVG(contig_no, i, fontsize));
	SVG.push("</svg>\n");
	return SVG;
    }
}
