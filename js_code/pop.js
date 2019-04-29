function individual (genome = [], father=null, mother=null, children=[]) {
    this.genome = genome;
    this.father = father;
    this.mother = mother;
    this.children = children;
    this.get_alleles(allels,chrom=0,l=0) {
	for (var i=0; i < this.genome.length; ++i) {
	    this.genome[i].get_alleles(alleles,chrom,l);
	}
    }
    this.generate_gamet() {
	var gamet_genome = new genome();
	var n_chrom = Number.MAX_VALUE;
	for (var i = 0; i < this.genome.length; ++i) {
	    if (this.genome[i].chromosomes.length < n_chrom) { n_chrom = this.genome[i].chromosomes.length; }
	}
	for (i=0; i < n_chrom; ++i) {
	    gamet_genome.push
	}
    }
}

function genome ( chromosomes = [] ) {
    this.chromosomes = chromosomes;
    this.get_alleles(allels, chrom=0, l=0) {
    	this.chromosomes[chrom].get_alleles(alleles,l);
    }
    this.add_chromosome(chrom) {
	this.chromosomes.push(chrom);
    }
}

function chromosome (loci=[]) {
    this.loci = loci;
    this.add_gene(allele) {
	this.loci.push(allele);
    }
    this.get_alleles(alleles, l) {
	alleles.add_allele(this.loci[l]);
    }
    this.clone() = 
}

function allele_counter () {
    this.alleles = {};
    this.add_allele(allele) {
	++this.alleles[allele];
    }
    this.get_alleles= function () {
	return_array = [];
	for (var k in this.alleles) {
	    if (this.alleles.hasOwnProperty(k)) {
		return_array.push(k);
	    }
	}
	return return_array;
    }
}

function generation (individuals = []); //, ploidy = 2, alleles = [1,2]) {
    this.individuals = individuals;
    this.start_generation = function (n=100, ploidy = 2, alleles = [1,2]) {
	this.individuals.length=0;
	for (var i = 0; i < n; ++i) {
	    gen = [];
	    for (var j=0; j < ploidy; ++j) {
		gen.push(new chromosome [alleles[Math.floor(Math.random()*alleles.length)]]);
	    }
	    this.individuals.push(new individual(gen));
	}
    }
    this.get_allele_freq = function () {
	var alleles = new allele_counter();
	for (var i=0; i < this.individuals.length; ++i) {
	    this.individual[i].add_alleles(alleles);
	}
	return alleles;
    }
    this.new_generation = function (parents, n=null, outcrosser = true, sex = null ) {
	if (n === null || n === undefined) {
	    n=parents.length;
	}
	this.individuals.length=0;
	var alleles = parents.get_allele_freq().get_alleles();
	for (var i = 0; i < n; ++i) {
	    var father;
	    var mother;
	    if (sex === null && outcrosser) {
		var one = Math.floor(Math.random()*parents.individuals.length);
		mother = parents.individuals[one];
		var two = Math.floor(Math.random()*parents.individuals.length-1);
		if (two >= one) { ++two; }
		father = parents.individuals[two];
	    }
	}

    }
}

function population (generations = []) {
    this.generations = generations;
    this.init_new = function (n=100, ploidy = 2, alleles = [1,2]) {
	this.generations.length = 0;
	var start = new generation ();
	start.start_generation(n, ploidy, alleles);
	this.generations[0] = start;
    }
    this.add_generation = function (parents, n=null) {
	

    }
}
