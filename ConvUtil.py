#!/usr/bin/python
#
# Amir Marcovitz
# 
# helper functions for mammalianCodingConvergence.py

import sys
import optparse
import gzip
import os

species_list = sorted(["hg38","ailMel1","bosTau8","calJac3","camFer1","canFam3","capHir1","cavPor3","cerSim1",
		"chiLan1","chlSab2","chrAsi1","conCri1","criGri1","dasNov3","echTel2","eleEdw1","eptFus1","equCab2",
		"eriEur2","felCat8","hetGla2","jacJac1","lepWed1","loxAfr3","macFas5","mesAur1","micOch1",
		"mm10","musFur1","myoDav1","myoLuc2","nomLeu3","ochPri3","octDeg1","odoRosDiv1","orcOrc1",
		"oryAfe1","oryCun2","otoGar3","oviAri3","panHod1","panTro5","papAnu2","ponAbe2","pteAle1","pteVam1","rheMac8",
		"rn6","saiBol1","sorAra2","speTri2","susScr3","triMan1","tupChi1","turTru2","vicPac2"])

legal_amino_acids = "ACDEFGHIKLMNPQRSTVWY"

def parse():
	'''
	  Read command-line arguments
	'''
        parser = optparse.OptionParser(description="Scan and detect convergent & divergent mutations in target species")
        parser.add_option('--tg', dest="target_groups", help="two comma-delimited species lists for targets", action="store")
        parser.add_option('--og', dest="outgroups", help="two comma-delimited species lists for outgroups", action="store")
        parser.add_option("--pc", dest="position_conservation", help="conservation minimal threshold (BBLS)", 
		type=float, action="store", default=0.9)
        args = parser.parse_args()
	if len(sys.argv) < 5:
		parser.print_help()
		sys.exit(1)
	species_assembly_map = read_species_assembly('data/species_list.txt.gz')
        target_groups = []
        data = args[0].target_groups.strip().split(":")
        for d in data:
		a = []
		for s in d.split(","):
			if s in species_assembly_map.values():
				a.append(s)
			elif s in species_assembly_map.keys():
				a.append(species_assembly_map[s])
			else:
				print('\n\tERROR: target species ' + s + ' not in screen\n\n')
				sys.exit(1)
		target_groups.append(a)
        outgroups = [] 
        data = args[0].outgroups.strip().split(":")
        for d in data:
		a = []
		for s in d.split(","):
			if s in species_assembly_map.values():
				a.append(s)
			elif s in species_assembly_map.keys():
				a.append(species_assembly_map[s])
			else:
				print('\n\tERROR: outgroup species ' + s + ' not in screen\n\n')
				sys.exit(1)
		outgroups.append(a)
        target_groups, outgroups = deep_sort(target_groups, outgroups)
        for tg in target_groups:
                for ts in tg:
                        assert ts in species_list
        for og in outgroups:
                for os in og:
                        assert os in species_list
	target_species = []
	for tg in target_groups:
		for ts in tg:
			target_species.append(ts)
			assert ts in species_list
	species_string = "_".join(target_species)
	return (target_groups,outgroups,args[0].position_conservation, target_species, species_string)


def read_species_assembly(filename):
	'''
	  Read the mapping between species names and genome assemblies
	'''
	d = {}
	with gzip.open(filename) as f:
		for line in f:
			(key, val) = line.split()
			d[key] = val
	return d


def deep_sort(targets, outgroups):
	'''
	  Output: Sorts the lower level of lists, and then sorts the
	         higher level. Useful for standardizing results.
	'''
	new_list = []
	for i in range(len(targets)):
		target_group = sorted(targets[i])
		outgroup = sorted(outgroups[i])
		new_list.append((target_group,outgroup))
	new_list = sorted(new_list)
	return [nl[0] for nl in new_list], [nl[1] for nl in new_list]


def load_RefGenes():
	'''
	  Load data for reference genes (human, GRCh38/hg38; Ensembl86)
	'''
	transcript_list = []
	transcript_to_gene = {}
	exon_locations = {}
	gene_symbols = {}
	# load list of transcripts to analyze
	with gzip.open('data/filteredTranscripts.hg38.txt.gz') as f:
		for line in f:
			if line[0]=="#":
				continue
			transcript_list.append(line.strip())	
	# load bed-file (genomic coordinates on GRCh38.hg38)
	genes_mapped = set([])
	with gzip.open('data/nonDuplicatedCanonical.hg38.bed.gz') as f:
		for line in f:
			data = line.strip().split()
			if data[3] in genes_mapped:
				if data[4] in transcript_to_gene:
					del genes_mapped[transcript_to_gene]			
			else:
				transcript_to_gene[data[4]] = data[3]
			genes_mapped.add(data[3])
	# load file with all exon locations
	with gzip.open('data/hg38.exonAlignments.bed.ref.gz') as f:
		for line in f:
			data = line.strip().split()
			chrom = data[0]
			start = int(data[1])
			end = int(data[2])
			transcript = data[3].split(".")[0]
			strand = data[5]
			if transcript not in exon_locations:
				exon_locations[transcript] = []
			# some of the genes in the input file
			# have the same exons listed multiple times.
			# Get rid of these duplicate exons.
			if (chrom, start, end, strand) not in exon_locations[transcript]:
				exon_locations[transcript].append((chrom, start, end, strand))
	# load mapping of gene-symbols to gene IDs
	with gzip.open('data/humanGeneSymbol.humanEnsembl.biomart86.NoSyn.map.gz') as f:
		for line in f:
			data = line.strip().split()
			gene_symbols[data[1]] = data[0]
	return transcript_list, transcript_to_gene, exon_locations, gene_symbols


def init_outputs(species_string, position_conservation):
	if not os.path.isdir(os.getcwd() + "/ConvDiv_sites"):
		os.mkdir(os.getcwd() + "/ConvDiv_sites")
	if not os.path.isdir(os.getcwd() + "/ConvDiv_sites/Background"):
		os.mkdir(os.getcwd() + "/ConvDiv_sites/Background")
	convergentFileOutput = "{0}/convergentMutations_{1}_{2}.txt".format(os.getcwd() + "/ConvDiv_sites", 
		species_string, position_conservation)
	backgroundFileOutput = "{0}/background_{1}_{2}.txt".format(os.getcwd() + "/ConvDiv_sites/Background", 
		species_string, position_conservation)
	divergentFileOutput = "{0}/divergentMutations_{1}_{2}.txt".format(os.getcwd() + "/ConvDiv_sites", 
		species_string, position_conservation)
	return convergentFileOutput, backgroundFileOutput, divergentFileOutput
	
	
















