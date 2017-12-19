#!/usr/bin/python
#
# Amir Marcovitz
# 
# helper functions for mammalianCodingConvergence.py

import sys
import optparse
import gzip
import os
import copy
import pprint
import subprocess
import operator

import BranchLengthScoring as BLS

species_list = sorted(["hg38","ailMel1","bosTau8","calJac3","camFer1","canFam3","capHir1","cavPor3","cerSim1",
		"chiLan1","chlSab2","chrAsi1","conCri1","criGri1","dasNov3","echTel2","eleEdw1","eptFus1","equCab2",
		"eriEur2","felCat8","hetGla2","jacJac1","lepWed1","loxAfr3","macFas5","mesAur1","micOch1",
		"mm10","musFur1","myoDav1","myoLuc2","nomLeu3","ochPri3","octDeg1","odoRosDiv1","orcOrc1",
		"oryAfe1","oryCun2","otoGar3","oviAri3","panHod1","panTro5","papAnu2","ponAbe2","pteAle1","pteVam1","rheMac8",
		"rn6","saiBol1","sorAra2","speTri2","susScr3","triMan1","tupChi1","turTru2","vicPac2"])

legal_amino_acids = "ACDEFGHIKLMNPQRSTVWY"

pp = pprint.PrettyPrinter(indent = 4)

paml_control_template = '''
    seqfile = ./paml4.8/convergence/control/{0}_{1}.aa * sequence data file name
    outfile = ./paml4.8/convergence/output/{0}_{1}.mp        * main result file
   treefile = ./paml4.8/convergence/control/{0}_{1}.trees  * tree structure file name

    seqtype = 2  * 0:nucleotides; 2:amino acids, 3:binary
      ncatG = 8  * # of categories in the dG model of rates
      nhomo = 0  * nonhomogeneous in calcualting P for branch
'''

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
	#with gzip.open('data/filteredTranscripts.hg38.trial.txt.gz') as f:
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
	'''
	  Create output dir (if not exists) and customed output filenames
	'''
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


def species_in_species_list(species):
	'''
	  Making sure the species is in our screen
	'''
	return species in species_list


def convert_alignment_to_sequence(alignment):
	'''
	  Outputs a clean amino-acid sequence
	'''
	sequence = ""
	for l in alignment:
		if l not in "0123":
			sequence += l
	return sequence.upper()


def get_codon_locations(exons):
	'''
	  Get the exon locations of a transcript (on human genome, GRCh38/hg38)
	'''
	codon_locations = []
	# Each exon tuple contains 7 bases of 5' padding and 6 bases
	# of 3' padding on the Watson strand, as displayed in the genome browser.
	if exons[0][3] == "+":
		# Set current exon and position at beginning of gene
		done = False
		current_exon = 0
		start_pos = exons[current_exon][1] + 7
		stop_pos = start_pos + 2
		# Loop through exons, computing codon positions.
		while True:
			# If we're off the end of the exon, move to the next one.
			if start_pos > exons[current_exon][2] - 6:
				offset = start_pos - (exons[current_exon][2] - 6) - 1
				while start_pos > (exons[current_exon][2] - 6):
					current_exon += 1
					# If there is no next exon, then we're done.
					# Return the list of codon locations.
					if current_exon >= len(exons):
						return codon_locations
					start_pos = exons[current_exon][1] + 7 + offset
					stop_pos = start_pos + 2
			# Add codon position to list
			codon_locations.append((exons[0][0], start_pos, stop_pos, current_exon+1))
			# Advance to next codon position
			start_pos += 3
			stop_pos += 3
	else:
		# Set current exon and position at beginning of gene
		done = False
		current_exon = 0
		stop_pos = exons[current_exon][2] - 6
		start_pos = stop_pos - 2
		# Loop through exons, computing codon positions.
		while True:
			# If we're off the end of the exon, move to the next one.
			if stop_pos < exons[current_exon][1] + 7:
				offset = (exons[current_exon][1] + 7) - stop_pos - 1
				while stop_pos < exons[current_exon][1] + 7:
					current_exon += 1
					# If there is no next exon, then we're done.
					# Return the list of codon locations.
					if current_exon >= len(exons):
						return codon_locations
					stop_pos = exons[current_exon][2] - 6 - offset
					start_pos = stop_pos - 2
			# Add codon position to list
			codon_locations.append((exons[0][0], start_pos, stop_pos, current_exon+1))
			# Advance to next codon position
			start_pos -= 3
			stop_pos -= 3
	
	
def present_in_sufficient_species(species_amino_acids, target_groups, outgroups, convergentSoft):
	'''
	  Make sure at least one member of each target_group is present.
	'''
	if convergentSoft:
		for tg in target_groups:
			present = False
			for ts in tg:
				if ts in species_amino_acids:
					present = True
					break
			if not present:
				return False
	else:
	# Make sure all target species are present.
		for tg in target_groups:
			for ts in tg:
				if ts not in species_amino_acids:
					return False
	# Make sure at least one member of each outgroup is present.
	for og in outgroups:
		present = False
		for os in og:
			if os in species_amino_acids:
				present = True
				break
		if not present:
			return False
	return True


def count_amino_acids(species_amino_acids):
	'''
	  Count how many times each amino acid appears at this position in the sequence.
	  Returns a dictionary of counts.
	'''
	AA_counts = {}
	for species in species_amino_acids.keys():
		AA = species_amino_acids[species]
		if AA not in AA_counts:
			AA_counts[AA] = 0
		AA_counts[AA] += 1
	return AA_counts


def get_BBLS_conservation(species_amino_acids, A0, target_groups, workdir):
	'''
	  This function computes a conservation score based on the Bayesian Branch Length scoring method. 
	  This method is accounting for bransh lengths between species and is more agnostic to phylogenetic topology.
	  Output: A score ranging from 0 to 1, showing the level of conservation at an amino-acid position.
	'''
	# Remove target species from list so as not to bias results
	trimmed_AA_index = copy.deepcopy(species_amino_acids)
	for tg in target_groups:
		for ts in tg:
			if ts in trimmed_AA_index:
				del trimmed_AA_index[ts]
	species_to_keep = []
	for species in trimmed_AA_index.keys():
		species_to_keep.append(species)
	# Load and prune tree
	text = os.popen("tree_doctor {0}/data/mammals_hg38.nh -P {1}".format(workdir, ",".join(species_to_keep))).read()
	tree = BLS.parse_tree(text)
	'''
	with open("./paml4.8/convergence.trees") as f:
		for i in range(2):
			f.readline()
		text = f.read().strip()
		tree = BLS.parse_tree(text)
	os.system(mammalian_tree_command)
	cmd = "tree_doctor pipe.tmp -n -P {0}".format(",".join(species_to_keep))
	p1 = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	text = p1.communicate()[0]
	tree = BLS.parse_tree(text)
	'''
	# Construct dictionary showing which species have the
	# canonical amino acid
	leaves = {}
	for key in trimmed_AA_index.keys():
		if trimmed_AA_index[key] == A0:
			leaves[key] = 1.0
		else:
			leaves[key] = 0.0
	# Compute BBLS and max possible BBLS
	bbls = BLS.BBLS(tree, leaves)
	max_bbls = BLS.getMaxBBLS(tree)
	return bbls / max_bbls


def groups_have_same_amino_acid(target_groups, species_amino_acids, convergentSoft):
	'''
	  Input: A list of 2 target groups, each of which is a list of species group of interest.

	  default: convergentSoft=True
		Can tolerate a missing sequence in some of the target species as long as:
			1) the 'convergent' amino acid is present at least once in all the targets
			2) no other amino acids are present in the target species

	  Output: True if above criteria are satisified, otherwise False.
	'''
	if convergentSoft:
		group1_AA=[]
		group2_AA=[]
		for species in target_groups[0]:
			if species in species_amino_acids:
				group1_AA.append(species_amino_acids[species])
		for species in target_groups[1]:
			if species in species_amino_acids:
				group2_AA.append(species_amino_acids[species])
		return list(set(group1_AA)), (list(set(group2_AA))==list(set(group1_AA)) and len(list(set(group1_AA)))==1 and len(list(set(group2_AA)))==1)
	else:
		for group in target_groups:
			for species in group:
				if species not in species_amino_acids or \
					species_amino_acids[species] != species_amino_acids[target_groups[0][0]]:
					return [''], False
		return species_amino_acids[target_groups[0][0]], True


def checkDivergent(target_groups, outgroups, species_amino_acids, A0, convergentSoft):
	'''
	  Check if a tested position is likely a divergent mutation site (different amino acids in the targets)
	  This subroutine makes sure that amino acids are different between the two target groups
	   and that the amino acid of each target is different than its outgroup
	'''
	aAcids1 = [];	aAcids2 = [];	aAcids1_og = [];	aAcids2_og = [];
	missingTargetAlignments = False	
	for s in target_groups[0]:
		if s in species_amino_acids:
			aAcids1.append(species_amino_acids[s])
		else:
			missingTargetAlignments = True
	for s in target_groups[1]:
		if s in species_amino_acids:
			aAcids2.append(species_amino_acids[s])
		else:
			missingTargetAlignments = True
	for s in outgroups[0]:
		if (s in species_amino_acids.keys()):
			aAcids1_og.append(species_amino_acids[s])
	for s in outgroups[1]:
		if (s in species_amino_acids.keys()):
			aAcids2_og.append(species_amino_acids[s])
	# Now, check that each outgroup is different than its target
	success_1 = (len(list(set(aAcids1) & set(aAcids1_og)))==0)
	success_2 = (len(list(set(aAcids2) & set(aAcids2_og)))==0)
	#if (success_1 and success_2 and len(list(set(aAcids1) & set(aAcids2)))==0):
	if convergentSoft:
		if (len(list(set(aAcids1) & set(aAcids2)))==0 and len(list(set(aAcids1) & set(A0)))==0 and len(list(set(aAcids2) & set(A0)))==0 and success_1 and success_2 and len(aAcids1)>0 and len(aAcids2)>0):
			return 1, '|'.join(aAcids1 + aAcids2)
		else:
			return 0, 'X'
	else:
		if (len(list(set(aAcids1) & set(aAcids2)))==0 and len(list(set(aAcids1) & set(A0)))==0 and len(list(set(aAcids2) & set(A0)))==0 and success_1 and success_2 and len(aAcids1)>0 and len(aAcids2)>0 and not missingTargetAlignments):
			return 1, '|'.join(aAcids1 + aAcids2)
		else:
			return 0, 'X'


def outgroups_have_different_amino_acids(target_amino_acid, outgroups, species_amino_acids):
	'''
	  Input: amino acid found in target species, list of outgroups,
	         dict of amino acids at this position in each species.
	  Output: True if for every outgroup, there exists a species
		 that does not have the target amino acid. False otherwise.
	'''
	for og in outgroups:
		success = False
		for species in og:
			if species in species_amino_acids and species_amino_acids[species] != target_amino_acid:
				success = True
				break
		if not success:
			return False
	return True


def trim_newick_tree(species_to_keep, species_string, position_conservation, workdir):
	'''
	  Call tree doctor and prune the mammalian Newick tree.
	'''
	os.system("printf '{0}\\t1\\n\\n' > {1}/paml4.8/convergence/control/{2}_{3}.trees".format(len(species_to_keep), workdir, species_string, position_conservation))
	os.system("tree_doctor ./data/mammals_hg38.nh -P {0} >> {1}/paml4.8/convergence/control/{2}_{3}.trees".format(",".join(species_to_keep), workdir, species_string, position_conservation))


def run_paml_pamp(species_string, position_conservation, workdir):
	'''
	  Function: run_PAML_pamp
	'''
	with open("{0}/paml4.8/convergence/control/{1}_{2}.ctl".format(workdir, species_string, position_conservation), "w") as w:
		w.write(paml_control_template.format(species_string, position_conservation))
	#os.chdir(workdir + "/paml4.8")
	run_command("{0}/paml4.8/pamp {0}/paml4.8/convergence/control/{1}_{2}.ctl".format(workdir, species_string, position_conservation))
	#os.chdir("/cluster/u/amirma/rot/mike/scripts/convergence")


def run_command(command):
	subprocess.check_call(command, shell=True)


def parse_pamp_results(total_aligned, species_string, position_conservation, workdir):
	'''
	  Read results of a PAML pamp run and extract the relevant info

	  Output: parent = dict representation of phylogenetic tree.
		  ancestral_seqs = amino acids inferred at each ancestral node. 
	'''
	parent = {}
	l = 1
	with open("{0}/paml4.8/convergence/output/{1}_{2}.mp".format(workdir, species_string, position_conservation)) as f:
		# First figure out tree structure from file
		line = f.readline()
		l += 1
		print l
		while line != "":
			if not line.startswith("(1) Branch lengths and substitution pattern"):
				line = f.readline()
				continue
			break
		tree = f.readline().strip().split()
		for t in tree:
			data = [int(d) for d in t.split("..")]
			parent[data[1]] = data[0]	
		while line != "":
			if not line.startswith("(3) Parsimony reconstructions"):
				line = f.readline()
				continue
			break
		for i in range(4):
			f.readline()
		line = f.readline()
		seqs = line.split("|")[-1].strip()
		if seqs == "":
			seqs = line.split("|")[0].split(":")[1].strip()
		try:
			confidence = float(seqs.split("(")[1].split(")")[0])
			print(confidence)
			seqs = seqs.split("(")[0].strip()
		except:
			confidence = 1.0
		ancestral_seqs = {}
		issue_detected = False
		for j in range(len(seqs)):
			if seqs[j] not in legal_amino_acids:
				issue_detected = True
				print(seqs[j])
				break
			ancestral_seqs[total_aligned + j + 1] = seqs[j]
		if issue_detected:
			# Throw an error if we have an issue; in that case
			# we'll just move on the next position in sequence.
			return (None, None, None)
	return (parent, ancestral_seqs, confidence)


def shows_convergence(key_ancestors, target_groups, ancestral_seqs, species_amino_acids, convergentSoft):
	'''
	  Input: list of indices of key ancestors, list of species groups in which we
	         are testing for convergence, list of ancestral sequences at internal
	         tree nodes, list of amino acids at each extant species.

	  Output: True if the target groups have converged,
		  otherwise False.
	'''
	# We just want to make sure that the ancestral amino acid is not present in any of the target groups, 
	# but we tolerate missing sequence
	if convergentSoft:
		for i in range(len(key_ancestors)):
			group_AA=[]
			for species in target_groups[i]:
				if species in species_amino_acids:
					group_AA.append(species_amino_acids[species])
			if ancestral_seqs[key_ancestors[i]] in group_AA:
				return False
		return True
	else:
		for i in range(len(key_ancestors)):
			for tg in target_groups[i]:
				if tg in species_amino_acids:
					if ancestral_seqs[key_ancestors[i]] == species_amino_acids[target_groups[i][0]]:
						return False
				else:
					return False
		return True


def get_conservation_score(sequences, position, padding):
	# Sanity check
	for seq in sequences:
		assert len(seq) == len(sequences[0])
	
	# Trim off the parts of the sequences that we actually care about
	cut_seqs = []
	for seq in sequences:
		# The extra math here ensures that we don't go beyond the edges of the sequence.
		cs = seq[max(0,position-padding):position] + seq[position+1:min(position+padding+1, len(sequences[0]))]
		cut_seqs.append(cs)
	consensus = compute_consensus_sequence(cut_seqs)
	# Compute score
	score = 0
	max_score = 0
	for i in range(len(cut_seqs[0])):
		for cs in cut_seqs:
			max_score += 1
			if cs[i] == consensus[i] and cs[i] not in "?-":
				score += 1
	return score*1.0/max_score


def compute_consensus_sequence(sequences):
	'''
	  Input: a list of sequences of the same length
	
	  Output: a single sequence showing the character that appears
		  most frequently at each position within the input sequences.
	'''
	# Sanity check
	for seq in sequences:
		assert len(seq) == len(sequences[0])
	
	# Build consensus sequence
	consensus = ""
	for i in range(len(sequences[0])):
		dCounts = {}
		for seq in sequences:
			dCounts[seq[i]] = dCounts.get(seq[i], 0) + 1
		consensus += sorted(dCounts.items(), reverse = True, key = operator.itemgetter(1))[0][0]

	return consensus
















