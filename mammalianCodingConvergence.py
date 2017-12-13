#!/usr/bin/python
#
# Amir Marcovitz
# 
# code base for:
#	"A novel unbiased test for molecular convergent evolution and discoveries in echolocating, aquatic and high-altitude mammals"
#	Amir Marcovitz, Yatish Turakhia, Michael Gloudemans, Benjamin A Braun, Heidi I Chen & Gill Bejerano
#	https://doi.org/10.1101/170985

"""
mammalianCodingConvergence - A code for coding convergent and divergent mutations scanning in mammalian genomes 

Usage:
 mammalianCodingConvergenc.py --tg {target_group1}:{target_group2} --og {outroup1}:{outgroup2} --pc 0.9
  where
  --tg	target_group    	is a comma-delimited list of target1 or 2 species group (can provide also genome assembly abbreviation)
  --og	outgroup        	is a comma-delimited list of outgroup1 or 2 species group (can provide also genome assembly abbreviation)
  --pc	position conservarion	is the Bayesian Branch Length Conservation Score (BBLS) of an amino acid position 

Example:
 mammalianCodingConvergence.py --tg bosTau7,capHir1,oviAri1,panHod1:pteAle1,pteVam1 --og orcOrc1,turTru2:eptFus1,myoDav1,myoLuc2 --pc 0.9
"""

import ConvUtil as conv
import os

legal_amino_acids = "ACDEFGHIKLMNPQRSTVWY"	# 20 amino acids
minimum_species_cutoff = 40			# Amino acid must be present in AT LEAST this many species for us to consider it

def run(target_groups, outgroups, position_conservation, target_species, species_string):

	transcript_list, transcript_to_gene, exon_locations, gene_symbols = conv.load_RefGenes()
	convergentFileOutput, backgroundFileOutput, divergentFileOutput = conv.init_outputs(species_string, position_conservation)

	# Main loop
	iter = -1
	with open(convergentFileOutput, "w") as w:
		with open(backgroundFileOutput, "w") as wTestCounts:
			wd = open(divergentFileOutput, "w")
			total_positions = 0		# Total positions at which humans have an amino acid.
			total_positions_tested = 0	# Total positions at which we applied convergence/divergence tests.
			total_converged_positions = 0	# Count of of total convergent positions
			for transcript in transcript_list:
				iter += 1
				bin = int(transcript[-2:])
				with gzip.open("{0}/{1}/{2}.txt.gz".format(os.getcwd() + "/data/protein_alignments_hg38", bin, transcript)) as alignment_file:
					sequence_list = {}
					data = alignment_file.readline().strip().split()
					assert transcript == data[0]
					sequence_list["hg38"] = data[1]
					species_seen = set([])
					for line in alignment_file:
						data = line.strip().split()
						if line[0]=="#":
							continue
						assert transcript == data[1]
						species = data[0]
						if not conv.species_in_species_list(species):
							continue
						# Discard species appearing multiple times
						if species in species_seen:
							try:
								del sequence_list[species]
							except:
								pass
							continue
						species_seen.add(species)
						sequence_list[species] = conv.convert_alignment_to_sequence(data[2])
						assert len(sequence_list[species]) == len(sequence_list["hg38"])
					codon_locations = conv.get_codon_locations(exon_locations[transcript])
					# Making sure that the number of codons in list should always be equivalent to number
					# of positions in amino acid sequence.
					if len(codon_locations) != len(sequence_list["hg38"]):
						continue
					number_tested = 0	# A count of how many positions we test in the gene (i.e., in the current transcript)
					# Now loop through the sequence and analyze one position at a time.
					for i in range(len(sequence_list["hg38"])):
						divergent = 0	# Divergent substitution
						# skeeeping over stop codons (marked as '*')
						if sequence_list["hg38"][i] == "*":
							continue
						# Human sequences should only contain amino acids. If not, quit and alert.
						if sequence_list["hg38"][i] not in legal_amino_acids:
							print("Error in sequence", sequence_list["hg38"], "at position", i)
							assert False
						total_positions += 1
						# Make map of amino acids found in every species.
						# Delete species with invalid amino acids.
						species_amino_acids = {}
						for species in species_list:
							if species in sequence_list:
								if sequence_list[species][i] in legal_amino_acids:
									species_amino_acids[species] = sequence_list[species][i]
						# If not enough species aligned here, then just move on to next AA.
						total_aligned = len(species_amino_acids.items())					
						if total_aligned < minimum_species_cutoff:
							continue
						# Accept position if at least one species from each target is represented (default: convergentSoft = True)
						if not conv.present_in_sufficient_species(species_amino_acids, target_groups, outgroups, convergentSoft=True):
							continue
						AA_counts = conv.count_amino_acids(species_amino_acids)
						AA_sort = sorted(AA_counts.items(), reverse = True, key = operator.itemgetter(1))
						A0, nA0 = AA_sort[0][0], AA_sort[0][1]
						# Compute amino-acid position conservation with BBLS, and make sure it's passing the threshold
						BBLS_conservation = conv.get_BBLS_conservation(species_amino_acids, A0, target_groups)
						if BBLS_conservation < position_conservation:
							continue						 
						number_tested += 1
						total_positions_tested += 1
						# Check to see whether our target species have the same amino acid at this position.
						# If so, then proceed with analysis. If not, check divergent substitution (default: convergentSoft = True)
						target_AA, same_AA = conv.groups_have_same_amino_acid(target_groups, species_amino_acids, convergentSoft=True)
						if not same_AA:
							divergent, divergentAA = conv.checkDivergent(target_groups, 
								outgroups, species_amino_acids, A0, convergentSoft=True)
							if not divergent:
								continue
						# If outgroups have identical amino acids to our target species, we will not call it
						# convergence/divergence. PAML analysis will be skipped and save some run time.
						if not conv.outgroups_have_different_amino_acids(target_AA[0], outgroups, species_amino_acids):
							continue
						# Now, infer ancestral sequences and determine whether convergent evolution has occurred.		
						aligned_species = species_amino_acids.keys()
						species_indices = {}
						with open("/cluster/u/amirma/rot/mike/bin/paml4.8/convergence/control/{0}_{1}.aa".format(species_string, position_conservation) ,"w") as wConv:
							wConv.write("{0} 1\n\n".format(len(species_amino_acids.keys())))
							index = 1
							for species in species_amino_acids.keys():
								wConv.write("{0}   {1}\n".format(species, species_amino_acids[species]))
								species_indices[species] = index
								index += 1
						# Use tree_doctor to create the corresponding Newick tree for analysis with PAML pamp
						conv.trim_newick_tree(species_amino_acids.keys(),species_string, position_conservation)
						conv.run_paml_pamp(species_string, position_conservation)
						parent, ancestral_seqs, confidence = conv.parse_pamp_results(total_aligned, species_string, position_conservation)
						if parent == None:	
							number_tested -= 1
							total_positions_tested -= 1
							continue
						# Determine which is outgroup branchpoint for each target group.
						key_ancestors = []
						for index in range(len(target_groups)):
							group = target_groups[index] + outgroups[index]
							ancestor_lines = []
							for species in group:
								# We have to do this check because not all outgroup species
								# need to be represented.
								if species not in species_indices:
									continue
								current = species_indices[species]
								line = [current]
								while current in parent:
									line.append(parent[current])
									current = parent[current]
								ancestor_lines.append(line)
							for j in range(len(ancestor_lines[0])):
								present_in_all = True
								for al in ancestor_lines:
									if ancestor_lines[0][j] not in al:
										present_in_all = False
								if present_in_all:
									key_ancestors.append(ancestor_lines[0][j+1])
									break
						if not conv.shows_convergence(key_ancestors, target_groups, ancestral_seqs, species_amino_acids, convergentSoft=True):
							continue
						if divergent==0:
							total_converged_positions += 1
						# Make a list of all species besides the target species that also have substitutions in their sequences
						nonconforming_species = []
						nonconforming_AAs = []
						conforming_species = []
						for species in species_list:
							if species in sequence_list and sequence_list[species][i] != A0 and \
							 species not in target_species and sequence_list[species][i] not in "-?":
								nonconforming_species.append(species)
								nonconforming_AAs.append(sequence_list[species][i])
							else:
								conforming_species.append(species)
						# Make sure we have at least something to output in every case, so we don't break the file format.
						nonconforming_count = len(nonconforming_species)
						if nonconforming_count == 0:
							nonconforming_species = ["NA"]
							nonconforming_AAs = ["NA"]
						gene_symbol = gene_symbols[transcript_to_gene[transcript]]
						if not divergent:
							w.write("\t".join([transcript_to_gene[transcript],
								  transcript,
								  gene_symbol,
								  codon_locations[i][0],
								  str(codon_locations[i][1]),
								  str(codon_locations[i][2]),
								  str(codon_locations[i][3]),
								  str(i),
								  str(len(sequence_list["hg38"])),
								  str(A0),
								  str(nA0),
								  str(total_aligned),
								  str(nA0 * 1.0 / total_aligned),
								  "{0}".format("|".join(target_species)),
								  target_AA[0],
								  "|".join(nonconforming_species),
								  "|".join(nonconforming_AAs),
								  "|".join(conforming_species),
								  str(nonconforming_count),
								  str(BBLS_conservation),
								  str(conservation_window_padding),
								  str(get_conservation_score(sequence_list.values(), i, 2)),
								  "2",
								  str(get_conservation_score(sequence_list.values(), i, 10)),
								  "10",
								  str(confidence)
								  ])
									+ "\n")
						else:
							wd.write("\t".join([transcript_to_gene[transcript],
								  transcript,
								  gene_symbol,
								  codon_locations[i][0],
								  str(codon_locations[i][1]),
								  str(codon_locations[i][2]),
								  str(codon_locations[i][3]),
								  str(i),
								  str(len(sequence_list["hg38"])),
								  str(A0),
								  str(nA0),
								  str(total_aligned),
								  str(nA0 * 1.0 / total_aligned),
								  "{0}".format("|".join(target_species)),
								  divergentAA,
								  "|".join(nonconforming_species),
								  "|".join(nonconforming_AAs),
								  "|".join(conforming_species),
								  str(nonconforming_count),
								  str(BBLS_conservation),
								  str(conservation_window_padding),
								  str(get_conservation_score(sequence_list.values(), i, 2)),
								  "2",
								  str(get_conservation_score(sequence_list.values(), i, 10)),
								  "10",
								  str(confidence)
								  ])
									+ "\n")


if __name__ == "__main__":
	target_groups, outgroups, position_conservation, target_species, species_string = conv.parse()
	run(target_groups, outgroups, position_conservation, target_species, species_string)
