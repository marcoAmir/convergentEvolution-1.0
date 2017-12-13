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


def run(target_groups, outgroups, position_conservation, target_species, species_string):

	transcript_list, transcript_to_gene, exon_locations, gene_symbols = conv.load_RefGenes()
	convergentFileOutput, backgroundFileOutput, divergentFileOutput = conv.init_outputs(species_string, position_conservation)

	# Main loop
	iter = -1
	with open(convergentFileOutput, "w") as w:
		with open(backgroundFileOutput, "w") as wTestCounts:
			wd = open(divergentFileOutput, "w")
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
						
				
				
	

if __name__ == "__main__":
	target_groups, outgroups, position_conservation, target_species, species_string = conv.parse()
	run(target_groups, outgroups, position_conservation, target_species, species_string)
