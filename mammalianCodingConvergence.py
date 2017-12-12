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
				
	

if __name__ == "__main__":
	if not os.path.isdir(os.getcwd() + "/ConvDiv_sites"):
		os.mkdir(os.getcwd() + "/ConvDiv_sites")
	if not os.path.isdir(os.getcwd() + "/ConvDiv_sites/Background"):
		os.mkdir(os.getcwd() + "/ConvDiv_sites/Background")
	target_groups, outgroups, position_conservation, target_species, species_string = conv.parse()
	run(target_groups, outgroups, position_conservation, target_species, species_string)
