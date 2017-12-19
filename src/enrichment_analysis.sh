#!/bin/bash -e
#
# Amir Marcovitz
# 
# a script for computing pathway enriched in convergent and/or divergent mutations
#

bold=$(tput bold)
normal=$(tput sgr0)

ontoFile=${PWD}/data/GeneSymbolAnnotations.MGI2016.NoSyn.onto.gz
ConvDivDir=${PWD}/ConvDiv_sites
srcDir=${PWD}/src

if [[ "$#" -ne 3 ]]; then
  echo -e "\nEnrichment analysis of a convergent evolution experiment and its controls\n  Usage:\n  $0 target_group_1 target_group_2 position_conservation\n"
  echo -e "  Example:\n  $0  dolphin,killer_whale  microbat,big_brown_bat,davids_myotis  0.9\n"
  exit 1
fi

if [ ! -d ${PWD}/enrichments ]; then
	mkdir ${PWD}/enrichments
fi
enrich_dir=${PWD}/enrichments


## I - Parsing the inputs:
target1=$1
target2=$2
position_conservation=$3
target1_assembly=`grep -f <(echo -e ${target1} | sed -e "s/,/\n/g") \
	<(zcat ./data/species_list.txt.gz) | cut -f2 | sort -u | awk '{printf $0","} END {print ""}'`
target2_assembly=`grep -f <(echo -e ${target2} | sed -e "s/,/\n/g") \
	<(zcat ./data/species_list.txt.gz) | cut -f2 | sort -u | awk '{printf $0","} END {print ""}'`
experiment=`cat <(echo $target1_assembly) <(echo $target2_assembly) | sort -u | sed -e "s/,/_/g" \
	| awk -v position_conservation=${position_conservation} '{printf $0} END {print position_conservation}'`
echo -e "\n\t...running enrichment analysis for ${bold}${experiment}${normal}"


## II - Look on genes and amino-acids that are in our background (conservation threshold dependant)
join -t$'\t' -1 3 -2 1 -o '1.1 1.2 1.3 1.4 2.1 2.2 2.3 2.4' <(zcat ${PWD}/ConvDiv_sites/Background/background_${experiment}.txt.gz | sort -t$'\t' -k3,3) \
	<(zcat $ontoFile | sort -u) | awk -F'\t' '{if($4>0) print $5"\t"$6"\t"$7"\t"$8"\t"$4}' | sort -u > ontoFileTrimmed.txt
cut -f2,3 ontoFileTrimmed.txt | sort -u > ontoTerms


## III - Iterate over all pathways (i.e., MGI-phenotype ontology terms) to collect convergent and divergent events
if [ ! -f ${enrich_dir}/ConvDiv_${experiment}.txt ]; then
	echo -e "\n\t...counting convergent and divergent substitutions per pathway"
	echo -e "termID\ttermDef\ttotalTermGenes\ttotalTermSizeAA\tnConvGenes\tnConvSites\tnDivGenes\tnDivSites\ttotalMutatedGenes\tfuncClust" \
	> ${enrich_dir}/ConvDiv_${experiment}.txt
	while read termDef
	do
		term=`echo -e $termDef | cut -d" " -f1`
		termText=`grep $term ontoTerms | cut -f2`
		grep -w $term ontoFileTrimmed.txt | sort -u > tmpTerm
		totalTermGenes=`cut -f1 tmpTerm | sort -u | wc -l`
		totalTermSizeAA=`sort -u tmpTerm | awk -F'\t' '{sum+=$5} END {print sum}'`
		funcClust=`zcat ${PWD}/data/*.manualFunctionalClustering.gz | grep $term  | cut -f3`
		nConvSites=`join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.10 2.15 1.1 1.4' <(cat tmpTerm) \
			<(zcat ${ConvDivDir}/convergentMutations_${experiment}.txt.gz | egrep -v "^#" | sort -t$'\t' -k3,3) | cut -f4-6 | sort -u | wc -l`
		nConvGenes=`join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.10 2.15 1.1 1.4' <(cat tmpTerm) \
			<(zcat ${ConvDivDir}/convergentMutations_${experiment}.txt.gz | egrep -v "^#" | sort -t$'\t' -k3,3) | cut -f3 | sort -u | wc -l`
		nDivSites=`join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.10 2.15 1.1 1.4' <(cat tmpTerm) \
			<(zcat ${ConvDivDir}/divergentMutations_${experiment}.txt.gz | egrep -v "^#" | sort -t$'\t' -k3,3) | cut -f4-6 | sort -u | wc -l`
		nDivGenes=`join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.10 2.15 1.1 1.4' <(cat tmpTerm) \
			<(zcat ${ConvDivDir}/divergentMutations_${experiment}.txt.gz | egrep -v "^#" | sort -t$'\t' -k3,3) | cut -f3 | sort -u | wc -l`
		totalMutatedGenes=`cat <(join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.10 2.15 1.1 1.4' \
			<(cat tmpTerm) <(zcat ${ConvDivDir}/divergentMutations_${experiment}.txt.gz | egrep -v "^#" | sort -t$'\t' -k3,3) | cut -f3) \
			<(join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.10 2.15 1.1 1.4' <(cat tmpTerm) \
			<(zcat ${ConvDivDir}/convergentMutations_${experiment}.txt.gz | egrep -v "^#" | sort -t$'\t' -k3,3) | cut -f3) | sort -u | wc -l`
		echo -e "$term\t$termText\t$totalTermGenes\t$totalTermSizeAA\t$nConvGenes\t$nConvSites\t$nDivGenes\t$nDivSites\t$totalMutatedGenes\t$funcClust" \
			>> ${enrich_dir}/ConvDiv_${experiment}.txt
		rm -rf tmpTerm
	done < ontoTerms
fi


## IV - Compute enrichments (of convergent and divergent mutations) per term (R)
cat <(echo -e "filename <- \""${enrich_dir}"/ConvDiv_"${experiment}".txt\"") \
	<(echo -e "ConvDiv <- read.table(filename, sep='\\\t', quote='', header=TRUE)") \
	<(echo -e "backgroundFILE <- \""${PWD}"/ConvDiv_sites/Background/background_"${experiment}".txt.gz\"") \
	<(echo -e "convergentFILE <- \""${PWD}"/ConvDiv_sites/convergentMutations_"${experiment}".txt.gz\"") \
	<(echo -e "divergentFILE <- \""${PWD}"/ConvDiv_sites/divergentMutations_"${experiment}".txt.gz\"\n") \
	<(echo -e "consThres <- "${position_conservation}"\n") \
	<(cat ${srcDir}/ontoTerms_ConvergentDivergentPlots.R) > ${PWD}/ontoTerms_ConvergentDivergentPlots.R

Rscript ${PWD}/ontoTerms_ConvergentDivergentPlots.R

rm -rf ${enrich_dir}/ConvDiv_${experiment}.txt

## V - Output enrichment results, if found
echo -e "\n   ${bold}Results:${normal}\n   ${bold}--------${normal}"
n_conv_enrichments=`awk -F'\t' '{if($25<0.05 && $17>2) print $0}' ${enrich_dir}/ConvDiv_${experiment}.termEnrichments | wc -l`
if [ ${n_conv_enrichments} -ne 0 ]; then
	n_conv_enrichments=`awk -F'\t' '{if($25<0.05 && $17>2 && $26>0.05 && $18<2) print $0}' ${enrich_dir}/ConvDiv_${experiment}.termEnrichments | wc -l`
	if [ ${n_conv_enrichments} -ne 0 ]; then
		echo -e "\n\t${bold}Convergent Pathway Found${normal}\n"
		paste <(echo -e "${bold}\tontology-term\n\t#Conv. Genes\n\t#Conv. Sites\n\t#Conv. q-value\n\t#Conv. Fold\n\t#Div. Genes\n\t#Div. Sites\n\t#Div. q-value\n\t#Div. Fold${normal}") <(awk -F'\t' '{if($25<0.05 && $17>2 && $26>0.05 && $18<2) print $0}' ${enrich_dir}/ConvDiv_${experiment}.termEnrichments | sort -t$'\t' -k 25,25g -k17rn,17 | head -1 | awk -F'\t' '{print "\t"$2" ("$1")\n\t"$5"\n\t"$6"\n\t"$25"\n\t"$17"\n\t"$7"\n\t"$8"\n\t"$26"\n\t"$18}') 
	else
		echo -e "\n\t${bold}Relaxation Pathway Found${normal}\n"
		paste <(echo -e "${bold}\tontology-term\n\t#Conv. Genes\n\t#Conv. Sites\n\t#Conv. q-value\n\t#Conv. Fold\n\t#Div. Genes\n\t#Div. Sites\n\t#Div. q-value\n\t#Div. Fold${normal}") <(awk -F'\t' '{if($25<0.05 && $17>2) print $0}' ${enrich_dir}/ConvDiv_${experiment}.termEnrichments | sort -t$'\t' -k 25,25g -k17rn,17 | head -1 | awk -F'\t' '{print "\t"$2" ("$1")\n\t"$5"\n\t"$6"\n\t"$25"\n\t"$17"\n\t"$7"\n\t"$8"\n\t"$26"\n\t"$18}') 
	fi
else
	echo -e "\n\tNo enrichments were found in ${experiment}\n\n"
fi


echo -e "\n\n"
