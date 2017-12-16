#!/bin/bash -e
#
# Amir Marcovitz
# 
# a script for computing pathway enriched in convergent and/or divergent mutations
#

ontoFile=${PWD}/data/GeneSymbolAnnotations.MGI2016.NoSyn.onto.gz
ConvDivDir=${PWD}/ConvDiv_sites

if [[ "$#" -ne 3 ]]; then
  echo -e "\nEnrichment analysis of a convergent evolution experiment and its controls\n  Usage:\n  $0 target_group_1 target_group_2 position_conservation\n"
  echo -e "  Example:\n  $0  dolphin,killer_whale  microbat,big_brown_bat,davids_myotis  0.9\n"
  exit 1
fi

if [ ! -d ${PWD}/enrichments ]; then
	mkdir ${PWD}/enrichments
fi
enrich_dir=${PWD}/enrichments


## I - parsing the inputs:
target1=$1
target2=$2
position_conservation=$3
target1_assembly=`grep -f <(echo -e ${target1} | sed -e "s/,/\n/g") \
	<(zcat ./data/species_list.txt.gz) | cut -f2 | sort -u | awk '{printf $0","} END {print ""}'`
target2_assembly=`grep -f <(echo -e ${target2} | sed -e "s/,/\n/g") \
	<(zcat ./data/species_list.txt.gz) | cut -f2 | sort -u | awk '{printf $0","} END {print ""}'`
experiment=`cat <(echo $target1_assembly) <(echo $target2_assembly) | sort -u | sed -e "s/,/_/g" \
	| awk -v position_conservation=${position_conservation} '{printf $0} END {print position_conservation}'`
echo -e $experiment


## II - look on genes and amino-acids that are in our background (conservation threshold dependant)
join -t$'\t' -1 3 -2 1 -o '1.1 1.2 1.3 1.4 2.1 2.2 2.3 2.4' <(zcat ${PWD}/ConvDiv_sites/Background/background_${experiment}.txt.gz | sort -t$'\t' -k3,3) \
	<(zcat $ontoFile | sort -u) | awk -F'\t' '{if($4>0) print $5"\t"$6"\t"$7"\t"$8"\t"$4}' | sort -u > ontoFileTrimmed.txt
cut -f2,3 ontoFileTrimmed.txt | sort -u > ontoTerms
echo -e "termID\ttermDef\ttotalTermGenes\ttotalTermSizeAA\tnConvGenes\tnConvSites\tnDivGenes\tnDivSites\ttotalMutatedGenes\tfuncClust" \
	> ${enrich_dir}/ConvDiv_${experiment}.txt


## III - Iterate over all pathways (i.e., MGI-phenotype ontology terms) to collect convergent and divergent events
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

#experimentsFile=$1
#consThres=$2
#Onto=$3
#pheno=`echo $experimentsFile | cut -d"_" -f2 | sed -e "s/.txt//g"`
#funcHighlight=`echo $pheno | sed -e 's/echolocation/hearing/g;s/moles2/vision/g;s/moles3/vision/g;s/aquatic/skin/g;s/StarNoseMoleWeddelSeal/respiratory/g;s/moles/vision/g;s/NakedMoleWeddelSeal/respiratory/g;s/GoldenMoleWeddelSeal/#respiratory/g;s/BipedalsHumanJerboa/skeleton/g;s/ManateeSheepGoat/vision/g;s/Tusk/craniofacial/g;s/ManateeCat/vision/g;s/SheepGoatCat/vision/g;s/halt/respiratory/g;s/hibernation/respiratory/g'`
#grep -v "^#" $experimentsFile > experimentList_${pheno}.txt	# removing commented experiments
#experimentsFile=experimentList_${pheno}.txt
#experimentsDir=${pheno}'_ConvDiv_'${consThres}'_'${Onto}
#ontoFile=$GENELOSS/data/Ontologies_hg38/GeneSymbolAnnotations.${Onto}2016.NoSyn.onto
#if [ $onto2016 -eq 0 ];then
#	#experimentsDir=$pheno'_2_pre2016'
#	ontoFile=$GENELOSS/data/Ontologies_hg38/GeneSymbolAnnotations.${Onto}.NoSyn.onto
#fi
#n=`wc -l $experimentsFile | cut -d" " -f1`	                # the number of experiments (target + control)

#if [ $convergentSoft -ne 0 ];then
#	convergentSoftExtension='.convergentSoft'
#	experimentsDir=${experimentsDir}_convergentSoft
#fi

#if [ ! -d $experimentsDir ]; then
#	mkdir ${experimentsDir}
#fi
##########################

# Iterate over the experiments:


	# filtering on the ontology to look at genes and amino acids that are in our background
	#join -t$'\t' -1 3 -2 1 -o '1.1 1.2 1.3 1.4 2.1 2.2 2.3 2.4' <(sort -t$'\t' -k3,3 $backgroundDir/background_${currExperiment}0.0_${consThres}.txt${convergentSoftExtension}) <(sort -u $ontoFile) | awk -F'\t' '{if($4>0) print $5"\t"$6"\t"$7"\t"$8"\t"$4}' | sort -u > tmpTerms; mv tmpTerms ${Onto}_ontoFileTrimmed.txt
	# trim the ontology file by min,max genes per terms:
	#join -t$'\t' -1 1 -2 2 -o '2.1 2.2 2.3 2.4 2.5' <(cut -f2 ${Onto}_ontoFileTrimmed.txt | sort | uniq -c | awk -v minGenesPerTerm=$minGenesPerTerm -v maxGenesPerTerm=$maxGenesPerTerm '{if($1>=minGenesPerTerm && $1<=maxGenesPerTerm) print $2}' | sort -u) <(sort -t$'\t' -k2,2 ${Onto}_ontoFileTrimmed.txt) | awk -F'\t' '{if($2!=$5) print $0}' | sort -u > tmpTerms; mv tmpTerms ${Onto}_ontoFileTrimmed.txt
	#cut -f2,3 ${Onto}_ontoFileTrimmed.txt | sort -u > ontoTerms
	#echo -e "termID\ttermDef\ttotalTermGenes\ttotalTermSizeAA\tnConvGenes\tnConvSites\tnDivGenes\tnDivSites\ttotalMutatedGenes\tfuncClust" > ${experimentsDir}/ConvDiv_${currExperiment}${consThres}
#	while read termDef
#	do
#		term=`echo -e $termDef | cut -d" " -f1`
#		termText=`grep $term ontoTerms | cut -f2`
#		grep -w $term ${Onto}_ontoFileTrimmed.txt | sort -u > tmpTerm
#		totalTermGenes=`cut -f1 tmpTerm | sort -u | wc -l`
#		totalTermSizeAA=`sort -u tmpTerm | awk -F'\t' '{sum+=$5} END {print sum}'`
#		funcClust=`grep $term ontoFunctionClustered/*.manualFunctionalClustering | cut -f3`
#		nConvSites=`join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.10 2.15 1.1 1.4' <(cat tmpTerm) <(sort -t$'\t' -k3,3 $convergentDir/convergentMutations_${currExperiment}0.0_${consThres}.txt${convergentSoftExtension} | egrep -v "^#") | cut -f4-6 | sort -u | wc -l`

#		nConvGenes=`join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.10 2.15 1.1 1.4' <(cat tmpTerm) <(sort -t$'\t' -k3,3 $convergentDir/convergentMutations_${currExperiment}0.0_${consThres}.txt${convergentSoftExtension} | egrep -v "^#") | cut -f3 | sort -u | wc -l`
#		nDivSites=`join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.10 2.15 1.1 1.4' <(cat tmpTerm) <(sort -t$'\t' -k3,3 $divergentDir/divergentMutations_${currExperiment}0.0_${consThres}.txt${convergentSoftExtension} | egrep -v "^#") | cut -f4-6 | sort -u | wc -l`
#		nDivGenes=`join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.10 2.15 1.1 1.4' <(cat tmpTerm) <(sort -t$'\t' -k3,3 $divergentDir/divergentMutations_${currExperiment}0.0_${consThres}.txt${convergentSoftExtension} | egrep -v "^#") | cut -f3 | sort -u | wc -l`
#		totalMutatedGenes=`cat <(join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.10 2.15 1.1 1.4' <(cat tmpTerm) <(sort -t$'\t' -k3,3 $divergentDir/divergentMutations_${currExperiment}0.0_${consThres}.txt${convergentSoftExtension} | egrep -v "^#") | cut -f3) <(join -t$'\t' -1 1 -2 3 -o '2.1 2.2 2.3 2.4 2.5 2.6 2.10 2.15 1.1 1.4' <(cat tmpTerm) <(sort -t$'\t' -k3,3 $convergentDir/convergentMutations_${currExperiment}0.0_${consThres}.txt${convergentSoftExtension} | egrep -v "^#") | cut -f3) | sort -u | wc -l`
#		echo -e "$term\t$termText\t$totalTermGenes\t$totalTermSizeAA\t$nConvGenes\t$nConvSites\t$nDivGenes\t$nDivSites\t$totalMutatedGenes\t$funcClust" >> ${experimentsDir}/ConvDiv_${currExperiment}${consThres}
#		rm -rf tmpTerm
#	done < ontoTerms
#	R_input="filename <- '${experimentsDir}\/ConvDiv_${currExperiment}${consThres}'"
#	R_func="funcHighlight <- '${funcHighlight}'"
#	R_pheno="pheno <- '${pheno}'"
#	R_background="backgroundFILE <- '\/cluster\/u\/amirma\/git\/forwardGenomics\/geneLoss\/convergentEvolution\/convergent_mutations_hg38\/background\/background_${currExperiment}0.0_${consThres}.txt${convergentSoftExtension}'"
#	R_convergent="convergentFILE <- '\/cluster\/u\/amirma\/git\/forwardGenomics\/geneLoss\/convergentEvolution\/convergent_mutations_hg38\/convergentMutations_${currExperiment}0.0_${consThres}.txt${convergentSoftExtension}'"
#	R_divergent="divergentFILE <- '\/cluster\/u\/amirma\/git\/forwardGenomics\/geneLoss\/convergentEvolution\/divergent_mutations_hg38\/divergentMutations_${currExperiment}0.0_${consThres}.txt${convergentSoftExtension}'"
#	sed -e "s/filename <- FILENAME/${R_input}/g;s/funcHighlight <- FUNCHIGHLIGHT/${R_func}/g;s/pheno <- PHENO/${R_pheno}/g;s/backgroundFILE <- BACKGROUNDFILE/${R_background}/g;s/convergentFILE <- CONVERGENTFILE/${R_convergent}/g;s/divergentFILE <- DIVERGENTFILE/${R_divergent}/g;s/CONSERVATIONTHRESHOLD/${consThres}/g;s/ONTOLOGYNAME/${Onto}/g" ontoTerms_ConvergentDivergentPlots.R > tmpR; mv tmpR ontoTerms_ConvergentDivergentPlots.R
#	Rscript ontoTerms_ConvergentDivergentPlots.R
#	sed -e "s/${R_input}/filename <- FILENAME/g;s/${R_func}/funcHighlight <- FUNCHIGHLIGHT/g;s/${R_pheno}/pheno <- PHENO/g;s/${R_background}/backgroundFILE <- BACKGROUNDFILE/g;s/${R_convergent}/convergentFILE <- CONVERGENTFILE/g;s/${R_divergent}/divergentFILE <- DIVERGENTFILE/g;s/consThres <- ${consThres}/consThres <- CONSERVATIONTHRESHOLD/g;s/${Onto}/ONTOLOGYNAME/g" ontoTerms_ConvergentDivergentPlots.R > tmpR; mv tmpR ontoTerms_ConvergentDivergentPlots.R
#done < $experimentsFile 

#rm -rf ontoTerms ${Onto}_ontoFileTrimmed.txt
