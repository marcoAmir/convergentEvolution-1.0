require(ggplot2)

filename <- CURR_FILENAME
ConvDiv <- read.table(filename, sep='\t', quote='', header=TRUE)
backgroundFILE <- BACKGROUNDFILE
convergentFILE <- CONVERGENTFILE
divergentFILE <- DIVERGENTFILE


applyOntoDagReduce <- TRUE

# upload the background, convergent and divergent data, and ontology file:
background <- read.table(backgroundFILE, col.names=c('ensgID','enstID','geneName','nPositions'))
convergent <- read.table(convergentFILE, col.names=c('ensgID','enstID','geneName','chr','codonStart',
							'codonEnd','exon','protPos','protLen','A0','nA0',
							'speciesAligned','proportionA0','targetSpecies','targetAA',
							'nonconforming_species','nonconforming_AAs','conforming_species',
							'nonconforming_count','BBLS_conservation','conservation_window_padding',
							'conservation_score','X1','X2','X3','confidence)'))
divergent <- read.table(divergentFILE, col.names=c('ensgID','enstID','geneName','chr','codonStart',
							'codonEnd','exon','protPos','protLen','A0','nA0',
							'speciesAligned','proportionA0','targetSpecies','targetAA',
							'nonconforming_species','nonconforming_AAs','conforming_species',
							'nonconforming_count','BBLS_conservation','conservation_window_padding',
							'conservation_score','X1','X2','X3','confidence)'))
ontoFile <- read.table('ontoFileTrimmed.txt', sep='\t', quote='', col.names=c('geneName','termID','termDef','ensgID','nPositions'))

# trim onto by genes/term:
genePerTerm <- data.frame(table(ontoFile$termID))
colnames(genePerTerm) <- c('termID','nGenes')

funcHighlight <- FUNCHIGHLIGHT
pheno <- PHENO
consThres <- CONSERVATIONTHRESHOLD

# Some normalizations (not used in current study or enrichments)
### Normalization1 - substituion per length ('length' of onto-term is the total AA length of its genes)
ConvDiv$normedConvergence <- (ConvDiv$nConvSites/ConvDiv$totalTermSizeAA)
ConvDiv$normedDivergence <- (ConvDiv$nDivSites/ConvDiv$totalTermSizeAA)
### Normalization2 - mutated genes per term size (term size = its number of genes)
ConvDiv$normedConvergenceGenes <- (ConvDiv$nConvGenes/ConvDiv$totalTermGenes)
ConvDiv$normedDivergenceGenes <- (ConvDiv$nDivGenes/ConvDiv$totalTermGenes)	
### Substitutions per gene
ConvDiv$convPerGene <- ConvDiv$nConvSites/ConvDiv$nConvGenes
ConvDiv$divPerGene <- ConvDiv$nDivSites/ConvDiv$nDivGenes
ConvDiv$convPerGene[is.nan(ConvDiv$convPerGene)] <- 0
ConvDiv$divPerGene[is.nan(ConvDiv$divPerGene)] <- 0

# compute Fold enrichments
ConvDiv$FoldChangeConv <- rep(1, nrow(ConvDiv))
ConvDiv$FoldChangeDiv <- rep(1, nrow(ConvDiv))
for(i in 1:nrow(ConvDiv)){
	ConvDiv$FoldChangeConv[i] <- round((ConvDiv$nConvSites[i]/ ConvDiv$totalTermSizeAA[i])/(nrow(convergent)/sum(background$nPositions)),1)
	ConvDiv$FoldChangeDiv[i] <- round((ConvDiv$nDivSites[i]/ ConvDiv$totalTermSizeAA[i])/(nrow(divergent)/sum(background$nPositions)),1)
}

minGenesPerTerm = 10
maxGenesPerTerm = 500

ConvDivTest <- subset(ConvDiv, termID %in% unique(as.character(subset(genePerTerm, nGenes>=minGenesPerTerm & nGenes<=maxGenesPerTerm)$termID)))
annotatedGenes <- unique(as.character(subset(ontoFile, termID %in% unique(as.character(subset(genePerTerm, nGenes>=minGenesPerTerm & nGenes<=maxGenesPerTerm)$termID)))$geneName))
background <- subset(background, nPositions>0)
background <- subset(background, geneName %in% annotatedGenes)
convergent <- subset(convergent, geneName %in% annotatedGenes)
divergent <- subset(divergent, geneName %in% annotatedGenes)
for(i in 1:nrow(ConvDivTest)){
	ConvDivTest$FoldChangeConv[i] <- round((ConvDivTest$nConvSites[i]/ ConvDivTest$totalTermSizeAA[i])/(nrow(convergent)/sum(background$nPositions)),1)
	ConvDivTest$FoldChangeDiv[i] <- round((ConvDivTest$nDivSites[i]/ ConvDivTest$totalTermSizeAA[i])/(nrow(divergent)/sum(background$nPositions)),1)
}

# compute hypergeometric p-values for divergence and convergence:
n_totalGenes <- nrow(background)
n_totalConvergentGenes <- length(unique(as.character(convergent$geneName)))
n_totalDivergentGenes <- length(unique(as.character(divergent$geneName)))
ConvDivTest$GenesHyperConv <- phyper(ConvDivTest$nConvGenes-1, n_totalConvergentGenes, n_totalGenes-n_totalConvergentGenes, ConvDivTest$totalTermGenes, lower.tail=FALSE)
ConvDivTest$GenesHyperDiv <- phyper(ConvDivTest$nDivGenes-1, n_totalDivergentGenes, n_totalGenes-n_totalDivergentGenes, ConvDivTest$totalTermGenes, lower.tail=FALSE)

# now, compute enrichment (hypergeometric p-value & fold) seperately for convergent and divergent substitutions for each term
ConvDivTest$hyperConv <- rep(1, nrow(ConvDivTest))
ConvDivTest$hyperDiv <- rep(1, nrow(ConvDivTest))
for(i in 1:nrow(ConvDivTest)){
	ConvDivTest$hyperConv[i] <- phyper(ConvDivTest$nConvSites[i]-1, nrow(convergent), sum(background$nPositions)-nrow(convergent), ConvDivTest$totalTermSizeAA[i], lower.tail=FALSE)
	ConvDivTest$hyperDiv[i] <- phyper(ConvDivTest$nDivSites[i]-1, nrow(divergent), sum(background$nPositions)-nrow(divergent), ConvDivTest$totalTermSizeAA[i], lower.tail=FALSE)
}

# correct for multiple testing and perform tests:
ConvDivTest$GenesHyperConv_correctBH <- signif(p.adjust(ConvDivTest$GenesHyperConv, method="BH"),3)
ConvDivTest$GenesHyperDiv_correctBH <- signif(p.adjust(ConvDivTest$GenesHyperDiv, method="BH"),3)
ConvDivTest$hyperConv_correctBH <- signif(p.adjust(ConvDivTest$hyperConv, method="BH"),3)
ConvDivTest$hyperDiv_correctBH <- signif(p.adjust(ConvDivTest$hyperDiv, method="BH"),3)
ConvDivTest$genePerTermRange <- rep(paste0(minGenesPerTerm,'-',maxGenesPerTerm),nrow(ConvDivTest))
write.table(ConvDivTest, file=paste0(filename, '.termEnrichments'), col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")


