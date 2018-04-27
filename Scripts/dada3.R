#Importing libraries
library(dada2)

oldw <- getOption("warn")
options(warn = -1)
args = commandArgs(trailingOnly=TRUE)

#Input is the path to the Filtered reads (R1 and R2 in the same folder)
path <- args[1]
Output <- args[2]
#These variables will hold the list with the full paths for each R of the paired-end reads
fnFs <- sort(list.files(path, pattern="_R1.*\\.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.*\\.fq", full.names = TRUE))
sample.names <- (basename(fnFs))

#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence. Dereplication substantially reduces computation time by eliminating redundant comparisons.

derepFs <- derepFastq(fnFs, verbose=TRUE)
derepRs <- derepFastq(fnRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#The DADA2 algorithm depends on a parametric error model (err) and every amplicon dataset has a different set of error rates. The  learnErrors method learns the error model from the data, by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution.
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errF <- dadaFs.lrn[[1]]$err_out
errR <- dadaRs.lrn[[1]]$err_out

#core sequence-variant inference algorithm to the dereplicated data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


#Spurious sequence variants are further reduced by merging overlapping reads. The core function here is mergePairs, which depends on the forward and reverse reads being in matching order at the time they were dereplicated.

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)

#The core dada method removes substitution and indel errors, but chimeras remain. Fortunately, the accuracy of the sequences after denoising makes identifying chimeras simpler than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as a bimera (two-parent chimera) from more abundant sequences.
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)

write.table(t(seqtab.nochim), file = Output, sep = "\t")









