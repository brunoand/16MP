#Importing libraries
library(dada2)
args = commandArgs(trailingOnly=TRUE)

#Input is the path to the Filtered reads (R1 and R2 in the same folder)
path <- args[1]
Output <- args[2]
#These variables will hold the list with the full paths for each R of the paired-end reads
fnFs <- sort(list.files(path, pattern="_R1.*\\.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.*\\.fastq", full.names = TRUE))
sample.names <- (basename(fnFs))
#write.table(sample.names, file = Output, sep = "\t")
#Filtering and Trimming

filt_path <- file.path((dirname(Output)), "filtered")
#if(!file_test("-d", file_path)) dir.create(file_path)
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
#write.table(filtFs, file = Output, sep = "\t")
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
#write.table(fnFs, file = Output, sep = "\t")
# Filter
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                    truncLen=c(240,160),
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE)
}
Outdir <- paste(Output, '/Output.tsv', sep="")
#write.table(filtFs, file = Output, sep = "\t")

#write.table(fnFs, file = Output, sep = "\t")
#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence. Dereplication substantially reduces computation time by eliminating redundant comparisons.
#filtFs
#derepFs <- derepFastq(filtFs, verbose=TRUE)
#derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
#names(derepFs) <- sample.names
#names(derepRs) <- sample.names
write.table(filtRs, file = Outdir, sep = "\t")
#The DADA2 algorithm depends on a parametric error model (err) and every amplicon dataset has a different set of error rates. The  learnErrors method learns the error model from the data, by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution.
#dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
#dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
#errF <- dadaFs.lrn[[1]]$err_out
#errR <- dadaRs.lrn[[1]]$err_out

#core sequence-variant inference algorithm to the dereplicated data
#dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
#dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


#Spurious sequence variants are further reduced by merging overlapping reads. The core function here is mergePairs, which depends on the forward and reverse reads being in matching order at the time they were dereplicated.

#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#seqtab <- makeSequenceTable(mergers)

#The core dada method removes substitution and indel errors, but chimeras remain. Fortunately, the accuracy of the sequences after denoising makes identifying chimeras simpler than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as a bimera (two-parent chimera) from more abundant sequences.
#seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)

#write.table(t(seqtab.nochim), file = Output, sep = "\t")









