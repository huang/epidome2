#!/usr/bin/env Rscript

args = commandArgs(TRUE)
library(dada2); packageVersion("dada2")
pathF = "cutadapted_16S"
pathR = "cutadapted_16S"
fastqFs <- sort(list.files(pathF, pattern="_R1.fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="_R2.fastq.gz"))
#output filepath
filtpathF <- file.path(pathF, "filtered_R1")
filtpathR <- file.path(pathR, "filtered_R2") 
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering  #should include  truncLen=c(278,262),  maxEE=2, minLen=138, 
out <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
                     rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                      maxEE=10, truncQ=2, maxN=0, rm.phix=TRUE, minLen=170,
                     compress=TRUE, verbose=TRUE, multithread=FALSE, matchIDs = TRUE)
#Alternative: cp *_R1.fastq.gz filtered_R1
#Alternative: cp *_R2.fastq.gz filtered_R2

filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_R"), `[`, 1) 
sample.namesR <- sapply(strsplit(basename(filtRs), "_R"), `[`, 1) 

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

set.seed(100)
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE, pool=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE, pool=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap=10, maxMismatch=1)
  mergers[[sam]] <- merger
}
#Processing: staggered-mock3-1 
#Sample 1 - 23535 reads in 14240 unique sequences.
#Sample 1 - 23535 reads in 14484 unique sequences.
rm(derepF); rm(derepR)
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "16S_seqtab_from_dada2.rds")
write.csv2(seqtab, "16S_seqtab_from_dada2.csv")
#14*263-->14*44
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
saveRDS(seqtab.nochim, "16S_seqtab_nochim.rds")
write.csv2(seqtab.nochim, "16S_seqtab_nochim.csv")

#16S_seqtab_from_dada2.rds
save.image("16S_seqtab_image.RData")

getN <- function(x) sum(getUniques(x)) 
track <- cbind(out, sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input_read_count", "filtered_and_trimmed_read_count", "merged_after_dada2_read_count", "non-chimeric_read_count")
rownames(track) <- sample.names
write.csv2(track, "16S_track.csv")
