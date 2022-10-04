#!/usr/bin/env Rscript

args = commandArgs(TRUE)
library(dada2); packageVersion("dada2")
pathF = args[1]
pathR = args[2]
filtpathF <- file.path(pathF, "filtered_R1")
filtpathR <- file.path(pathR, "filtered_R2") 
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering
out <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
                     rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                     truncLen=c(278,262), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=FALSE, matchIDs = TRUE)


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
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, args[3])
write.csv2(seqtab, args[4])

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

save.image(args[5])

getN <- function(x) sum(getUniques(x)) 
track <- cbind(out, sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input_read_count", "filtered_and_trimmed_read_count", "merged_after_dada2_read_count", "non-chimeric_read_count")
rownames(track) <- sample.names
write.csv2(track, args[6])