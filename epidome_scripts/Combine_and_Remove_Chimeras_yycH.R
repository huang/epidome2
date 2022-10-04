#!/usr/bin/env Rscript

args = commandArgs(TRUE)
library(dada2); packageVersion("dada2")
library(dplyr)
library(tibble)

sts <- list.files(path = args[1], pattern = "*yycH_seqtab_from_dada2.rds", full.names = TRUE)
st.all <- mergeSequenceTables(tables = sts)
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

saveRDS(seqtab.nochim, args[2])
write.csv2(seqtab.nochim, args[3])