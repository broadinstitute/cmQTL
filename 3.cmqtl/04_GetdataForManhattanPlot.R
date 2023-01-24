library(data.table)
library(readr)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

raw.files <- Sys.glob(file.path(paste("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/", args[2],"/",args[2],"_chr*_pheno",args[1],".LR.fastGWA.gz", sep="")))

results1 <- suppressWarnings(rbindlist(lapply(raw.files, function(filename) {suppressMessages(readr::read_tsv(filename, na = c("NaN")))})))
results1=na.omit(results1)
results1$P <- as.numeric(results1$P)

sig <- subset(results1, P < 5e-8)

#results1.1 <- results1 %>% filter(P <= 1e-2)
#results1.2 <- results1 %>% filter(P > 1e-2) %>% sample_frac(0.01)
#
#results <- rbind(results1.1, results1.2)

write.table(sig, paste("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/", args[2],"/",args[2],"_pheno",args[1],".lowP4manhattan.tmp", sep=""), quote=F, col.names=T, row.names=F, sep="\t")

#write.table(results, paste("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/", args[2],"/",args[2],"_pheno",args[1],".DownSampled.tmp", sep=""), quote=F, col.names=T, row.names=F, sep="\t")

