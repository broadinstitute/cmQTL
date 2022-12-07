library(qqman)
library(tidyverse)
library(data.table)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

base.files <- Sys.glob(file.path(paste("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/", args[1],"/",args[1],"*.DownSampled.tmp.gz", sep="")))

results1 <- suppressWarnings(rbindlist(lapply(base.files, function(filename) {suppressMessages(readr::read_tsv(filename, na = c("NaN")))})))
results1=na.omit(results1)

results1.1 <- results1 %>% filter(P <= 1e-3)
results1.2 <- results1 %>% filter(P > 1e-3) %>% sample_frac(0.05)

results <- rbind(results1.1, results1.2)

##fdr <- fread(paste("gzip -dc FDRsigPvalue_", args[1], ".txt.gz", sep=""))
##x <- nrow(fdr)
##if(x == 0){
##sug <- -log10(5e-8/908)
##} else {
##sug <- -log10(max(fdr$V1))
##}

png(file=paste(args[1], "_AllPheno_downSampled.png", sep=""), width=900, height=400)
manhattan(results, chr="CHR", bp="POS", snp="SNP", p="P", suggestiveline = -log10(5e-8), genomewideline = -log10(5e-8/246) ,col=c("coral","blueviolet"), ylim=c(0, -log10(min(results$P))+1))
dev.off()
