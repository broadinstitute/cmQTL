library(readr)
library(stringr)
library(dplyr)
library(data.table)

pheno <- read.table("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/isolate/isolate.phenotypeNames.txt", header=F)
colnames(pheno) <- c("Pheno", "Feature")
pheno$Pheno <- paste0("pheno", pheno$Pheno)

raw.files <- Sys.glob(file.path(
  "/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/*/*_pheno*.lowP4manhattan.tmp"
))

data <- suppressWarnings(rbindlist(lapply(raw.files, function(filename) {
    suppressMessages(readr::read_tsv(filename, col_names=FALSE))}), idcol = "File"))
data[, File := factor(File, labels = basename(raw.files))]
data <- data %>% filter(X1 != "CHR")
data$File <- str_remove(data$File, ".lowP4manhattan.tmp")

colnames(data) <-  c("File", "CHR", "SNP", "POS", "A1", "A2", "N", "AF1", "BETA","SE","P")
data$Cell <- str_split_fixed(data$File, "_", 2)[,1]
data$Pheno <- str_split_fixed(data$File, "_", 2)[,2]

data <- merge(data, pheno, by="Pheno")
data$P <- as.numeric(data$P)
data <- data %>% filter(P <= 4.1e-8)
write.table(data, "SuggestiveGWAS_results.txt", quote=F, col.names=T, row.names=F, sep="\t")
