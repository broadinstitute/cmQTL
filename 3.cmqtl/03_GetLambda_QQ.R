library(qqman)
library(data.table)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

raw.files <- Sys.glob(file.path(paste("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/", args[2],"/",args[2],"_chr*_pheno",args[1],".LR.fastGWA.gz", sep="")))

results <- suppressWarnings(rbindlist(lapply(raw.files, function(filename) {suppressMessages(readr::read_tsv(filename, na = c("NaN")))})))

results=na.omit(results)
pvalues=results$P

### qchisq(0.5,1) is the median of a chi-squared distribution with one degree of freedom
z <- as.numeric(results$BETA)/as.numeric(results$SE) 
lambda2 = round(median(z^2) / qchisq(0.5,1), 4)

write.table(t(c(args[1], lambda2)), paste("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/QQ/", args[2],"/",args[2],"_pheno",args[1],".QQ.txt", sep=""), quote=F, col.names=F, row.names=F)

#png(file=paste("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/QQ/", args[2],"/",args[2],"_pheno",args[1],".QQ.png", sep=""), width=400, height=400)
#qqman::qq(pvalues, cex = 0.5, cex.axis = 1,main = paste(args[2], "pheno", args[1], "\n lambda=", lambda2, sep=" "))
#dev.off()
