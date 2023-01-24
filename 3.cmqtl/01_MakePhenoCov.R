library(tidyr, lib="/apps/source/R/3.6.3/lib64/R/library/")
library(fansi, lib="/apps/source/R/3.6.3/lib64/R/library/")
library(dplyr, lib="/apps/source/R/3.6.3/lib64/R/library/")

###module load R/3.6.3
###inputs are:
###intcol
###isolate

args <- commandArgs(trailingOnly=TRUE)

id_map <- read.csv("../SecondRound_JatinFeatures/wgs.map.ModifiedSamira.csv")

features <- readRDS(paste("/data/srlab/jatin/cmqtl/gwas_data/final/donor_level/",args[1],".final.rds", sep=""))

Covs <- c("FID",
        "IID",
        "Metadata_Plate",
        "Metadata_Sex")

qCovs <- c("FID",
	"IID", 
	"PC1",
        "PC2",
        "PC3",
        "PC4",
	"Cells_Neighbors_NumberOfNeighbors_Adjacent",
	"Metadata_onEdge")

ExcludeColumns <- c("Metadata_line_ID",
	"Metadata_ID",
	"Metadata_Plate",
	"Metadata_Clinical_Diagnosis",
	"Metadata_Age",
	"Metadata_Sex",
	"Metadata_iPSC_Origin",
	"Metadata_Doubling_Hour",
	"Metadata_Disease_Status" ,
	"Metadata_Disease_Category",
	"Metadata_CellCount_avg",
	"PC1",
	"PC2",
	"PC3",
	"PC4",
	"PC5",
	"Cells_Neighbors_NumberOfNeighbors_Adjacent",
	"Metadata_onEdge")

### for isolate i used the following lists instead
#qCovs <- c("FID",
#        "IID",
#	"PC1",
#        "PC2",
#        "PC3",
#        "PC4",
#	"Metadata_onEdge")
#
#ExcludeColumns <- c("Metadata_line_ID",
#        "Metadata_ID",
#        "Metadata_Plate",
#        "Metadata_Clinical_Diagnosis",
#        "Metadata_Age",
#        "Metadata_Sex",
#        "Metadata_iPSC_Origin",
#        "Metadata_Doubling_Hour",
#        "Metadata_Disease_Status" ,
#        "Metadata_Disease_Category",
#	"Metadata_CellCount_avg",
#        "PC1",
#        "PC2",
#        "PC3",
#        "PC4",
#	"PC5",
#	"Metadata_onEdge")
#
## order data like fam
data <- merge(id_map, features, by=c("Metadata_line_ID","Metadata_ID"))

fam <- read.table("../SecondRound_JatinFeatures/PlinkFiles/V2F_chr1.nodup.fam", header=F)
fam <- fam[,1:2]
colnames(fam) <- c("FID", "IID")

new_data <- merge(fam, data, by.x="FID", by.y="wgs_id", all.x=T)
new_data <- new_data %>% group_by(Metadata_ID) %>% sample_n(1)
new_data <- new_data[match(fam$IID, new_data$IID),]
new_data <- subset(new_data, !is.na(Metadata_ID))

### make pheno and cov files
new_data <- as.data.frame(new_data)
pheno <-  new_data %>% dplyr::select(-c(all_of(ExcludeColumns)))
covar <- new_data %>% dplyr::select(c(all_of(Covs)))
qcovar <- new_data %>% dplyr::select(c(all_of(qCovs)))

write.table(colnames(pheno[,-c(1:2)]), paste("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/",args[1],"/",args[1],".phenotypeNames.txt", sep=""), quote=F, row.names=T, col.names=F)
write.table(colnames(covar[,-c(1:2)]), paste("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/",args[1],"/",args[1],".CatCovarNames.txt", sep=""), quote=F, row.names=T, col.names=F)
write.table(colnames(qcovar[,-c(1:2)]), paste("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/",args[1],"/",args[1],".QuanCovarNames.txt", sep=""), quote=F, row.names=T, col.names=F)

write.table(pheno, paste("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/",args[1],"/",args[1],".pheno", sep=""), quote=F, row.names=F, col.names=F)
write.table(qcovar, paste("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/",args[1],"/",args[1],".quantitative.cov", sep=""), quote=F, row.names=F, col.names=F)
write.table(covar, paste("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/FourthdeRound_TwoGroupEdgeCov/",args[1],"/",args[1],".categorical.cov", sep=""), quote=F, row.names=F, col.names=F)

