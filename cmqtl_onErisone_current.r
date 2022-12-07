#------------------------------------------------------------------------------------------------------------------
# jarora1@bwh.harvard.edu
#------------------------------------------------------------------------------------------------------------------
# NOTES
# 1. rds.objects/../processed data files have cells with > 5% 0s across features removed
# 2. rds.objects/../qc_ed data files have cells with any 0 across features removed
#------------------------------------------------------------------------------------------------------------------
rm(list = ls())
graphics.off()
options(scipen=999)
options(warn=1)
cores=30
#--------------- Packages and donor metadata ----------------
if(T){
args = commandArgs(trailingOnly = T)
print(args[1])
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({
  library("stringr")
  library("RColorBrewer")
  library("umap")
  library("ggplot2")
  library("ggfortify")
  library("lme4")
  library("dplyr")
  library("dbplyr")
  library("caret")
  library("EnvStats")
  library("CDFt")
  library("reshape")
})

# colors
my.cols = c(brewer.pal(8,"Dark2")[c(6)], brewer.pal(8,"Paired")[c(4,6)], rev(brewer.pal(11,"Spectral")))

# inter-replicate correlation for feature
replicate_correlation_values = readr::read_csv("bekar.cheeze/replicate_correlation_values.csv")

# blacllisted features
blacklist = c("Nuclei_Correlation_Manders_AGP_DNA","Nuclei_Correlation_Manders_AGP_ER","Nuclei_Correlation_Manders_AGP_Mito","Nuclei_Correlation_Manders_AGP_RNA","Nuclei_Correlation_Manders_DNA_AGP","Nuclei_Correlation_Manders_DNA_ER","Nuclei_Correlation_Manders_DNA_Mito","Nuclei_Correlation_Manders_DNA_RNA","Nuclei_Correlation_Manders_ER_AGP","Nuclei_Correlation_Manders_ER_DNA","Nuclei_Correlation_Manders_ER_Mito","Nuclei_Correlation_Manders_ER_RNA","Nuclei_Correlation_Manders_Mito_AGP","Nuclei_Correlation_Manders_Mito_DNA","Nuclei_Correlation_Manders_Mito_ER","Nuclei_Correlation_Manders_Mito_RNA","Nuclei_Correlation_Manders_RNA_AGP","Nuclei_Correlation_Manders_RNA_DNA","Nuclei_Correlation_Manders_RNA_ER","Nuclei_Correlation_Manders_RNA_Mito","Nuclei_Correlation_RWC_AGP_DNA","Nuclei_Correlation_RWC_AGP_ER","Nuclei_Correlation_RWC_AGP_Mito","Nuclei_Correlation_RWC_AGP_RNA","Nuclei_Correlation_RWC_DNA_AGP","Nuclei_Correlation_RWC_DNA_ER","Nuclei_Correlation_RWC_DNA_Mito","Nuclei_Correlation_RWC_DNA_RNA","Nuclei_Correlation_RWC_ER_AGP","Nuclei_Correlation_RWC_ER_DNA","Nuclei_Correlation_RWC_ER_Mito","Nuclei_Correlation_RWC_ER_RNA","Nuclei_Correlation_RWC_Mito_AGP","Nuclei_Correlation_RWC_Mito_DNA","Nuclei_Correlation_RWC_Mito_ER","Nuclei_Correlation_RWC_Mito_RNA","Nuclei_Correlation_RWC_RNA_AGP","Nuclei_Correlation_RWC_RNA_DNA","Nuclei_Correlation_RWC_RNA_ER","Nuclei_Correlation_RWC_RNA_Mito","Nuclei_Granularity_14_AGP","Nuclei_Granularity_14_DNA","Nuclei_Granularity_14_ER","Nuclei_Granularity_14_Mito","Nuclei_Granularity_14_RNA","Nuclei_Granularity_15_AGP","Nuclei_Granularity_15_DNA","Nuclei_Granularity_15_ER","Nuclei_Granularity_15_Mito","Nuclei_Granularity_15_RNA","Nuclei_Granularity_16_AGP","Nuclei_Granularity_16_DNA","Nuclei_Granularity_16_ER","Nuclei_Granularity_16_Mito","Nuclei_Granularity_16_RNA")

# get empty wells from plate map
plate.map = NULL
for (file in list.files(path = "platemap", pattern = "cmQTL.*.txt", full.names = T)){
  dat = read.delim(file)
  plate.map = rbind(plate.map, dat)
}
empty.wells = plate.map[plate.map$line_ID == 0, c("plate_map_name", "well_position")]
filled.wells = plate.map[plate.map$line_ID != 0, c("plate_map_name", "well_position")]

# load plate-batch data
plate.batch = read.delim("platemap/plate.batch.map", header = F, col.names = c("plate", "plate_map_name"))
filled.wells.2 = merge(x = filled.wells, y = plate.batch, by.x = c("plate_map_name"), by.y = c("plate_map_name"))

# 1. ----- meta data about all donors -----
# cell line ids of donors
donor.ids = read.delim("metadata/metaData_cellLine_DonorID.tab", col.names = c("Metadata_Project.Alias","Metadata_ID"))

# metadata about donors
donor.meta = read.delim("metadata/metaData_all.tab") %>% dplyr::select(one_of("ID","Clinical.Diagnosis","Age","Sex","iPSC.Origin","Doubling.Time..h.")) %>% 
  dplyr::rename("Doubling_Hour"="Doubling.Time..h.") %>% dplyr::rename_with(~ paste0("Metadata_", gsub("\\.+", "_", .x))) %>%
  inner_join(donor.ids, by = "Metadata_ID") %>% 
  mutate("Metadata_Disease_Status" = ifelse(grepl("Control", Metadata_Clinical_Diagnosis, ignore.case=T), "Control", "Disease"))

# disease category
donor.meta = donor.meta %>% mutate("Metadata_Disease_Category" = ifelse( grepl("Phelan|Depression|Alzheimer|Intellectual|Epilepsy|Autism|Cerebral",Metadata_Clinical_Diagnosis), yes = "Neuro", 
                                                                          no = ifelse( grepl("Control",Metadata_Clinical_Diagnosis), yes = "Control", 
                                                                                       no = ifelse( grepl("Eye",Metadata_Clinical_Diagnosis), yes = "Opth", 
                                                                                                    no = ifelse( grepl("Heart",Metadata_Clinical_Diagnosis), yes = "Cardio", 
                                                                                                                 no = ifelse( grepl("Diab",Metadata_Clinical_Diagnosis), yes = "Endo", 
                                                                                                                              no = ifelse( grepl("Liver",Metadata_Clinical_Diagnosis), yes = "Liver", 
                                                                                                                                           no = Metadata_Clinical_Diagnosis) ) ) ) ) ) )

donor.meta = donor.meta %>% mutate("Metadata_iPSC_Origin" = ifelse(grepl("No",Metadata_iPSC_Origin), "B-Lymphocyte", Metadata_iPSC_Origin))
paste("number of donors having metadata",length(unique(donor.meta$Metadata_ID)))

# 2. ----- meta data about donors with wgs -----
# metadata about wgs -- short description
wgs.map = read.delim("bekar.cheeze/participant_july2020.tsv") %>% mutate("meta_id_1" = sub("_P[1-9]{2}.*", "", collaborator_participant_id)) %>% 
  mutate("meta_id_final" = sub("([A-Z])\\1[0-9]","", meta_id_1)) %>% dplyr::select(one_of("entity.participant_id","collaborator_participant_id","meta_id_final"))
paste("number of donors reported in terra",length(unique(wgs.map$meta_id_final)))

###
wgs.map.big.meta = read.delim("metadata/sample.tsv") %>% dplyr::select(one_of("collaborator_sample_id","collaborator_participant_id")) %>% 
  mutate("to_sample_id" = sub("_P[0-9]{2}_.*", "", sub("([A-Z])\\1[0-9]","", collaborator_participant_id))) %>%
  inner_join(wgs.map %>% mutate("to_sample_id" = ifelse(!grepl("SCBB",entity.participant_id), meta_id_final, entity.participant_id)), by = "to_sample_id") %>%
  dplyr::select(one_of("collaborator_sample_id","to_sample_id","meta_id_final"))

donor.meta = inner_join(donor.meta, wgs.map.big.meta, by = c("Metadata_ID"="meta_id_final"))
paste("number of donors in terra having metadata",length(unique(donor.meta$Metadata_ID)))

}
#--------------------- Filter donors ------------------------
# lines having issues
# no colony and overall very few cells, CW40187
# line 53 CW60334 removed, low number of cells
donor.meta = donor.meta %>% filter(!(Metadata_ID %in% c("CW20009","CW50058","CW60029","CW70272","CW30178","CW40187")))
# remove donors with missing age
# donor.meta = donor.meta %>% filter(!grepl("\\?",Metadata_Age)) %>% mutate(Metadata_Age = as.numeric(Metadata_Age)) #%>% filter(!grepl("Not", Metadata_Doubling_Hour, ignore.case=T))
print(paste("donors taken after filtering", length(unique(donor.meta$Metadata_ID))))
#------------------- Specific disease -----------------------
if(F){
print(table(donor.meta$Metadata_Disease_Category))
donor.meta = donor.meta %>% filter(Metadata_Disease_Category %in% c("Control","Neuro"))
print("new affected status done")
print(table(donor.meta$Metadata_Disease_Category))

# table(donor.meta.wgs$Metadata_donor_disease.name)
# donor.meta %>% filter(Metadata_donor_disease.name == "POAG") %>% pull("Metadata_donor_Description") %>% unique()
}
#-------------- Genotype PCs and Cell count -----------------
if(T){
# PCs from genotypes
# pcs = read.table("bekar.cheeze/SCBB_CIRM_WGScallset_Aug2020.bi.commonCM.eigenvec", col.names = c("raw_id",paste0("PC",1:10))) %>%
#   mutate("wgs_id" = toupper(gsub("\\_|\\.","-", raw_id)))

# pcs calculated by samira
pcs = read.table("bekar.cheeze/V2F_AllChr.nodup.qc.pruned.noMHC.pc.eigenvec", col.names = c("raw_id","raw_id_2",paste0("PC",1:20))) %>%
  dplyr::select(-one_of("raw_id_2",paste0("PC",6:20))) %>% mutate("wgs_id" = toupper(gsub("\\_|\\.","-", raw_id)))

# combine pcs with cell line data
pcs.terra = donor.meta %>% mutate("wgs_id" = toupper(gsub("\\_|\\.","-",collaborator_sample_id))) %>%
  dplyr::select(c(wgs_id,Metadata_Project.Alias)) %>% inner_join(pcs, by = "wgs_id") %>% 
  dplyr::select(one_of("Metadata_Project.Alias",paste0("PC",1:5)))
pcs.terra = pcs.terra[!duplicated(pcs.terra$Metadata_Project.Alias), ]
print(paste("genotype pcs loaded for donors", length(unique(pcs.terra$Metadata_Project.Alias))))

# get average cell count per donor or well
cell.count = list.files(path = "profiles", pattern = "*_count.csv", full.names = T) %>% purrr::map_df(read.csv)

cell.count.summ.line = read.csv("bekar.cheeze/plate.metadata.csv") %>% inner_join(cell.count, by=c("Metadata_Plate","Metadata_Well")) %>%
  group_by(Metadata_line_ID) %>% mutate("Metadata_CellCount_avg" = mean(Count_Cells)) %>%
  dplyr::select(Metadata_line_ID,Metadata_CellCount_avg) %>% unique()

cell.count.summ.well = read.csv("bekar.cheeze/plate.metadata.csv") %>% inner_join(cell.count, by=c("Metadata_Plate","Metadata_Well")) %>%
  dplyr::rename("Metadata_CellCount_avg" = "Count_Cells") %>% dplyr::select(Metadata_line_ID,Metadata_Plate,Metadata_Well,Metadata_CellCount_avg) %>% unique()
}
#---------------------- Functions ---------------------------
if(T){
# plot association
plot_association <- function(feature, features.values, gene.values, frame.to.plot, interaction.feat, model, round, level, plot=T){
  
  # p values of associations
  pval.frame = frame.to.plot %>% filter(feat == feature)
  data = features.values[!duplicated(features.values$Metadata_line_ID), ] %>% as.data.frame() %>%
    dplyr::select(c(Metadata_line_ID, !!as.name(interaction.feat), !!as.name(feature))) %>%
    inner_join(gene.values %>% dplyr::select( -one_of(interaction.feat) ), by = c("Metadata_line_ID"="Metadata_Project.Alias")) %>% 
    dplyr::select(c(!!as.name(feature),unique(pval.frame$variable),Metadata_line_ID,!!as.name(interaction.feat))) %>%
    reshape2::melt(id = c(feature,"Metadata_line_ID",interaction.feat)) %>% 
    mutate(value=plyr::mapvalues(value, c(0,1), c("no","yes")))
  
  # prepare pdf file
  if(plot == T){
    if(nrow(pval.frame) == 1) wid = 2.5 else wid = 2*nrow(pval.frame)
    if(model == "interaction") pdf(paste("plots/p", round, level, feature, model, interaction.feat, "pdf", sep = "."), width=wid, height=3.5)
    if(model == "baseline") pdf(paste("plots/p", round, level, feature, model, "pdf", sep = "."), width=wid, height=3.5)
  }
  
  if(round == "intcol") title="Colony cells"
  if(round == "isolate") title="Isolated cells"
  
  # make plot, !!as.name(interaction.feat)
  p = ggplot(data, aes(x = factor(value), y = !!as.name(feature))) +
    geom_jitter(width=0.15, stroke=0, alpha=0.45, size=2, color=feat.cat.cols[unique(sub(".*_(.*?)_.*", "\\1", feature))]) + 
    geom_boxplot(outlier.size=0, alpha=0.75, width=0.35, color=feat.cat.cols[unique(sub(".*_(.*?)_.*", "\\1", feature))]) +
    ggthemes::scale_color_tableau(name = sub("Metadata_","",interaction.feat)) +
    facet_wrap(~variable, scales = "free", nrow = 1) + theme_bw() +
    theme(panel.grid.major.x = element_blank(), strip.background = element_blank()) +
    theme(axis.title = element_text(size=12), axis.text = element_text(size=12), strip.text = element_text(size=12), 
          title = element_text(size=10)) + 
    xlab("Rare variant burden") + ylab(feature) + ggtitle(title) + guides(fill=F, color=F) +
    geom_text(data = pval.frame, mapping = aes(x = Inf, y = Inf, label = paste("beta ==", round(est,2))), hjust = 1.1, vjust = 1.5, size = 4, parse = T) +
    geom_text(data = pval.frame, mapping = aes(x = Inf, y = Inf, label = paste("p =", p)), hjust = 1.1, vjust = 3.5, size = 4, parse = F)
  
  if(plot == T) print(p) else return(p)
  dev.off()
  print("plot done")
  
}
  
# variation in features accounted by fixed and random effects
variance_analysis <- function(data.allPlates, features.to.test, round, level){
  
  # if single-cell data, then add a term to account for multiple measurements from the same well
  if(level == "single") data.allPlates = data.allPlates %>% tidyr::unite("Metadata_PlateWell", Metadata_Plate:Metadata_Well, remove = F)
  
  # convert line ID to factors
  data.allPlates = data.allPlates %>% mutate(Metadata_line_ID = factor(Metadata_line_ID))
    
  # variance in feature accounted by fixed and random effects
  var.random = NULL
  for(feature in features.to.test){
    
    # prepare mixed effects Metadata_image_batch, Metadata_CellCount_avg
    print(feature)
    covars = c("Metadata_onEdge", "Metadata_Plate", "Metadata_iPSC_Origin", "Metadata_Sex", "Metadata_Age", "Metadata_Disease_Status", "Metadata_line_ID", paste0("PC",1:4))
    if(round != "isolate") covars = union("Cells_Neighbors_NumberOfNeighbors_Adjacent", covars)
    if(level == "well") covars = union(covars, "Metadata_Well")
    if(level == "single") covars = union(covars, c("Metadata_Well", "Metadata_PlateWell"))
    if(level %in% c("well","single")) effects = paste(c(sub("(Meta.*)","(1|\\1)", grep("Age|Neighbors|CellCount", covars, value=T, invert=T)), grep("Age|Neighbors|Origin|CellCount", covars, value=T)), collapse = " + ")
  
    # calculate mixed effect models
    model.effc = lmerTest::lmer(formula = as.formula(paste(feature,"~",effects)), data = data.allPlates, 
                                control = lmerControl(optimizer = "bobyqa"), REML = F)
    
    # pvalue for random effects
    rand.p.values = lmerTest::ranova(model.effc) %>% as.data.frame() %>% tibble::rownames_to_column("var") %>% 
      dplyr::select(one_of("var","Pr(>Chisq)")) %>% filter(!grepl("none",var)) %>% 
      mutate("grp" = sub("\\(.*\\| (.*)\\)","\\1",var), "p" = round(`Pr(>Chisq)`,10)) %>% dplyr::select(one_of("grp","p"))
      
    # variation by fixed effects
    # MuMIn::r.squaredGLMM(model.effc)
    betas = r2glmm::r2beta(model.effc, method="nsj", partial=T)
    r2.neigh = betas[betas$Effect == "Cells_Neighbors_NumberOfNeighbors_Adjacent", "Rsq"]
    r2.age = betas[betas$Effect == "Metadata_Age", "Rsq"]
    # r2.cct = betas[betas$Effect == "Metadata_CellCount_avg", "Rsq"]
  
    # variation accounted by random effects
    vari = VarCorr(model.effc) %>% data.frame() %>% mutate("pvcov" = vcov) %>% dplyr::select(grp,pvcov) %>% 
      full_join(rand.p.values, by="grp") %>% mutate(p = ifelse(grepl("Residual", grp), 0, p)) %>%
      tibble::add_row("grp"="Metadata_Age", "pvcov"=r2.age, "p"=coefficients(summary(model.effc))["Metadata_Age","Pr(>|t|)"])
    
    # variation by number of cell neighbors
    if(round != "isolate") vari = vari %>% tibble::add_row("grp"="Cell_neighbors", "pvcov"=r2.neigh, "p"=round(coefficients(summary(model.effc))["Cells_Neighbors_NumberOfNeighbors_Adjacent","Pr(>|t|)"],10))
    
    # store results
    var.random = rbind(var.random, vari %>% mutate("feat"=feature, pvcov.2 = round(pvcov/sum(pvcov),3)))
  }
  
  # remove NAs
  var.random = var.random[complete.cases(var.random), ]
  print(paste("mixed model done for features", length(unique(var.random$feat))))
  return(var.random)
}

# p value adjustment - modification of real function
p.adjust <- function (p, method = p.adjust.methods, n = length(p)){
  method <- match.arg(method)
  if (method == "fdr") 
    method <- "BH"
  nm <- names(p)
  p <- as.numeric(p)
  p0 <- setNames(p, nm)
  if (all(nna <- !is.na(p))) 
    nna <- TRUE
  p <- p[nna]
  lp <- length(p)
  #stopifnot(n >= lp)
  if (n <= 1) 
    return(p0)
  if (n == 2 && method == "hommel") 
    method <- "hochberg"
  p0[nna] <- switch(method, bonferroni = pmin(1, n * p), holm = {
    i <- seq_len(lp)
    o <- order(p)
    ro <- order(o)
    pmin(1, cummax((n + 1L - i) * p[o]))[ro]
  }, hommel = {
    if (n > lp) p <- c(p, rep.int(1, n - lp))
    i <- seq_len(n)
    o <- order(p)
    p <- p[o]
    ro <- order(o)
    q <- pa <- rep.int(min(n * p/i), n)
    for (j in (n - 1L):2L) {
      ij <- seq_len(n - j + 1L)
      i2 <- (n - j + 2L):n
      q1 <- min(j * p[i2]/(2L:j))
      q[ij] <- pmin(j * p[ij], q1)
      q[i2] <- q[n - j + 1L]
      pa <- pmax(pa, q)
    }
    pmax(pa, p)[if (lp < n) ro[1L:lp] else ro]
  }, hochberg = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    pmin(1, cummin((n + 1L - i) * p[o]))[ro]
  }, BH = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    pmin(1, cummin(n/i * p[o]))[ro]
  }, BY = {
    i <- lp:1L
    o <- order(p, decreasing = TRUE)
    ro <- order(o)
    q <- sum(1/(1L:n))
    pmin(1, cummin(q * n/i * p[o]))[ro]
  }, none = p)
  p0
}

# function to load single cell data, it calls pre-process function inside
load_sql_data <- function(file, good.wells){
  
  # connect to sqlite database
  con = DBI::dbConnect(RSQLite::SQLite(), dbname = file)
  print(src_dbi(con))
  print("connection with db opened")
  
  # take cell line info from aggregate level data
  aug.data = read.delim(list.files(path = "profiles", pattern = paste0(sub(".*[/](.*)[.]sqlite", "\\1", file), "_augmented.csv"), full.names = T), header = T, sep = ",") %>%
    dplyr::select(one_of("Metadata_Plate","Metadata_Well","Metadata_line_ID"))
  print("aggregate plate data loaded for cell line info")
  
  # collect data for different compartments
  ## head(., 10000) %>%
  im.res = tbl(con, "Image") %>% dplyr::select(c("TableNumber","ImageNumber","Metadata_Plate","Metadata_Well")) %>% collect()
  ce.res = tbl(con, "Cells") %>% collect()
  nc.res = tbl(con, "Nuclei") %>% collect()
  cy.res = tbl(con, "Cytoplasm") %>% collect()
  print("single cell data retrieved from db")
  print(dim(im.res))
  print(dim(nc.res))
  # 
  # # merge tables for 3 cell compartments and colony level data (for metaline info)
  # data.single = ce.res %>%
  #   inner_join(., nc.res, by = c("TableNumber", "ImageNumber", "ObjectNumber")) %>%
  #   inner_join(., cy.res, by = c("TableNumber", "ImageNumber", "ObjectNumber")) %>%
  #   inner_join(., im.res, by = c("TableNumber", "ImageNumber")) %>%
  #   inner_join(., aug.data, by = c("Metadata_Plate","Metadata_Well"))
  # print("data from cells, nuclei, cytoplasm, images and cell line id combined")
  # print(dim(data.single))
  # 
  # rm(list = c("im.res","ce.res","nc.res","cy.res"))
  # 
  # # remove number columns as no need for them downstream
  # data.single = data.single[, stringr::str_subset(string = colnames(data.single), pattern = "^TableNumber|^ObjectNumber|^ImageNumber|^plate_map_name", negate = T)]
  # print("table, object, image number cols removed")
  # print(dim(data.single))
  
  # # ----- temp -----
  # # saveRDS(data.single, "rds.objects/p7.all.rds")
  # # print("intermediate data saved to rds")
  data.single = readRDS("/data/srlab/jatin/cmqtl/qced.data/misc/p7.all.rds")
  print("intermediate data loaded form rds")
  print(dim(data.single))
  # # ----- temp -----
  
  # preprocess data
  processed.data = pre_process_data(my.data = data.single, data.level = "singlecell", good.wells = good.wells)
  print("pre-processing done")
  print(dim(processed.data))
  
  # save data
  saveRDS(object = processed.data, file = paste0("rds.objects/singleCell/", sub(".*[/](.*)[.]sqlite", "\\1", file), ".processed.rds"), compress = "xz")
  print(paste("rds object for processed data saved for", file))

  # disconnect db connection and clear env
  DBI::dbDisconnect(conn = con)
  rm(list = c("data.single", "im.res", "ce.res", "nc.res", "cy.res"))
  print("connection with db closed")
  
}

# robust z-score data per plate
normalize_data <- function(data.single){
  # call cytominer to z-score data per plate
  data.norm = cytominer::normalize(
    population = data.single,
    variables = stringr::str_subset(string = colnames(data.single), pattern = "Metadata", negate = T),
    strata = c("Metadata_Plate"),
    sample = data.single,
    operation = "robustize"
  )
  print("normalisation done")
  print(dim(data.norm))
  
  # remove cols containing NA
  data.norm = data.norm[ , colSums(is.na(data.norm)) == 0]
  print("columns with na removed")
  print(dim(data.norm))
  
  return(data.norm)
}

# pre process data, remove bad samples and features
pre_process_data <- function(my.data, data.level, good.wells){
  print("filtering...")
  # remove empty wells
  if(data.level == "bulk"){
    my.data = merge(x = my.data, y = good.wells, by.x = c("Metadata_Plate_Map_Name","Metadata_well_position"), by.y = c("plate_map_name","well_position"), all = F)
    print("empty wells removed from bulk level")
  } else if (data.level == "singlecell"){
    my.data = merge(x = my.data, y = good.wells, by.x = c("Metadata_Plate","Metadata_Well"), by.y = c("plate","well_position"), all = F)
    print("empty wells removed from singlecell level")
  }
  print(dim(my.data))

  # remove blacklisted, costes, correlation features
  features.remove = unique(c("plate_map_name", blacklist, stringr::str_subset(colnames(my.data), "_Costes_"), stringr::str_subset(colnames(my.data), "_Correlation_")))
  cols.non.num = colnames(my.data %>% dplyr::select_if(purrr::negate(is.numeric)))
  cols.non.num = grep(pattern = "Metadata", x = cols.non.num, invert = T, value = T)
  print(paste("number of non-numeric features", length(cols.non.num)))
  my.data = my.data %>% dplyr::select(-one_of(union(cols.non.num, features.remove)))
  print("blacklisted and non-numeric features filtered out")
  print(dim(my.data))

  # features with all 0s
  feat.dat = my.data[, grep("Metadata", colnames(my.data), value = T, invert = T)]
  col.total = apply(feat.dat, 2, function(x){sum(abs(x))})
  my.data = my.data %>% dplyr::select(-one_of(names(which(col.total == 0))))
  print("columns with all 0s removed")
  print(dim(my.data))
  rm(feat.dat)
  
  # ----- temp for plate 7 diagnostic -----
  saveRDS(my.data[sample(x = 1:nrow(my.data), size = 25000), ], "rds.objects/p7.1.sub.rds")
  print("intermediate sampled data saved to rds")
  print("quitting...")
  quit()
  # ----- temp for plate 7 diagnostic -----
  
  # remove cells which miss certain % of features (the number of columns)
  feat.dat = my.data[, grep("Metadata", colnames(my.data), value = T, invert = T)]
  # --- for plate 7 ---
  row.nas = apply(feat.dat, 1, function(x){length(which(is.na(x)))})
  # --- for plate 7 ---
  row.zeros = apply(feat.dat, 1, function(x){length(which(x == 0))})
  my.data = my.data %>% dplyr::filter(!row_number() %in% union(which(row.nas > 0), which(row.zeros > ncol(feat.dat)*0.05)))
  print("rows with missing data removed")
  print(dim(my.data))
  rm(feat.dat)
  
  # remove the columns containing missing data
  my.data = my.data[ ,colSums(is.na(my.data)) == 0]
  print(paste("columns with missing data removed"))
  print(dim(my.data))
  
  return(my.data)
  
}

# function to regression feature on genes
regress_feature_allele <- function(feature, gene, covars, data, model, level, inter.var){
  # well level
  if(level != "bulk"){
    # factorize cell line id
    data = data %>% mutate(Metadata_line_ID = factor(Metadata_line_ID))
  
    # prepare models
    pattern = "Age|Cells_Neighbors|Sex|cellCat"
    effects = paste(c(sub("(Meta.*)","(1|\\1)", c(grep(pattern, covars, value=T, invert=T), "Metadata_line_ID")), grep(pattern, covars, value = T)), collapse = " + ")
    baseline = as.formula(paste(feature,"~", gene,"+", effects))

    # execute baseline model, optimizer = "nloptwrap"
    lm.b = lmerTest::lmer(formula = baseline, data = data, REML=F, control = lmerControl(optimizer = "bobyqa", calc.derivs = F))
    
    if(model == "baseline"){
      summ = coefficients(summary(lm.b))
      if(length(grep(gene, rownames(summ))) != 0) return(data.frame("var" = gene, "est" = summ[gene,"Estimate"], "std.err" = summ[gene,"Std. Error"], "p" = summ[gene,"Pr(>|t|)"], row.names = NULL))
    
    } else if(model == "interaction"){
      # alternate model
      lm.a = lmer(formula = alternative, data = data, control = lmerControl(optimizer = "bobyqa"), REML=F)
      model.imp = anova(lm.b, lm.a, test="LRT")$`Pr(>Chisq)`[2]
      
      # get summary of models
      summ = coefficients(summary(lm.a))
      index = grep(paste0(gene,":"),rownames(summ), value=T)
      
      # return results
      if(length(index) > 0 & !is.na(model.imp) & model.imp < 0.05){
        c = data.frame("var" = paste0(gene,"_inter"), "est" = summ[index,"Estimate"], "std.err" = summ[index,"Std. Error"], "p" = summ[paste0(gene,":",inter.var),"Pr(>|t|)"], row.names = NULL)
      } else{
        c = data.frame("var" = paste0(gene,"_inter"), "est" = 0, "std.err" = 1, "p" = 1, row.names = NULL)
      }
      
      return( data.frame(c,model.imp) )
    }
    
  }
  # donor level
  if(level == "bulk"){
    # prepare formulas
    baseline = as.formula(paste(feature,"~", gene,"+", paste(covars, collapse = "+")))
    alternative = as.formula(paste(feature,"~", gene,"*", inter.var, "+", paste(grep(inter.var, covars, value=T, invert=T), collapse = "+")))
    
    # baseline linear regression
    lm.b = lm(formula = baseline, data = data)
    
    # return results
    if(model == "interaction"){
      # alternate model
      lm.a = lm(formula = alternative, data = data)
      model.imp = anova(lm.b, lm.a, test="LRT")$`Pr(>Chi)`[2]
      # get summary of regression
      summ = coefficients(summary(lm.a))
      index = grep(paste0(gene,":"),rownames(summ), value = T)
      
      # return results
      if(length(index) > 0 & !is.na(model.imp) & model.imp < 0.05){
        c = data.frame("var" = paste0(gene,"_inter"), "est" = summ[index,"Estimate"], "std.err" = summ[index,"Std. Error"], "p" = summ[index,"Pr(>|t|)"], row.names = NULL)
      } else{
        c = data.frame("var" = paste0(gene,"_inter"), "est" = 0, "std.err" = 1, "p" = 1, row.names = NULL)
      }
      return( data.frame(c,model.imp) )
      
    } else if(model == "baseline"){
      summ = coefficients(summary(lm.b))
      if(length(grep(gene, rownames(summ))) != 0) return(data.frame("var" = gene, "est" = summ[gene, "Estimate"], "std.err" = summ[gene, "Std. Error"], "p" = summ[gene, "Pr(>|t|)"], row.names = NULL))
    }
    
  }
}

# ---- temp, test for residues ----
regress_with_residues <- function(feature, genes.to.test, covars, data, inter.var){
  
  effects = paste(c(sub("(Meta.*)","(1|\\1)", grep(paste0("Age|line_ID|",inter.var), covars, value = T, invert = T)), grep("Age", covars, value = T)), collapse = " + ")
  baseline.nogene = as.formula(paste(feature,"~", inter.var, "+", effects))
  lm.b.1 = lmer(formula = baseline.nogene, data = data, control = lmerControl(optimizer = "bobyqa"), REML=F)
  
  res.allgenes = NULL
  for(gene in genes.to.test$gene.2){
    data.2 = data %>% dplyr::select(one_of(gene,"Metadata_line_ID")) %>% mutate("residuals" = residuals(lm.b.1))
    lm.b.2 = lmer(formula = as.formula(paste("residuals","~", gene, "+ (1|Metadata_line_ID)")), data = data.2, control = lmerControl(optimizer = "bobyqa"), REML=F)
    chiq = data.frame(car::Anova(lm.b.2))
    summ = coefficients(summary(lm.b.2))
    if(length(grep(gene, rownames(summ))) != 0) res.allgenes = rbind(res.allgenes, data.frame("var" = gene, "est" = summ[gene,"Estimate"], "std.err" = summ[gene,"Std. Error"], "p" = chiq[gene,"Pr..Chisq."], row.names = NULL))
  }
  
  return(res.allgenes)
}
# ---- temp, test for residues ----
}
#------------------------------------------------------------
#------------------------------------------------------------
seed = as.numeric(args[1])
level = "bulk" # well, bulk, single
round = "isolate" # intcol, isolate, anycells
sum.or.present = "present" # sum, present
variants = "all" # all, high, moderate
variant.source = "rareCM" # rareCM, rareCM.Gnomad, 5pCM, rareCM.rareGn
phenotype = "features" # features, pca
reg.model = "baseline" # baseline, interaction
inter.var = "Metadata_iPSC_Origin" # Metadata_Disease_Status, Metadata_iPSC_Origin
gene.freq = 0.2
print(paste("rare variant burden for", round, ", on phenotype", phenotype, ", aggregate by", sum.or.present, 
            ", variants effect", variants, ", model", reg.model, ", level", level, ", variant from", variant.source,
            ", interaction with", inter.var, ", gene with variants", gene.freq))
suffix = paste(round,level,phenotype,sum.or.present,gene.freq,seed,reg.model,variant.source,"wellMean",sep = ".")
#------------------------------------------------------------
#--------------------- pre-process data ---------------------
#------------------------------------------------------------
# plate check
if(F){
sampled_cells <- readr::read_csv("cmQTLplate7-7-22-20_sampled.csv.gz")

no_check_na_features <- 
  c("TableNumber", "ImageNumber", 
    "Cells_AreaShape_Area",
    "Cytoplasm_AreaShape_Area",
    "Nuclei_AreaShape_Area")

check_na_features <- 
  setdiff(
    names(sampled_cells),
    no_check_na_features
  )

na_frequency <- 
  sampled_cells %>% 
  summarize_at(check_na_features, ~sum(is.na(.))) %>%
  pivot_longer(everything(), values_to = "number_of_na")

na_frequency %>%
  arrange(desc(number_of_na)) %>%
  filter(number_of_na >= 0) %>%
  readr::write_csv("plate7_na_features.csv")

na_frequency %>%
  filter(number_of_na >= 10) %>%
  arrange(desc(number_of_na)) %>%
  show_table %>% head(20)

na_flag <- sampled_cells %>% 
  mutate_at(check_na_features, ~as.double(is.na(.)))

na_flag %>% 
  summarize_at(check_na_features, sum) %>%
  pivot_longer(everything(), values_to = "number_of_na") %>%
  filter(number_of_na >= 10) %>%
  arrange(desc(number_of_na)) %>% 
  filter(!str_detect(name, "_Correlation_")) %>%
  show_table %>% head(20)

data_matrix <- na_flag %>%
  mutate(Nuclei_Correlation_Costes_AGP_Mito = scale(Nuclei_Correlation_Costes_AGP_Mito)) %>%
  select(one_of(c(no_check_na_features, "Nuclei_Correlation_Costes_AGP_Mito"))) %>%
  select(-TableNumber, -ImageNumber)

model <- lm(Nuclei_Correlation_Costes_AGP_Mito ~ ., data_matrix)
summary(model)

# intermediate sampled data from plate 7
data.p7 = readRDS("rds.objects/p7.1.sub.rds")
dim(data.p7)

row.nas = apply(data.p7, 1, function(x){length(which(is.na(x)))})
col.nas = apply(data.p7, 2, function(x){length(which(is.na(x)))})

# --- temp ---
data.p7$NA_count = data.p7 %>% rowwise() %>% is.na() %>% rowSums()
data.p7 %>% ggplot(aes(NA_count)) + geom_histogram(bins = 50)
data.p7 %>% dplyr::select(Metadata_Plate,Metadata_Well,Metadata_line_ID,NA_count) %>% readr::write_csv(path = "rds.objects/p7.1.sub.csv")
# --- temp ---

data.p7 %>% filter(row_number() %in% which(row.nas > 0)) %>% pull(Metadata_Well) %>% table() %>% sort() %>% data.frame() %>%
  dplyr::rename("Line"=".") %>% arrange(Freq) %>% filter(Freq > 5) %>%
  ggplot(aes(Line,Freq)) + geom_col(width = 0.5) + labs(x = "well", y = "NA count")

data.frame("nas" = row.nas) %>% ggplot(aes(nas)) + geom_histogram(bins = 50) + xlab("number of NAs") + theme_bw() + ggtitle("across cells")
data.frame("nas" = col.nas) %>% filter(nas > 100) %>% nrow()

}
# plate summary of processed cells by cell group
if(F){
# select file
file = list.files(path = "/data/srlab/jatin/cmqtl/qced.data", pattern = "*.processed.rds", full.names = T)[seed]

# load processed data
data.sampled = readRDS(file)
print(paste("processed single cell level data loaded for", file))
print(dim(data.sampled))

# # summary
# data.sampled %>% dplyr::select(Metadata_line_ID,Cells_Neighbors_NumberOfNeighbors_Adjacent) %>% 
#   mutate( "cell.cat" = ifelse(Cells_Neighbors_NumberOfNeighbors_Adjacent == 0, "isolate", 
#                               ifelse(Cells_Neighbors_NumberOfNeighbors_Adjacent %in% 1:3, "inter", "colony")) ) %>% 
#   group_by(Metadata_line_ID,cell.cat) %>% count() %>% 
#   readr::write_delim(path=paste0("metadata/cell.category.", sub(".*[/](.*)[.]processed.rds", "\\1",  file), ".summ"), delim="\t")

# summary
data.sampled %>% dplyr::select(Metadata_line_ID,Cells_Neighbors_NumberOfNeighbors_Adjacent) %>% 
  mutate( "cell.cat" = ifelse(Cells_Neighbors_NumberOfNeighbors_Adjacent == 0, "isolate", "non_isolate") ) %>% 
  group_by(Metadata_line_ID,cell.cat) %>% count() %>% 
  readr::write_delim(path=paste0("metadata/cell.cat.", sub(".*[/](.*)[.]processed.rds", "\\1",  file), ".summ"), delim="\t")

print("quitting...")
quit()
}
# preprocess single cell data for each plate
if(F){
# go over sqlite files
files = list.files(path = "/data/srlab/jatin/cmqtl/sqlite", pattern = "*.sqlite", full.names = T)
print(paste("loading and processing data for", files[seed]))
 
# call function
load_sql_data(file = files[seed], good.wells = filled.wells.2)

}
# subset cell lines from processed data
if(F){
  # subset lines per plate
  if(F){
  # get a file to be worked on
  print(paste("doing cell lines sampling for", round))
  file = list.files(path = "rds.objects/singleCell", pattern = "*.qc_ed.rds", full.names = T)[seed]
  
  # load qc-ed data
  data.raw = readRDS(file)
  print(paste("processed data loaded for", file))
  print(dim(data.raw))
  
  # subset data for 2 lines
  data.raw.subset = subset(x = data.raw, subset = Metadata_line_ID %in% sample(unique(data.raw$Metadata_line_ID), 2))
  print(table(data.raw.subset$Metadata_line_ID))
  
  # save rds object of subset data
  filename = paste0("rds.objects/singleCell/subset.data/", sub(".*[/](.*)[.]processed.rds", "\\1",  file), ".", round, ".2lines.rds")
  saveRDS(object = data.raw.subset, file = filename)
  print("rds of qc-ed data for 2 lines saved")
  }
  # combine across plates
  if(F){
  # yo start it
  data.sampled = NULL
  for(file in list.files(path = "rds.objects/singleCell/subset.data", pattern = ".2lines.rds", full.names = T)){
    # load data
    data.sampled.plate = readRDS(file)
    print(paste("sampled lines loaded for", file))
    
    # to have common columns with previous plate
    if(!is.null(data.sampled)){
      common.cols = intersect(colnames(data.sampled), colnames(data.sampled.plate))
      data.sampled = rbind(data.sampled[, common.cols], data.sampled.plate[, common.cols])
    } else data.sampled = data.sampled.plate
    
    print("sampled lines appended")
    print(dim(data.sampled))
    
    rm(list = c("data.sampled.plate"))
  }
  
  saveRDS(object = data.sampled, file = paste0("rds.objects/singleCell/subset.data/", round, ".12lines.rds"))
  print(paste("rds saved for combined sampled lines for", round))
  print(dim(data.sampled))
  
  }
  # gaussianize data
  if(T){
  data.sampled = readRDS(paste0("rds.objects/singleCell/subset.data/", round, ".12lines.rds"))
  print("data for 12 lines loaded")
  print(dim(data.sampled))
  print(table(data.sampled$Metadata_line_ID))
  
  # remove non-variant features
  ftrs = data.sampled[, grep("Metadata", colnames(data.sampled), invert = T, value = T)]
  nzv = nearZeroVar(x = ftrs, saveMetrics = F, freqCut = 90/10, allowParallel = T, names = T)
  print(paste("number of near zero variance features", length(nzv)))
  data.sampled = data.sampled %>% dplyr::select(-one_of(nzv))
  print("near zero variance features removed")
  print(dim(data.sampled))
  
  # gaussianize each feature
  # --- per plate ---
  # data.sampled.gaussian = data.sampled %>% group_by(Metadata_Plate) %>% mutate_at(grep("Meta", colnames(data.sampled), value = T, invert = T), RNOmni::rankNorm)
  # --- across plates ---
  data.sampled.gaussian = data.sampled %>% mutate_at(grep("Meta", colnames(data.sampled), value = T, invert = T), RNOmni::rankNorm)
  
  print("features have been gaussised")
  print(dim(data.sampled.gaussian))
  
  # combine donor information
  data.sampled.gaussian = merge(x = data.sampled.gaussian, y = donor.meta, by.x = "Metadata_line_ID", by.y = "Metadata_donor_Project.Alias", sort = F)
  print("donor info added to data with gaussianized features")
  print(dim(data.sampled.gaussian))
  print(paste("number of cell lines", length(unique(data.sampled.gaussian$Metadata_line_ID))))
  
  ####################################
  # save data for samira
  saveRDS(object = data.sampled.gaussian, file = paste0("rds.objects/singleCell/subset.data/", round, ".12lines.gaussianAcrossPlate.rds"))
  print("gaussianized data saved")
  ####################################
  
  }
}
# sample cells from whole data
if(F){
# get a file to work on
print(paste("doing cell sampling for round", round))
file = list.files(path = "/data/srlab/jatin/cmqtl/qced.data", pattern = "*.processed.rds", full.names = T)[seed]
print(file)

# load qc-ed data
data.raw = readRDS(file)
print(paste("processed data loaded for", file))
print(dim(data.raw))

# take cells which are part of colony or isolate
if(grepl("intcol", round)){
  data.raw = data.raw[data.raw$Cells_Neighbors_NumberOfNeighbors_Adjacent > 0, ]
  print("taking cells from colonies")
  print(dim(data.raw))
} else if(grepl("isolate", round)){
  data.raw = data.raw[data.raw$Cells_Neighbors_NumberOfNeighbors_Adjacent == 0, ]
  print("taking cells from isolate")
  print(dim(data.raw))
} else {
  print("taking any cells")
}

# sample cells from each well
data.sampled.plate = NULL
for(well in unique(data.raw$Metadata_Well)){
  if(nrow(data.raw[data.raw$Metadata_Well %in% well, ]) <= 63) sample.size=nrow(data.raw[data.raw$Metadata_Well %in% well, ]) else sample.size=63
  data.sampled.plate = rbind(data.sampled.plate, data.raw[data.raw$Metadata_Well %in% well, ] %>% sample_n(., size = sample.size))
}
print("sampling from wells done")
print(dim(data.sampled.plate))

# save rds object of sampled data
filename = paste0("/data/srlab/jatin/cmqtl/sampled/sampled_single_cells.", sub(".*[/](.*)[.]processed.rds", "\\1",  file), ".", round, ".rds")
saveRDS(object = data.sampled.plate, file = filename)
print("rds of sampled data saved")

print("quitting...")
quit()
}
# combine sampled data from all plates
if(F){
# iterate over each plate
# combine summarized from plates

###
# all.data.files = list.files("/data/srlab/jatin/cmqtl/sampled_single_cells", pattern = paste0(round, ".rds"), full.names = T)
# data.sampled.plate = all.data.files %>% purrr::map_df(readRDS)
# ###
# index = purrr::map_lgl(data.sampled.plate, ~ any(is.na(.)))
# data.sampled.plate = data.sampled.plate %>% dplyr::select(-one_of(names(which(index == T))))
# print("sampled single cell loaded for all plates")
# print(dim(data.sampled.plate))

# data.sampled = NULL
# for(file in list.files("/data/srlab/jatin/cmqtl/sampled", pattern = paste0(round,".rds"), full.names=T)){
#   # load data
#   data.sampled.plate = readRDS(file)
#   print(paste("sampled data loaded for", file))
#   
#   # to have common columns with previous plate
#   if(!is.null(data.sampled)){
#     common.cols = intersect(colnames(data.sampled), colnames(data.sampled.plate))
#     data.sampled = rbind(data.sampled[, common.cols], data.sampled.plate[, common.cols])
#   } else data.sampled = data.sampled.plate
#   
#   print("sampled data appended")
#   print(dim(data.sampled))
#   rm(list = c("data.sampled.plate"))
# }

# ---- temp ----
# combine sampled cells from all plates
data.files = list.files("/data/srlab/jatin/cmqtl/sampled", pattern = paste0(round,".rds"), full.names=T)
print(length(data.files))
data.sampled = data.files %>% purrr::map_df(readRDS)
print(dim(data.sampled))

# take common columns across all plates
index = purrr::map_lgl(data.sampled, ~ any(is.na(.)))
data.sampled = data.sampled %>% dplyr::select(-one_of(names(which(index == T))))

print("sampled data from all plates combined")
print(dim(data.sampled))
# ---- temp ----

saveRDS(object = data.sampled, file = paste0("rds.objects/singleCell/sampled.data/", round, ".sampled.rds"))
print(paste("rds saved for combined sampled data for", round))

print("quitting...")
quit()
}
# summarize data from single cell to bulk level
if(F){
# get the correct file
print(paste("converting single cells to", level, "level for", round))
file = list.files(path = "/data/srlab/jatin/cmqtl/qced.data", pattern = "*.processed.rds", full.names = T)[seed]

# load processed data
data.sampled = readRDS(file)
print(paste("processed single cell level data loaded for", file))
print(dim(data.sampled))

# take isolate cells
if(round == "isolate") data.sampled = data.sampled %>% filter(Cells_Neighbors_NumberOfNeighbors_Adjacent == 0)
if(round == "intcol") data.sampled = data.sampled %>% filter(Cells_Neighbors_NumberOfNeighbors_Adjacent > 0)
if(round == "colony") data.sampled = data.sampled %>% filter(Cells_Neighbors_NumberOfNeighbors_Adjacent >= 4)
if(round == "intermediate") data.sampled = data.sampled %>% filter(Cells_Neighbors_NumberOfNeighbors_Adjacent > 0 & Cells_Neighbors_NumberOfNeighbors_Adjacent < 4)
print(head(data.sampled$Cells_Neighbors_NumberOfNeighbors_Adjacent))

# make bulk data from single cells
if(level == "bulk") data.sampled.bulk = as.data.frame(data.sampled %>% group_by(Metadata_Plate, Metadata_line_ID) %>% summarise_at(grep("Meta", colnames(data.sampled), invert=T, value=T), mean))
if(level == "well") data.sampled.bulk = as.data.frame(data.sampled %>% group_by(Metadata_Plate, Metadata_Well, Metadata_line_ID) %>% summarise_at(grep("Meta", colnames(data.sampled), invert=T, value=T), mean))
print(dim(data.sampled.bulk))
print(paste("number of cell lines in bulk is", length(unique(data.sampled.bulk$Metadata_line_ID))))

# save bulk data
saveRDS(object = data.sampled.bulk, paste0("rds.objects/",level,"/",round,"_", sub(".*[/](.*)[.]processed.rds", "\\1",  file),"_donor.rds"))
print(paste(level,"level data rds saved"))

print("quitting...")
quit()
}
#------------------------------------------------------------
# load bulk data and get features to test
if(T){
# combine summarized from plates
# if(level == "single") all.data.files = list.files(path = "/data/srlab/jatin/cmqtl/sampled_single_cells", pattern = paste0("anycells.rds"), full.names = T)
if(level != "single") all.data.files = list.files(path = paste0("rds.objects/well/archive"), pattern = paste0(round,".*donor.rds"), full.names = T)
data.bulk = all.data.files %>% purrr::map_df(readRDS)
index = purrr::map_lgl(data.bulk, ~ any(is.na(.)))
data.bulk = data.bulk %>% dplyr::select(-one_of(names(which(index == T))))
print("raw well-level qc-ed data loaded from all plates")
print(dim(data.bulk))

# combine donor information
data.bulk = merge(x = data.bulk, y = donor.meta, by.x="Metadata_line_ID", by.y="Metadata_Project.Alias", sort=F)
print("donor info added to data with all features")
print(dim(data.bulk))

# gaussianize each bulk feature across all plates
data.sampled.bulk.gaussian = data.bulk %>% mutate_at(.vars = grep("Meta|Sample", colnames(data.bulk), value=T, ignore.case=T, invert=T), .funs = RNOmni::RankNorm)
print("features gaussianized")
print(dim(data.sampled.bulk.gaussian))

# remove certain columns
cols.to.remove = grep("sample_|Sample_|_Children_|_Number_|_Parent_|_Center_|_Location_|_Count_|Granularity_1[4-6]|_Orientation|Euler", colnames(data.sampled.bulk.gaussian), value = T)
data.sampled.bulk.gaussian = data.sampled.bulk.gaussian %>% dplyr::select(-one_of(cols.to.remove))
print("certain columns removed")
print(dim(data.sampled.bulk.gaussian))
print(head(data.sampled.bulk.gaussian$Cells_Neighbors_NumberOfNeighbors_Adjacent))

# remove features with near zero variance
ftrs = data.sampled.bulk.gaussian[, grep("Meta", colnames(data.sampled.bulk.gaussian), value = T, invert = T)]
nzv = nearZeroVar(x = ftrs, saveMetrics = F, freqCut = 90/10, allowParallel = T, names = T)
data.sampled.bulk.gaussian = data.sampled.bulk.gaussian %>% dplyr::select(-one_of(nzv))
print("near zero variance features removed")
print(dim(data.sampled.bulk.gaussian))

# on plate edge or not
data.sampled.bulk.gaussian = data.sampled.bulk.gaussian %>% 
  mutate(Metadata_onEdge = sub("[0-9]+", "", Metadata_Well) %in% c("A","P") + sub("[A-Z]+", "", Metadata_Well) %in% c("01","24")) %>%
  mutate(Metadata_onEdge = ifelse(Metadata_onEdge == 0, 0, 1))
print("cell lines on plate edge marked")
print(dim(data.sampled.bulk.gaussian))

grep("^(Cells|Cytoplasm|Nuclei)_",colnames(data.sampled.bulk.gaussian), value = T, invert = F) %>% grep("Neighbo", ., value = T, invert = T) %>%
  data.frame("feature"=.) %>% readr::write_delim("features.all.tab")
  
# summarize gaussianised data from well to donor level, Metadata_cellCat
if(level != "well"){
  data.sampled.bulk.gaussian = data.sampled.bulk.gaussian %>% group_by(Metadata_line_ID,Metadata_Plate) %>% dplyr::select(-Metadata_Well) %>%
    mutate_at(.vars=c("Metadata_onEdge",grep("Meta", colnames(data.sampled.bulk.gaussian), value=T, invert=T, ignore.case=T)), .funs=mean) %>% 
    unique() %>% arrange(Metadata_line_ID) %>% as.data.frame()
  print("features summarized from well to donor level using mean")
  print(dim(data.sampled.bulk.gaussian))
}

# get uncorrelated features for burden testing
if(F){  
  
  ftrs = data.sampled.bulk.gaussian %>% dplyr::select(one_of( replicate_correlation_values %>% filter(median > 0 & !grepl("Meta|_Neighbors_|Sample_|_Orientation", variable)) %>% pull(variable) ))
  assign(x=paste0("dat.", round), value=cor(ftrs))

  # a = dat.isolate
  # b = dat.intcol
  # c = dat.intermediate
  # dat.isolate = a
  # dat.intcol = b
  # dat.intermediate = c

  # get correlated traits
  get_correlation_summ <- function(cor.mat, round){
    cor.stat = apply(cor.mat, 1, function(x){length(which(abs(x) > 0.9))-1})
    tmp = sort(cor.stat, decreasing=T) %>% as.data.frame() %>% tibble::rownames_to_column("feat")
    colnames(tmp) = c("feat", round)
    return(tmp)
  }
  # remove correlated traits from cor mat
  update_cor_mat <- function(cor.mat, feat){
    return(cor.mat[names(which(abs(cor.mat[feat, ]) < 0.9)), names(which(abs(cor.mat[feat, ]) < 0.9))])
  }

  features.to.test = NULL
  while(T){
    # get correlation stats per feature
    iso = get_correlation_summ(dat.isolate, "isolate")
    col = get_correlation_summ(dat.intcol, "colony")
    #int = get_correlation_summ(dat.intermediate, "intermediate")
    
    # get most correlated feature
    feat = iso %>% inner_join(col, by="feat") %>%
      #inner_join(int, by="feat") %>%
      mutate("total"=isolate+colony) %>%
      arrange(-total) %>% head(1) %>% pull(feat)
    if(length(feat) == 0) break
    features.to.test = c(features.to.test, feat)

    # update correlation matrix
    dat.isolate = update_cor_mat(cor.mat=dat.isolate, feat)
    dat.intcol = update_cor_mat(cor.mat=dat.intcol, feat)
    #dat.colony = update_cor_mat(cor.mat=dat.colony, feat)
    #dat.intermediate = update_cor_mat(cor.mat=dat.intermediate, feat)
    print(length(features.to.test))
  }

  # write traits to a file
  features.to.test %>% data.frame("feature"=.) %>% readr::write_delim(path="features.uncorrelated.common2groups.tab")

  # correlation between correlated traits between cell groups
  iso = get_correlation_summ(dat.isolate, "isolate")
  col = get_correlation_summ(dat.intcol, "colony")
  int = get_correlation_summ(dat.intermediate, "intermediate")
  d = iso %>% inner_join(col, by="feat") %>% inner_join(int, by="feat") %>% mutate("total"=isolate+colony+intermediate) %>%
    arrange(-total)

  ###
  cr = cor.test(d$isolate, d$colony)$est
  d[sample(1:nrow(d)), ] %>% mutate("feat.cat" = sub(".*_(.*?)_.*", "\\1", feat), "compt" = sub("_.*","", feat)) %>%
    ggplot(aes(isolate, intermediate)) + geom_point(alpha=0.5, aes(color=feat.cat)) + theme_bw() +
    annotate(geom="text", x=Inf, y=Inf, label=paste("r =",round(cr,2)), hjust=1.5, vjust=2, color="red")

  ###
  a[features.to.test,features.to.test] %>% data.frame() %>% summarise_at(.vars=features.to.test, .funs=median) %>% t() %>% data.frame("r"=.) %>%
    ggplot(aes(r)) + geom_density() + theme_bw() + scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + ggtitle("isolate")

  # histogram of pairwise correlation among selected test
  mat = dat.intcol[features.to.test,features.to.test]
  hist(mat[lower.tri(mat)], breaks=30, main="intermediate")
  
  # get correlated traits
  dat.isolate %>% as.data.frame() %>% dplyr::select(Cells_Intensity_MaxIntensityEdge_DNA) %>% 
    filter(Cells_Intensity_MaxIntensityEdge_DNA > 0.9)
  
}

# colors for feature category
feat.cat = unique(sub(".*_(.*?)_.*", "\\1", grep("Meta",colnames(data.sampled.bulk.gaussian),value = T,invert = T)))
feat.cat.cols = my.cols[1:length(feat.cat)]
names(feat.cat.cols) = feat.cat

# # add image batch information
# data.sampled.bulk.gaussian = data.sampled.bulk.gaussian %>% mutate(Metadata_image_batch = ifelse(Metadata_Plate %in% c("BR00106708","BR00106709"), "a", 
#                                                                                                  ifelse(Metadata_Plate %in% c("BR00107339","BR00107338"), "b", 
#                                                                                                         ifelse(Metadata_Plate %in% c("cmQTLplate7-7-22-20"), "d", "c"))) )
# 
# print("metadata image batch added")
# print(dim(data.sampled.bulk.gaussian))

}
# pca on feature values
if(F){
# PCA to group features
print(paste("doing pca at", level, "for", round, "..."))
if(level == "bulk"){
  pca.res = data.for.pca %>% mutate("Meta_id" = paste(Metadata_Plate,Metadata_line_ID,sep = "_")) %>% tibble::column_to_rownames("Meta_id") %>%
    dplyr::select(one_of(grep("Meta|Sample", colnames(data.sampled.bulk.gaussian), invert = T, value = T))) %>% prcomp()
} else if(level == "well"){
  data.for.pca = data.sampled.bulk.gaussian %>% mutate("Meta_id" = paste(Metadata_Plate,Metadata_Well,Metadata_line_ID,sep = "_"))
  pca.res = data.for.pca[!duplicated(data.for.pca$Meta_id), ] %>% 
    tibble::remove_rownames() %>% tibble::column_to_rownames("Meta_id") %>%
    dplyr::select(one_of(grep("Meta|Sample", colnames(data.for.pca), invert=T, value=T))) %>% prcomp()
}

# PCA embeddings and appending donor metadata
pca.embeds = pca.res$x %>% as.data.frame() %>% dplyr::select(one_of(paste0("PC",1:10))) %>% setNames(object = ., nm = paste0("feat_PC",1:10)) %>%
  tibble::rownames_to_column("Meta_id") %>% mutate("Metadata_line_ID" = sub(".*_", "", Meta_id), Metadata_Plate = sub("_.*", "", Meta_id), Meta_id = NULL) %>%
  mutate("Metadata_line_ID" = as.numeric(Metadata_line_ID)) %>% 
  inner_join(donor.meta, by = c("Metadata_line_ID"="Metadata_Project.Alias"))
  # mutate("Metadata_Age" = as.numeric(Metadata_Age)) %>%
  # inner_join(data.bulk %>% dplyr::select(Metadata_line_ID,Cells_Neighbors_NumberOfNeighbors_Adjacent) %>% unique(), by="Metadata_line_ID")
  # inner_join(cell.count.summ, by = "Metadata_line_ID")

# plot pca
plots.meta = NULL
for(meta in c("Metadata_Sex","Metadata_iPSC_Origin","Metadata_Disease_Status","Metadata_Plate","Metadata_Disease_Category","Cells_Neighbors_NumberOfNeighbors_Adjacent")){
  p = pca.embeds %>% ggplot(aes(feat_PC1, feat_PC2)) + geom_point(aes(color = !!as.name(meta))) + theme_bw() + ggtitle(meta)
  if(!grepl("Neigh", meta)){
    plots.meta[[meta]] = p + ggthemes::scale_color_tableau(name="")
  } else {
    plots.meta[[meta]] = p + ggthemes::scale_color_gradient_tableau(name="")
  }
}
cowplot::plot_grid(plotlist = plots.meta, nrow = 3)


# umap from pca embeddings
umap.defaults$input = "data"
pca.2.umap = umap(d = pca.embeds %>% dplyr::select(paste0("feat_PC",1:10)), config = umap.defaults)
colnames(pca.2.umap$layout) = paste("UMAP", 1:2, sep = "")

umap.embeds = cbind(pca.2.umap$layout, pca.embeds %>% dplyr::select(grep(pattern="Meta", colnames(pca.embeds), value=T)))
umap.embeds %>% ggplot(aes(UMAP1, UMAP2)) + geom_point(aes(color = factor(!!as.name(meta)))) + theme_bw() + ggtitle(meta) +
  guides(color=F)

#filter(!grepl("not", Metadata_Doubling_Hour, ignore.case = T)) %>%
#mutate("Metadata_Doubling_Hour" = as.numeric(Metadata_Doubling_Hour)) %>%
}
# summary of features and donor metadata
if(F){
# summarize correlation between traits
ftrs = data.sampled.bulk.gaussian %>% dplyr::select(one_of( replicate_correlation_values %>% filter(median > 0 & !grepl("Meta|_Neighbors_|Sample_|_Orientation", variable)) %>% pull(variable) ))
res = NULL
for(cor.cut in seq(0.5,0.9,0.1)){
  # take uncorrelated features for burden testing
  cv = findCorrelation(x = cor(ftrs), cutoff = cor.cut, names = T, exact = F)
  features.count = data.frame(feature = setdiff(colnames(ftrs), cv)) %>% 
    filter(!grepl("Meta|_Neighbors_|Sample_|_Orientation", feature)) %>% pull(feature) %>% length()
  res = rbind(res, data.frame("cor" = cor.cut, "nfeat" = features.count))
}
res = rbind(res, data.frame("cor" = 1, "nfeat" = ncol(ftrs)))
ggplot(res, aes(factor(cor), nfeat)) + geom_bar(stat = "identity", width = 0.35, fill = "blue", alpha=0.8) + theme_bw() + 
  theme(panel.grid.major.x = element_blank()) + geom_text(aes(label=nfeat), vjust=-0.5) + xlab("correlation (r)") + 
  ylab("number of traits") + ggtitle(ifelse(round=="intcol", "Colony cells", "Isolated cells")) + 
  scale_y_continuous(expand = expansion(mult = c(0,.1)))

# ----- comparison of variance across cell bins -----  
# cells in colony
data.col = list.files(path = paste0("rds.objects/",level,"/archive"), pattern = paste0("intcol.*donor.rds"), full.names = T) %>% 
  purrr::map_df(readRDS)
index = purrr::map_lgl(data.col, ~ any(is.na(.)))
data.col = data.col %>% dplyr::select(-one_of(names(which(index == T))))
print(dim(data.col))

# cells in others
data.iso = list.files(path = paste0("rds.objects/",level,"/archive"), pattern = paste0("isolate.*donor.rds"), full.names = T) %>% 
  purrr::map_df(readRDS)
index = purrr::map_lgl(data.iso, ~ any(is.na(.)))
data.iso = data.iso %>% dplyr::select(-one_of(names(which(index == T))))
print(dim(data.iso))

# variance per feature
a = data.col %>% summarise_at(.vars=grep("AreaShape|Granularity|Intensity|Texture|RadialDistribution",colnames(data.col), value=T), .funs=var) %>% 
  t() %>% data.frame("var"=.) %>% tibble::rownames_to_column("feat")
b = data.iso %>% summarise_at(.vars=grep("AreaShape|Granularity|Intensity|Texture|RadialDistribution",colnames(data.iso), value=T), .funs=var) %>% 
  t() %>% data.frame("var"=.) %>% tibble::rownames_to_column("feat")

# compare variance per feature between cells in colony and isolate, 3.5x11
feats_remove = "sample_|Sample_|_Children_|_Number_|_Parent_|_Center_|_Location_|_Count_|Granularity_1[4-6]|_Orientation|Euler"
frame = a %>% inner_join(b, by="feat", suffix = c(".col", ".iso")) %>% filter(var.col != 0 & var.iso != 0) %>% 
  filter(!grepl(feats_remove,feat)) %>%
  mutate("feat_comp" = sub("_.*", "", feat), "feat_cat" = sub(".*_(.*?)_.*", "\\1", feat))

# wilcox.test(frame$var.col, frame$var.iso)
# frame %>% dplyr::select(var.col, var.iso) %>% reshape2::melt() %>%
#   ggplot(aes(variable, -log10(value))) + geom_boxplot()

pdf(file="plots/tmp.pdf", width=11, height=3.5)
zoom.in = 1
frame[sample(1:nrow(frame)),] %>% filter(var.col < zoom.in) %>% ggplot(aes(var.col,var.iso)) + geom_point(aes(color = feat_cat)) + theme_bw() + geom_abline(slope=1) + 
  facet_wrap(~feat_comp, scales="free") + ggthemes::scale_color_tableau(name="feature class") + 
  theme(strip.background=element_blank()) +
  labs(x="feature variance in colony cells", y="feature variance in isolate cells", title=zoom.in)
dev.off()

# ----------- summary of donors' metadata -----------
# "Metadata_Disease_Status", "Metadata_Plate", 
meta.plots=NULL
dat = donor.meta %>% dplyr::select(one_of("Metadata_iPSC_Origin","Metadata_Sex","Metadata_Disease_Category"))
meta.plots[["all"]] = reshape2::melt(apply(dat, 2, function(x) {table(x)})) %>% group_by(L1) %>% mutate("color" = paste0("color",1:n())) %>%
  ggplot(aes(sub("Metadata_","",L1), value, fill = color)) + geom_bar(position=position_stack(), stat="identity", width=0.4, alpha=0.85) + 
  theme_bw() + theme(panel.grid.major=element_blank()) + labs(x="Donor's metadata",y="count") + 
  scale_y_continuous(expand = expansion(mult = c(0,.01))) + theme(axis.text=element_text(size=10), axis.title=element_text(size=13)) +
  geom_text(aes(label = x), position=position_stack(vjust=0.5), color="black", size=3) +
  scale_fill_manual(values=rep( ggsci::pal_jco()(10)[3:2], 3), guide=F)

meta.plots[["Metadata_Age"]] = donor.meta %>% dplyr::select(one_of("Metadata_Project.Alias","Metadata_Age")) %>% 
  unique() %>% pull(Metadata_Age) %>% as.numeric() %>% 
  data.frame("x" = .) %>% ggplot(aes(x)) + geom_histogram(bins=30, fill=ggsci::pal_jco()(10)[3], alpha=0.75) + 
  theme_bw() + labs(x = "Donor's age") + scale_y_continuous(expand = expansion(mult = c(0,.01))) +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=13))

# plot
cowplot::plot_grid(plotlist = meta.plots, nrow = 1, rel_widths=c(1.25,1))
# pdf("plots/donor.meta.pdf", width=13, height=2.25)
# print(comb.p)
# dev.off()

# ---- ----------- other summaries ------------------
# gaussian distribution of a feature
data.sampled.bulk.gaussian %>% ggplot(aes(Cells_AreaShape_Eccentricity)) + geom_histogram(bins=50, aes(fill=..count..)) + theme_bw() + 
  ggthemes::scale_fill_gradient_tableau(guide=F) + scale_y_continuous(expand = expansion(mult = c(0,.1)))

# summary of feature categories
col.map = c(feat.cat.cols, setNames(object=rep(ggsci::pal_jco()(10)[3],3), nm=c("Cells","Cytoplasm","Nuclei")))
to.plot = "feat.cat"
# colnames(data.sampled.bulk.gaussian)
p = features.to.test %>% data.frame("feature"=.) %>%
  filter(!grepl("Meta|_Neighbors_|Sample_|_Orientation", feature)) %>%
  mutate("feat.cat" = sub(".*_(.*?)_.*","\\1", feature), "compt" = sub("_.*","", feature)) %>% 
  group_by(!!as.name(to.plot)) %>% count() %>%
  ggplot(aes(!!as.name(to.plot), n)) + labs(x="Category", y="Number of traits") + coord_flip() +
  theme_bw() + geom_bar(stat = "identity", position = position_dodge2(), width = 0.3, aes(fill = !!as.name(to.plot))) + 
  geom_text(aes(label=n), vjust=0.5, size=4.5, hjust=-.1) + scale_fill_manual(values=col.map, name="trait class", guide=F) +
  scale_y_continuous(expand = expansion(mult=c(0,.2))) + scale_x_discrete(limits = rev) +
  theme(axis.text = element_text(size=11), axis.title=element_text(size=12), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
png(filename = "plots/traits.selected.summ.categories.png", width = 4.5, height = 3, units = "in", res = 300)
print(p)
dev.off()

}
#------------------------------------------------------------
#--------------------- variance component -------------------
#------------------------------------------------------------
# variance component analysis
if(F){
cores=10
print(paste("doing variance component analysis done for", round, "at", level, "level", "for seed", seed))

# features to test, features.uncor
features.to.test = colnames(data.sampled.bulk.gaussian) %>% data.frame("feature"=.) %>%
  filter(!grepl("Meta|_Neighbors_|Sample_|_Orientation", feature)) %>% pull(feature)
print(paste(length(features.to.test), "features to test for variance component analysis"))

# parallelize variance component analysis
block = split(1:length(features.to.test), cut_interval(1:length(features.to.test), n = cores))
feature.this.block = features.to.test[block[[seed]]]
print(paste("testing block of",length(feature.this.block), "features"))

# merge genotype PCs with traits
data.allPlates = data.sampled.bulk.gaussian %>% 
  inner_join(pcs.terra, by=c("Metadata_line_ID"="Metadata_Project.Alias")) %>%
  mutate_at(.vars = "Metadata_Age", .funs = RNOmni::RankNorm)
print("gt pcs merged with features and age gaussianized")
print(dim(data.allPlates))

# variance component analysis
var.res = variance_analysis(data.allPlates = data.allPlates, features.to.test = feature.this.block, round = round, level = level)
var.res %>% mutate(p.adj = p.adjust(p=p, method="bon", n=length(features.to.test))) %>% 
  readr::write_delim(paste("varcov/vc", round, level, seed, "allFeatures.nov.2022.tab", sep = "."), delim = "\t")
print("variance component analysis done")

print("quitting...")
quit()
}
# plot results
if(F){
# load var component analysis res
var.res = list.files(path="varcov", pattern="vc.anycells.well.*.nov.2022.tab", full.names=T) %>% purrr::map_df(read.delim) %>%
  mutate("feat.cat" = sub(".*_(.*?)_.*", "\\1", feat), "pvcov.2" = ifelse(p.adj < 0.05, pvcov.2, 0)) %>%
  mutate("xticks" = Hmisc::capitalize(sub("Metadata_","",grp))) %>% 
  mutate(xticks = replace(xticks, xticks == "Line_ID", "Difference_Donors")) %>%
  mutate(xticks = replace(xticks, xticks == "OnEdge", "Plate_Edge")) %>%
  mutate("xticks" = factor(xticks, levels=sort(unique(xticks)))) %>%
  filter(xticks != "Residual")

# difference between component per feature class
a = var.res %>% filter(grp == "Metadata_line_ID" & feat.cat == "AreaShape" & p.adj < 0.05) %>% dplyr::select(-one_of("p.adj","p"))
b = var.res %>% filter(grp == "Metadata_line_ID" & feat.cat != "AreaShape" & p.adj < 0.05) %>% mutate(feat.cat = "Not_AreaShape") %>% dplyr::select(-one_of("p.adj","p"))
test = wilcox.test(a$pvcov,b$pvcov)
rbind(a,b) %>% ggplot(aes(feat.cat,pvcov)) + geom_jitter(width=0.1, alpha=0.1) + geom_boxplot(width=0.5, outlier.shape=NA) + 
  annotate(geom="text", x=Inf, y=Inf, hjust=1.05, vjust=1.5, label=paste("p =", format(test$p.value,sc=T,digits=1))) + 
  theme_bw() + labs(x="trait category", y = "variation in traits accounted by\n difference among donors")

# selected features
feat.labels = c("Cytoplasm_AreaShape_Zernike_9_3","Cytoplasm_Granularity_3_RNA","Cells_RadialDistribution_RadialCV_Mito_1of4")

# plot individual component for selected features
png(filename = "plots/var.comp.selected.png", width = 11.5, height = 2.75, units = "in", res = 300)
p = var.res %>% filter(feat %in% feat.labels) %>%
  ggplot(aes(xticks, pvcov.2)) + geom_col(width = 0.45, aes(fill = feat.cat)) + theme_bw() + 
  ylab("Accounted variance in traits") + xlab("") + facet_wrap(~feat, nrow=1) + coord_flip() + 
  scale_fill_manual(values=feat.cat.cols, guide=F) + scale_x_discrete(limits = rev) +
  theme(strip.background = element_blank(), strip.text=element_text(size=11)) +
  theme(axis.text = element_text(size=11), axis.title=element_text(size=12)) +
  scale_y_continuous(expand = expansion(mult = c(0,.02)))
print(p)
dev.off()

# plot overall variation, box plots
p = var.res %>% mutate("color" = ifelse(p.adj < 0.05, feat.cat, "grey")) %>%
  ggplot(aes(xticks, pvcov)) + 
  #geom_point(position = position_jitterdodge(jitter.width=0.075, dodge.width=0.75), 
  #           alpha=0.75, size=0.25, aes(group=feat.cat, color=color)) +
  geom_boxplot(outlier.shape=NA, width=0.5, aes(color=color), position=position_dodge(preserve="total",width=0.75)) +
  scale_color_manual(values=c(setNames("grey", "grey"), feat.cat.cols), breaks=names(feat.cat.cols)) + 
  labs(x="", y="Accounted variance in traits", color="Trait category") + theme_bw() +
  theme(axis.text = element_text(size=11), axis.title=element_text(size=12)) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), panel.grid.major=element_blank())
png(filename = "plots/var.comp.overall.png", width = 7, height = 4, units = "in", res = 300)
print(p)
dev.off()
  
# residual vs cell line component
pdf("plots/var.cellid.res.pdf", height=3.5, width=5)
p = var.res %>% filter(p.adj < 0.05) %>% group_by(feat) %>% mutate(explained.var = sum(pvcov.2)) %>% 
  filter(grp %in% c("Metadata_line_ID","Residual")) %>% 
  dplyr::select(feat,grp,pvcov.2,feat.cat) %>% reshape2::dcast(., feat + feat.cat ~ grp, value.var="pvcov.2") %>%
  ggplot(aes(Metadata_line_ID,Residual)) + geom_point(aes(color=feat.cat), size=3, alpha=0.6, stroke=0) +
  ggrepel::geom_text_repel(aes(label = ifelse(feat %in% feat.labels, feat, "")), direction="both", nudge_y=0.05, nudge_x=-0.3, size=3, force=3) +
  scale_color_manual(values=feat.cat.cols) + theme_bw() +
  labs(x="Feature variance accounted by\ndifference across cell lines", y="Residual variance", color="Feature class") +
  theme(axis.text = element_text(size=13), axis.title=element_text(size=13))
print(p)
dev.off()

# combine variance component and association analysis
var.res %>% filter(p.adj < 0.05 & feat %in% sig.assoc.p.colony.bulk$feat & grp %in% c("Metadata_line_ID","Residual")) %>% 
  group_by(feat) %>% mutate(explained.var = sum(pvcov.2)) %>% arrange(-explained.var)

# intersect variance components with association testing
var.res %>% filter(grp %in% c("Metadata_line_ID","Residual") & p.adj < 0.05) %>% 
  dplyr::select(grp,pvcov.2,feat) %>% reshape2::dcast(formula=feat~grp, value.var="pvcov.2") %>%
  filter(Metadata_line_ID > 0.4 & Residual > 0.2) %>%
  inner_join(sig.assoc.p.intermediate.bulk, by="feat")

}
#------------------------------------------------------------
#----------------------- extreme values ---------------------
#------------------------------------------------------------
# freq spectrum of rare gnomad variant in cmqtl
if(F){
setwd("~/Documents/cmQTL")
gnomad.in.cmqtl = read.delim("other/cmqtl.gnomad.rare.frq", row.names = NULL)
setwd("/Volumes/ja286/cmQTL")

# plot 5x3
gnomad.in.cmqtl %>% mutate(maf = sub(".*:", "", X.ALLELE.FREQ.)) %>% filter(as.numeric(maf) < 0.25) %>% 
  ggplot(aes(as.numeric(maf))) + geom_histogram(bins = 25) + geom_vline(xintercept = 0.01, color = "blue") + theme_bw() + 
  ggtitle('freq. of rare gnomad variants\nin cmqt dataset')

###
cmqtl.gnomadAnn = read.delim("wgs/rareCM.GnomadAnn.tab")
data.to.plot = cmqtl.gnomadAnn[complete.cases(cmqtl.gnomadAnn), ] %>% filter(AF <= 0.01) %>% dplyr::select(one_of("AF","AF_amr","AF_afr","AF_nfe","AF_eas","AF_sas","AF_ami"))

data.to.plot %>% melt(data = ., id.vars = "AF") %>% filter(value > 0.01) %>%
  ggplot(aes(value)) + geom_histogram(bins = 50, aes(fill = variable)) + ggthemes::scale_fill_tableau()

data.to.plot[sample(1:nrow(data.to.plot), size = nrow(data.to.plot)*0.1), ] %>% melt(data = ., id.vars = "AF") %>%
  ggplot(aes(AF, value)) + geom_point(aes(color = variable), alpha = 0.5, size = 1.25) + facet_wrap(~variable) + #geom_smooth(method = "lm") +
  theme_bw() + ggthemes::scale_color_tableau() + labs(x = "maf in cmqtl", y = "maf in gnomad") #+ geom_hline(yintercept = 0.01, color = "red", size = 0.1)

#pdf(file = "plots/freq_cmqtlRare_GnomadPops.pdf", width = 7.5, height = 3.5, compress = T)
#dev.off()
}
# genes to test
if(T){
# load genotype of variants for all donors
donor.gt = read.delim(paste0("wgs/",variant.source,".high.tab"), header=F, strip.white=F, col.names=read.table("bekar.cheeze/tmp.txt")$V1, stringsAsFactors=F) %>% slice(-1)
if(variants == "high") donor.gt = donor.gt %>% filter(IMPACT == "HIGH")
if(variants == "moderate") donor.gt = donor.gt %>% filter(IMPACT == "MODERATE")
if(variants == "missense") donor.gt = donor.gt %>% filter(grepl("missense_variant", EFFECT))
print(paste(nrow(donor.gt), "variants taken with impact", variants))

# recode genotypes to 0,1,2
donor.gt.re = cbind(donor.gt[,c(1:10)], apply(donor.gt[,-c(1:10)], 2, function(x){ as.numeric(sub("0/1|1/0", "1", sub("1/1", "2", sub("0/0|\\./\\.", "0", x)))) }))

# summarize variant count per gene
if(sum.or.present == "sum"){
  donor.gt.re.sum = donor.gt.re %>% group_by(GENE) %>% summarise_at(vars(-one_of(colnames(donor.gt.re)[1:10])), .funs = sum) %>% 
    tibble::column_to_rownames("GENE") %>% as.data.frame() %>% t()
} else if(sum.or.present == "present"){
  myop <- function(x){
    return(ifelse(sum(x) > 0, 1, 0))
  }
  donor.gt.re.sum = donor.gt.re %>% group_by(GENE) %>% summarise_at(vars(-one_of(colnames(donor.gt.re)[1:10])), myop) %>% 
    tibble::column_to_rownames("GENE") %>% as.data.frame() %>% t()
}

# total number of variants
donor.total.vars = donor.gt.re %>% summarise_at(vars(-one_of(colnames(donor.gt.re)[1:10])), sum) %>% t() %>% data.frame("Metadata_total_vars" = .) %>% 
  tibble::rownames_to_column("wgs_id") %>% mutate("wgs_id" = sub("^X", "", toupper(gsub("\\_|\\.", "-", wgs_id))))

# combine donor metadata and genotype
gene.vars = donor.gt.re.sum %>% as.data.frame() %>% tibble::rownames_to_column("wgs_id") %>% mutate("wgs_id" = sub("^X","",toupper(gsub("\\_|\\.","-",wgs_id)))) %>%
  inner_join(donor.meta %>% mutate("wgs_id" = toupper(gsub("\\_|\\.","-",collaborator_sample_id))) %>% dplyr::select(c(wgs_id,Metadata_Project.Alias,Metadata_Disease_Status,!!as.name(inter.var))), by = "wgs_id") %>%
  inner_join(donor.total.vars, by = "wgs_id") %>% dplyr::select(-one_of("wgs_id"))
gene.vars = gene.vars[!duplicated(gene.vars$Metadata_Project.Alias), ]
print("summarized variant count per gene per donor prepared")
print(dim(gene.vars))

# take genes to test sum(.) <= ceiling(nrow(gene.vars)*gene.freq)
if(file.exists("rds.objects/genes.to.test.june2021.rds") & variants == "all" & variant.source == "rareCM"){
  genes.to.test = readRDS("rds.objects/genes.to.test.june2021.rds")
} else if(reg.model == "baseline"){
  genes.to.test = gene.vars %>% group_by(Metadata_iPSC_Origin) %>% 
    summarise_at(grep("Meta|Sample", colnames(.), value=T, invert=T), ~nnzero(.x)) %>%
    select_if( function(.) is.numeric(.) && all(. >= 1) && sum(.) > ceiling(nrow(gene.vars)*0.02) && sum(.) < floor(nrow(gene.vars)*0.98) ) %>% 
    colnames() %>% unique() %>% data.frame("gene" = .) %>% 
    filter(!grepl('\\.[0-9]|^RP[0-9]|^RPL|&|-[0-9]|^PGB|^AP[0-9]{2,}|^MTX|^Meta', gene)) %>% 
    mutate(gene.2 = sub("-", "_", gene)) %>% dplyr::select(gene.2)
}
print(paste("gene to test", nrow(genes.to.test)))

}
# prepare data for variant burden test
if(T){
# select feature values and features to test
features.values = data.sampled.bulk.gaussian
features.to.test = read.delim("features.uncorrelated.common2groups.tab") %>% 
  filter(!grepl("Meta|_Neighbors_|Sample_|_Orientation", feature)) %>% pull(feature) %>% unique()
print(paste("feature values and",length(features.to.test),"features to test taken"))

# remove cell lines which are on >1 plates
if(level == "bulk"){
  features.values = features.values[!duplicated(features.values$Metadata_ID), ]
  print("duplicated lines removed")
}

# # gaussianize age
# features.values = features.values %>% mutate_at(.vars = "Metadata_Age", .funs = RNOmni::RankNorm)
# print("donor age gaussianized")

# add cell count per cell line
if(level == "bulk") features.values = features.values %>% inner_join(cell.count.summ.line, by="Metadata_line_ID") %>%
  mutate_at(.vars = "Metadata_CellCount_avg", .funs = RNOmni::RankNorm)
print("cell count appended and gaussianized")

# gaussianise doubling hour
if(length(grep("Not", features.values$Metadata_Doubling_Hour, ignore.case=T)) == 0) features.values = features.values %>% mutate_at(.vars = "Metadata_Doubling_Hour", .funs = as.numeric)

# overall
print(paste("final feature values has donors", nrow(features.values)))

# split features for parallel processing
if(reg.model == "baseline" & variants == "all" & level %in% c("bulk","multibulk","well")){
  block = split(1:length(features.to.test), cut_interval(1:length(features.to.test), n = cores))
  feature.this.block = features.to.test[block[[seed]]]
  print(paste("testing block of",length(feature.this.block), "features"))

} else if(reg.model == "baseline" & variants == "all" & level == "well"){
  features.to.test = readr::read_delim(paste0("extreme.values/finals/",round,".bulk.tab"), delim="\t") %>% pull(feat) %>% unique()
  cores = length(features.to.test)
  print(paste(length(features.to.test),"features to test"))
  block = split(1:length(features.to.test), cut_interval(1:length(features.to.test), n = cores))
  feature.this.block = features.to.test[block[[seed]]]
  print(paste("testing",length(feature.this.block), "features"))
}

# # prepare data for gwas by Samira
# data = features.values %>% dplyr::select(one_of(grep("Meta",colnames(features.values),value=T), "Cells_Neighbors_NumberOfNeighbors_Adjacent", features.to.test)) %>%
#   inner_join(pcs.terra, by=c("Metadata_line_ID"="Metadata_Project.Alias"))
# saveRDS(object=data, file=paste0("rds.objects/", level, "/", round, ".final.rds"))
# print("data saved for gwas")
# print(dim(data))

}
# summary of variants
if(F){
# ------------------------------------------------------
# localisation of variants of specific genes in variant quality space
donor.gt %>% mutate(marked = ifelse(GENE %in% c("WASF2","ZNF576","PRLR","TSPAN15"),"y","n")) %>% 
    dplyr::select(GENE,MQ,DP,marked) %>% arrange(marked) %>%
    ggplot(aes(as.numeric(MQ),log10(as.numeric(DP)))) + geom_point(aes(color=marked), size=0.5) +
    scale_color_manual(values=c("grey","blue"), guide=F) + theme_bw() + 
    labs(x="mapping quality", y="log10(number of supporting reads)") + ggtitle("quality of variants")
    
# ------------------------------------------------------
# cell count per replicate grouped by wihout and with variants
gene = "TSPAN15"
frame = gene.vars %>% dplyr::select(Metadata_Project.Alias,!!as.name(gene)) %>%
  inner_join(cell.count.summ.well, by =c("Metadata_Project.Alias"="Metadata_line_ID")) %>% 
  mutate( "col" = ifelse(!!as.name(gene) != 0, Metadata_Project.Alias, "wt") )

wtest = wilcox.test(frame %>% filter(col == "wt") %>% pull(Metadata_CellCount_avg), frame %>% filter(col != "wt") %>% pull(Metadata_CellCount_avg))
frame %>% ggplot(aes(factor(!!as.name(gene)), Metadata_CellCount_avg)) + geom_boxplot(width=0.5, outlier.shape=NA) +
  geom_jitter(width=0.2, size=1, alpha=0.75, aes(color=col)) + ggthemes::scale_color_tableau() + theme_bw() +
  annotate(geom="text", x=Inf, y=Inf, label=paste0("p = ",round(wtest$p.value,5)), hjust=1.1, vjust=1.5) +
  ggtitle(gene) + xlab("rare variant(s) presence")
  
# ------------------------------------------------------
# correlation between variant count and gene length
library("biomaRt")
gene.var.count = donor.gt.re %>% group_by(GENE) %>% count() %>% filter(!grepl("-[0-9|A-Z]{5,}[\\.]",GENE)) %>% 
  filter(!grepl("&",GENE))

# get positon of genes
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
attributes = c("start_position","end_position","hgnc_symbol","chromosome_name")
genes.pos = getBM(attributes=attributes, filters=c("hgnc_symbol"), values=list(hgnc_symbol=unique(gene.var.count$GENE)), mart=mart, uniqueRows=T)

# merge variant count with gene length
gene.var.count.all = gene.var.count %>% inner_join(genes.pos, by=c("GENE"="hgnc_symbol")) %>% filter(chromosome_name %in% 1:22) %>% 
  mutate(len = end_position-start_position)

# correlation and plot
cr = cor.test(gene.var.count.all$n, gene.var.count.all$len)
gene.var.count.all %>% #filter(len < 10000) %>%
  ggplot(aes(len,n)) + geom_point(alpha=0.75, size=0.5) + theme_bw() + labs(x="gene length", y="#rare variants") +
  geom_smooth(method="lm", se=F) +
  annotate(geom="text", x=Inf, y=Inf, hjust=1.5, vjust=2, label = paste("r =",round(cr$estimate,2))) + 
  annotate(geom="text", x=Inf, y=Inf, hjust=1.25, vjust=4, label = "p < 0.00001") 

# ------------------------------------------------------
# overall total number of variants (rare+common)
sample.stats = read.delim("bekar.cheeze/SCBB_Aug2020_Sample_Stats_PSC.vchk", comment.char = "#") %>% 
  mutate("snpsTotal" = nNonRefHom+nHets, "wgs_id" = sub("^X", "", toupper(gsub("\\_|\\.", "-", sample)))) %>%
  mutate("wgs_id" = ifelse(wgs_id == "SCBB-1285", "CW20020-P12-DH-12-4-17", wgs_id))

# distribution of variants per donor
donor.gt.impact.tab = donor.gt.re %>% group_by(IMPACT) %>% summarise_at(vars(-one_of(colnames(donor.gt.re)[1:10])), sum) %>% tibble::column_to_rownames("IMPACT") %>% 
  t() %>% data.frame() %>% tibble::rownames_to_column("wgs_id") %>% mutate("wgs_id" = sub("^X", "", toupper(gsub("\\_|\\.", "-", wgs_id)))) %>%
  inner_join(donor.meta %>% mutate("wgs_id" = ifelse(test = grepl("CW20020",to_sample_id), "CW20020-P12-DH-12-4-17", toupper(gsub("^X", "", toupper(gsub("\\_|\\.", "-", collaborator_sample_id)))))), by = "wgs_id") %>% 
  inner_join(sample.stats, by = "wgs_id") %>% dplyr::select(-one_of("wgs_id","collaborator_sample_id","to_sample_id")) %>%
  inner_join(data.sampled.bulk.gaussian %>% dplyr::select(Metadata_ID,Metadata_Plate), by = c("Metadata_ID")) %>%
  inner_join(pcs.terra, by = c("Metadata_Project.Alias"="Metadata_Project.Alias")) %>% 
  mutate("disease" = ifelse(Metadata_Disease_Category == "Control", 0, 1)) %>%
  mutate("TOTAL" = HIGH + MODERATE, Metadata_Age = as.numeric(Metadata_Age))
donor.gt.impact.tab = donor.gt.impact.tab[!duplicated(donor.gt.impact.tab$Metadata_Project.Alias), ]

# samples in PC space
meta = "TOTAL"
p1 = donor.gt.impact.tab %>% ggplot(aes(PC1,PC2)) + geom_point(aes(color = !!as.name(meta))) + theme_bw() + 
  ggtitle(meta) + guides(color=F)
p2 = donor.gt.impact.tab %>% ggplot(aes(!!as.name(meta))) + geom_histogram(bins=50, aes(fill=..count..)) + theme_bw() + 
  ggtitle(meta) + guides(fill=F)
pdf(file="plots/tmp.pdf", width=6, height=3.25)
gridExtra::grid.arrange(p2, p1, nrow = 1)
dev.off()

# distribution of variants per gene
dat = gene.vars %>% summarise_at(grep("Meta|Sample", colnames(.), value=T, invert=T), ~nnzero(.x)) %>% 
  colSums() %>% data.frame("n.vars" = .) %>% mutate("freq" = n.vars/nrow(gene.vars)) %>% mutate(freq = ifelse(freq > 0.5, 1-freq, freq))
p1 = dat %>% ggplot(aes(n.vars)) + geom_histogram(bins=50, aes(fill=..count..)) + theme_bw() + 
  xlab("number of donors with variants") + guides(fill=F)
p2 = dat %>% filter(n.vars < nrow(gene.vars)*0.2) %>% ggplot(aes(n.vars)) + geom_histogram(bins=60, aes(fill=..count..)) + theme_bw() + 
  xlab("number of donors with variants") + guides(fill=F) + scale_x_continuous(breaks=seq(0,nrow(gene.vars)*0.2,5))
gridExtra::grid.arrange(p1, p2, nrow = 1)

# ------------------------------------------------------
# if metadata can predict disease status of donors
plots.list = NULL
plots.pcs.list = NULL
for(meta in c("Metadata_Sex","Metadata_iPSC_Origin","HIGH","MODERATE","TOTAL","snpsTotal")){
  
  print(paste("checking",meta,"..."))
  my.form = as.formula(paste("disease ~", meta, "+ Metadata_Plate +", paste0("PC",1:5,collapse = "+")))
  model.res = glm(formula = my.form, family = binomial(link = "logit"), data = donor.gt.impact.tab)
  p = coefficients(summary(model.res)) %>% data.frame() %>% filter(grepl(meta, rownames(.))) %>% pull(Pr...z..) %>% round(2)
  
  # plot the results in bar or box plot
  if(meta %in% c("HIGH","MODERATE","TOTAL","snpsTotal")){
    plot = ggplot(donor.gt.impact.tab, aes(factor(disease), !!as.name(meta))) + geom_boxplot(outlier.size = 0, width = 0.45, aes(color = factor(disease))) +
      geom_jitter(width = 0.4, size = 1, alpha = 0.5, aes(color = factor(disease))) + ggthemes::scale_color_tableau(guide = F) +
      annotate(geom = "text", x=Inf, y=Inf, size=4, hjust=1.15, vjust=3, label = paste("P =", p))
  } else{
    plot = ggplot(donor.gt.impact.tab, aes(factor(disease))) + geom_bar(stat = "count", width = 0.6, aes(fill = !!as.name(meta))) + 
      ggthemes::scale_fill_tableau(name = sub(".*_","",meta)) + annotate(geom = "text", x=1, y=Inf, size=4, hjust=0.25, vjust=2, label = paste("P =", p))
  }
  plot = plot + theme_bw() + xlab("disease status") + ggtitle(meta)
  plots.list[[meta]] = plot
  
  # distibution of meta in PC space  
  plot.pca = ggplot(donor.gt.impact.tab, aes(PC1,PC2)) + geom_point(size = 2, aes(color = !!as.name(meta))) + theme_bw() + ggtitle(sub(".*_","",meta))
  if(!(meta %in% c("HIGH","MODERATE","TOTAL","snpsTotal"))) plot.pca = plot.pca + ggthemes::scale_color_tableau() else plot.pca = plot.pca + ggthemes::scale_colour_gradient_tableau("Blue-Teal") 
  plots.pcs.list[[meta]] = plot.pca
}

# plots
cowplot::plot_grid(plotlist = plots.list, nrow = 2)
cowplot::plot_grid(plotlist = plots.pcs.list, nrow = 2)

# samples in PC space with 1kg samples
pcs.kg = read.table("extreme.values/gt.1kg.pcs", col.names = c("IID","FID",paste0("PC",1:20)))
pcs.kg.map = read.delim("metadata/GSA_Terrasamples_Meta_Fixed_V2.txt") %>% unique()
pcs.outliers = pcs.kg.map %>% filter(Linking.True.Donor.ID %in% (donor.gt.impact.tab %>% arrange(-TOTAL) %>% head(7) %>% pull("Metadata_ID"))) %>% 
  pull("Chipwell.barcode")

# locate outliers in PC space
pcs.kg[sample(nrow(pcs.kg)),] %>% mutate("color" = ifelse(IID %in% pcs.outliers, "outlier", "normal")) %>% 
  mutate("label" = ifelse(IID %in% pcs.outliers, "outlier", "")) %>%
  ggplot(aes(PC1,PC2)) + geom_point(aes(color = color)) + theme_bw() + scale_color_manual(values = c("grey","red")) + 
  ggrepel::geom_text_repel(aes(label = label))

}
# do rare variant burden test
if(T){
# iterate over block of features broke in previous section
burden.res = NULL
for(feature in feature.this.block){
  # start
  print(paste("calculating variant burden for", feature))
  
  # covariates - the cell cycle is Nuclei_Intensity_IntegratedIntensity_DNA, Metadata_image_batch, Metadata_CellCount_avg
  # "Metadata_iPSC_Origin", "Metadata_Disease_Status", "Metadata_Age", "Metadata_cellCat", Metadata_CellCount_avg
  covars = c("Metadata_Plate", "Metadata_Sex", paste0("PC",1:4),"Metadata_onEdge")
  if(level == "well") covars = union(covars, c("Metadata_Well"))
  if(round != "isolate") covars = union("Cells_Neighbors_NumberOfNeighbors_Adjacent", covars)
  
  # subset a given feature, filter cell lines based on wgs samples, and merge with pcs info
  feature.gene.var = features.values %>% dplyr::select(c(Metadata_line_ID, grep("PC[0-9]|total_vars",covars,invert=T,value=T), !!as.name(feature))) %>%
    inner_join(gene.vars %>% dplyr::select(-one_of(inter.var, "Metadata_Disease_Status")), by = c("Metadata_line_ID"="Metadata_Project.Alias")) %>%
    inner_join(pcs.terra, by = c("Metadata_line_ID"="Metadata_Project.Alias"))
  colnames(feature.gene.var) = sub("[-]", "_", colnames(feature.gene.var))
  print(paste("data prepared", nrow(feature.gene.var)))
  
  # original
  # regress_feature_allele(feature = feature, gene = gene, covars = covars, data = feature.gene.var, model = reg.model, level = level, inter.var = inter.var)
  res = do.call("rbind", apply(genes.to.test, 1, function(x){ regress_feature_allele(feature = feature, gene = x, covars = covars, data = feature.gene.var, model = reg.model, level = level, inter.var = inter.var) }))
  res = res %>% mutate(feat = feature, cells = round, variants = sum.or.present, p.adj = p.adjust(p, "bon", nrow(genes.to.test)*length(features.to.test)), data = "orig")
  print("regression done")

  # store results
  burden.res = rbind(burden.res, res)
  print("burden results stored")
}

# save results
if(reg.model == "baseline") write.table(x = burden.res, file = paste("extreme.values/burden.out.11.2022/rv", suffix, "tab", sep = "."), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
if(reg.model == "interaction") write.table(x = burden.res, file = paste("extreme.values/burden.out.11.2022/rv", suffix, inter.var, "tab", sep = "."), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
print(paste("burden test results written for", seed))
}
#---------------------- permutation test --------------------
# permutation test
if(F){
print(paste("calculating variant burden with permutation, round", seed))
# select those traits which were not taken previously
# features.to.test = setdiff(features.to.test, features.uncor$feature)
print(paste("--> testing features", length(features.to.test)))

# permute feature values
my.seed = sample(1:5000000,1)
set.seed(my.seed)
if(round != "isolate") features.values = features.values %>% mutate_at(.vars=c("Metadata_line_ID","Cells_Neighbors_NumberOfNeighbors_Adjacent"), .funs=sample)
if(round == "isolate") features.values = features.values %>% mutate_at(.vars=c("Metadata_line_ID"), .funs=sample)
print("feature matrix permuted")

burden.res = NULL
for(feature in features.to.test){
  # select covariates
  print(feature)
  # covars = c("Metadata_Plate", "Metadata_iPSC_Origin", "Metadata_Sex", "Metadata_Disease_Status", "Metadata_Age", paste0("PC",1:5))
  covars = c("Metadata_Plate", "Metadata_Sex", paste0("PC",1:5))
  if(level == "well") covars = union(covars, c("Metadata_Well", "Metadata_line_ID"))
  if(round != "isolate") covars = union("Cells_Neighbors_NumberOfNeighbors_Adjacent", covars)
  
  # prepare data
  feature.gene.var = features.values %>% dplyr::select(c(Metadata_line_ID, grep("PC[0-9]|total_vars",covars,invert=T,value=T), !!as.name(feature))) %>%
    inner_join(gene.vars %>% dplyr::select(-one_of(inter.var,"Metadata_Disease_Status")), by = c("Metadata_line_ID"="Metadata_Project.Alias")) %>%
    inner_join(pcs.terra, by = c("Metadata_line_ID"="Metadata_Project.Alias"))
  colnames(feature.gene.var) = sub("[-]", "_", colnames(feature.gene.var))
  print(paste("data prepared", nrow(feature.gene.var)))

  # burden test
  # regress_feature_allele(feature = feature, gene = gene, covars = covars, data = feature.gene.var, model = reg.model, level = level, inter.var = inter.var) %>% pull(p) %>% format(sc=T)
  res.perm = do.call("rbind", apply(genes.to.test, 1, function(x){ regress_feature_allele(feature = feature, gene = x, covars = covars, data = feature.gene.var, model = reg.model, level = level, inter.var = inter.var) }))
  res.perm = res.perm %>% mutate(feat = feature, cells = round, perm.seed = my.seed)
  
  # store results
  burden.res = rbind(burden.res, res.perm %>% arrange(p))
  print("burden results stored")
}
print("permutation done")

# save results
write.table(x = burden.res, file = paste0("/data/srlab/jatin/cmqtl/permut/", round, "/rv.featPermut.baselineSig.", suffix, ".newCovar.tab"), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
print(paste("burden test with permutation results written for", seed))

print("quitting...")
quit()
}
# combine permutation results
if(F){
###
cores=5
print("--> combining permutation results...")
all.res.files = list.files(paste0("/data/srlab/jatin/cmqtl/permut/",round), pattern=paste(level,phenotype,sum.or.present,gene.freq,"*.baseline",variant.source,"*newCovar.tab",sep="."), full.names=T)
print(length(all.res.files))
print(head(all.res.files,2))

block = split(1:length(all.res.files), cut_interval(1:length(all.res.files), n=cores))
files = all.res.files[block[[seed]]]
print(paste("working on",length(files)))

# summary from permutations
for(file in files){
  print(file)
  # 1.
  readr::read_delim(file, delim="\t") %>% filter(feat %in% features.to.test & var %in% genes.to.test$gene.2) %>% group_by(feat) %>%
    # summarise(quantile(p, c(0.01))) %>%
    # slice_max(order_by=c(-p), n=1) %>% arrange(p) , summarize(median(p))
    mutate(p.min = min(p)) %>% dplyr::select(feat,p.min,perm.seed) %>% unique() %>%
    readr::write_delim(sub(".*[/](.*).tab","extreme.values/permut.out/\\1.02.minp",file), delim="\t")
  
  # 2. trait by gene matrix filled with p vals
  # assoc.mat.rand = readr::read_delim(file, delim="\t") %>% filter(feat %in% features.to.test & var %in% genes.to.test$gene.2) %>%
  #   reshape2::dcast(feat~var, value.var="p")
  # prop.res = NULL
  # 
  # # number of associations below certain p value
  # for(p in c(10^-3,10^-4,10^-5,10^-6,10^-6.5,10^-7,10^-7.5,10^-8)){
  #   prop.res = rbind( prop.res, 
  #                     data.frame("p.th"=p, "assoc"=sum(assoc.mat.rand < p), "total"=nrow(assoc.mat.rand)*ncol(assoc.mat.rand)) )
  # }
  # prop.res %>% readr::write_delim(sub(".*[/](.*).tab","extreme.values/permut.out/\\1.02.summ",file), delim="\t")
  
}

print(paste("files saved for",seed,", quitting..."))
quit()
}
# select significant associations
if(F){
# load results for baseline
all.res.files = list.files("extreme.values/burden.out.11.2022/", pattern = paste("rv",round,level,phenotype,sum.or.present,gene.freq,"*.baseline",variant.source,"wellMean.tab",sep="."), full.names = T)
length(all.res.files)
head(all.res.files,2)
burden.res.bulk <- all.res.files %>% purrr::map_df(read.delim) %>% filter(data == "orig")

# effective number of traits and fdr
# ftrs = data.sampled.bulk.gaussian %>% dplyr::select(one_of(features.to.test))
# effective.traits = poolr::meff(R=cor(ftrs), method="gao")
sig.assoc.p = burden.res.bulk %>% mutate("p2" = p.adjust(p,"bon",n=nrow(genes.to.test)*length(features.to.test))) %>% 
  arrange(p) %>% filter(p2 < 0.05) %>% mutate("p.emp" = round(p2,3), "vjust" = 1.2)

# save associations for other analysis
potent.assoc = burden.res.bulk %>% filter(var %in% genes.to.test$gene.2) %>% arrange(p) %>% dplyr::select(feat,var,est,p,std.err)
assign(paste("potent.assoc",round,level,variant.source,sep="."), potent.assoc)

burden.res.bulk %>% filter(p < 10^-6) %>% dplyr::rename("gene"="var", "trait"="feat") %>% dplyr::select(trait,gene,est,p) %>% mutate(cells="isolate") %>%
  write.table(file = paste0("assoc_suggestive_evidence_",round,".tab"), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
  
burden.res.bulk %>% filter(p < 0.05 & var == "PRLR") %>% dplyr::rename("gene"="var", "trait"="feat") %>%
  dplyr::select(trait,gene,est,p) %>% mutate(cells="isolate") %>%
  write.table(file = paste0("supp_rare_assoc_PRLR_nominal_",round,".tab"), append = F, quote = F, sep = "\t", row.names = F, col.names = T)

# ----- 1. permutation-derived trait-specific threshold -----
# ------------------
# ------ temp ------
# load permuted associations (2 files only)
burden.res.perm.med <- list.files("tmp", pattern=paste(round,level,phenotype,sum.or.present,gene.freq,"33.*baseline",variant.source,sep="."), full.names = T) %>%
  purrr::map_df(read.delim) %>% filter(feat %in% features.to.test)
colnames(burden.res.perm.med) = c("feat","p.summ")
burden.res.perm.med %>% mutate("feat.cat" = sub(".*_(.*?)_.*", "\\1", feat)) %>% #filter(feat.cat == "Granularity") %>%
  
  # all summary files from permuted data    
  burden.res.perm.med <- list.files("extreme.values/permut.out", pattern=paste(round,level,phenotype,sum.or.present,gene.freq,"*baseline",variant.source,"*02.summ",sep="."), full.names = T) %>%
  purrr::map_df(read.delim)

# trait by gene matrix filled with p vals
assoc.mat = burden.res.bulk %>% filter(var %in% genes.to.test$gene.2) %>% reshape2::dcast(feat~var, value.var="p")

# trait by gene matrix filled with p vals
feat.more.p = which(table(burden.res.perm.med$feat)>10000) %>% names()
assoc.mat.rand = burden.res.perm.med %>% filter(var %in% genes.to.test$gene.2 & !feat %in% feat.more.p) %>% 
  reshape2::dcast(feat~var, value.var="p")

# estimate fdr
fdr.res = NULL
for(p in c(10^-3,10^-4,10^-5,10^-6,10^-6.5,10^-7,10^-7.5,10^-8)){
  print(-log10(p))
  assoc.len.fun <- function(x, pt=as.name(p)){
    return( length(which(x < pt)) )
  }
  
  # by theory
  # fdr.res = rbind( fdr.res, data.frame(p, "fdr" = (nrow(assoc.mat)*ncol(assoc.mat)*p)/sum(assoc.mat < p)) )
  
  # by association
  # a = sum(assoc.mat.rand < p)/(nrow(assoc.mat.rand)*ncol(assoc.mat.rand))
  a = (burden.res.perm.med %>% filter(p.th == p) %>% pull(assoc) %>% sum())/(nrow(burden.res.perm.med)/8)
  b = sum(assoc.mat < p)
  fdr.res = rbind( fdr.res, data.frame(p, "fdr" = a/b) )
  
  # by gene
  # assoc.len = assoc.mat %>% summarise_at(.vars=genes.to.test$gene.2, .funs=assoc.len.fun)
  # assoc.len.rand = assoc.mat.rand %>% summarise_at(.vars=genes.to.test$gene.2, .funs=assoc.len.fun)
  
  # # by feature
  # assoc.len = assoc.mat %>% tibble::column_to_rownames("feat") %>% t() %>% data.frame() %>%
  #   summarise_at(.vars=features.to.test, .funs=assoc.len.fun)
  # assoc.len.rand = assoc.mat.rand %>% tibble::column_to_rownames("feat") %>% t() %>% data.frame() %>%
  #   summarise_at(.vars=intersect(assoc.mat.rand$feat,features.to.test), .funs=assoc.len.fun)
  
  # store results
  # fdr.res = rbind( fdr.res, data.frame(p, "fdr" = length(which(assoc.len.rand[1,] > 0)) / length(which(assoc.len[1,] > 0))) )
}

fdr.res %>% ggplot(aes(-log10(p),fdr)) + geom_point() + geom_hline(yintercept=0.05, color="blue") + 
  ggtitle("permut:1000 rounds, all associations")

# assoc.mat[1:10,1:10]
# histogram(-log10(assoc.mat[,"WASF2"]), breaks=10)
# assoc.mat %>% readr::write_delim(path="~/Google Drive/cmQTL_Jatin_Samira/intermediate.trait.gene.p.mat", delim="\t")

# distribution of number of traits per gene associated with certain p vals
assoc.len.fun <- function(x, p=10^-5){
  return( length(which(x < p)) )
}
assoc.len.genes = assoc.mat %>% summarise_at(.vars=genes.to.test$gene.2, .funs=assoc.len.fun)
assoc.len.genes %>% t() %>% as.data.frame() %>% arrange(-V1) %>% filter(V1 > 0) %>%
  ggplot(aes(V1)) + geom_histogram(bins=50) + labs(x="number of linked traits", y="number of genes") + ggtitle("p < 10^-5")

burden.res.bulk %>% summarise_at(.vars=genes.to.test$gene.2, .funs=assoc.len)
# ------ temp ------
# ------------------


# load permutation results
all.res.files = list.files("extreme.values/permut.out", pattern=paste(round,level,phenotype,sum.or.present,gene.freq,"*baseline",variant.source,"*newCovar.02.minp",sep="."), full.names = T)
length(all.res.files)
head(all.res.files,2)
burden.res.perm.pmin <- all.res.files %>% purrr::map_df(read.delim) %>% filter(feat %in% features.to.test & var %in% genes.to.test$gene.2)


# ------------------
# ------ temp ------
# evaluate global threshold and class-specific
# feat.more.p = which(table(burden.res.perm.pmin$feat)>1000) %>% names()
# permut.x = burden.res.perm.pmin %>% filter(feat %in% feat.more.p) %>% group_by(feat) %>% sample_n(size=900)

# lowest p value per permutation round
burden.res.perm.all.pmin = burden.res.perm.pmin %>% mutate("feat.cat" = sub(".*_(.*?)_.*", "\\1", feat)) %>% #filter(feat.cat == "AreaShape") %>%
  #filter(!feat %in% c("Cytoplasm_Texture_DifferenceVariance_DNA_5_00","Cells_Intensity_LowerQuartileIntensity_Brightfield","Cytoplasm_Granularity_3_RNA")) %>%
  group_by(perm.seed) %>% slice_max(order_by=c(-p.min), n=1) %>% arrange(p.min)

# p value threshold
burden.res.perm.all.pmin %>% arrange(p.min) %>% head(floor(nrow(.)*0.05)) %>% tail(1) %>% pull(p) %>% format(sc=T)

# number of events for lowest p values
burden.res.perm.all.pmin$nvars = 0
for(i in 1:nrow(burden.res.perm.all.pmin)){
  gene = sub("_","-",as.character(burden.res.perm.all.pmin[i,"var"]))
  burden.res.perm.all.pmin[i,"nvars"] = gene.vars %>% dplyr::select(!!as.name(gene)) %>% filter(!!as.name(gene) != 0) %>% nrow()
}

# bar plot for lowest p values
burden.res.perm.all.pmin %>% filter(p < 10^-7) %>% #pull(feat) %>% table() %>% as.data.frame()
  ggplot(aes(feat)) + geom_bar(stat="count", width=0.5) + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) 

# histogram 
burden.res.perm.pmin %>% #filter(p < 10^-7) %>%
  ggplot(aes(-log10(p))) + geom_histogram(bins=50) + scale_x_continuous(breaks = scales::extended_breaks(n=5)) +
  geom_vline(xintercept=-log10(p.thresh)) + theme_bw()

# median or 1st percentile p value from permutation
burden.res.perm.med <- list.files("extreme.values/permut.out", pattern=paste(round,level,phenotype,sum.or.present,gene.freq,"*baseline",variant.source,"*02.01stp",sep="."), full.names = T) %>%
  purrr::map_df(read.delim) %>% filter(feat %in% features.to.test)
colnames(burden.res.perm.med) = c("feat","p.summ")
burden.res.perm.med %>% mutate("feat.cat" = sub(".*_(.*?)_.*", "\\1", feat)) %>% #filter(feat.cat == "Granularity") %>%
  ggplot(aes(feat, -log10(p.summ))) + geom_jitter(alpha=0.05, size=0.1) + geom_boxplot(outlier.shape=NA, alpha=0.5, color="#3182bd") + 
  theme(axis.text.x=element_blank()) + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  geom_hline(yintercept=3, color="red")
# ------ temp ------
# ------------------


# check observed p and est against background distribution
sig.assoc.p = NULL
for(feature in intersect(potent.assoc$feat, burden.res.perm.pmin$feat)){
  pthresh = burden.res.perm.pmin %>% filter(feat == !!feature) %>% arrange(p) %>% head(floor(nrow(.)*0.048)) %>% tail(1)
  sig.assoc.p = rbind( sig.assoc.p, potent.assoc %>% filter(feat == !!feature & p < pthresh$p) %>% mutate("pthresh" = pthresh$p) )
}

# p value derived from permutation
for(i in 1:nrow(sig.assoc.p)){
  sig.assoc.p[i,"p.emp"] = burden.res.perm.pmin %>% filter(feat == !!sig.assoc.p[i,"feat"] & p.min < !!sig.assoc.p[i,"p"]) %>% 
    nrow()/length(unique(burden.res.perm.pmin$perm.seed))
}

# save associations for other analysis
assign(paste("sig.assoc.p",round,level,sep="."), sig.assoc.p)

# format significant associations for plot
sig.assoc.p = sig.assoc.p %>% arrange(p.emp) %>% mutate("p.emp" = round(p.emp,3), "vjust" = 1.2)
sig.assoc.p[duplicated(sig.assoc.p$feat), "vjust"] = 2.5

# plot hist of p value against background
burden.res.perm.pmin %>% filter(feat %in% sig.assoc.p$feat) %>%
  ggplot(aes(-log10(p.min))) + geom_histogram(bins=50, aes(fill=..count..)) + theme_bw() + facet_wrap(~feat, scales="free") +
  geom_segment(data=sig.assoc.p %>% mutate("x"=-log10(p)), mapping=aes(x=x, y=0, xend=x, yend=5), arrow=arrow(length=unit(0.03,"npc"), ends="first"), color="red") +
  ggthemes::scale_fill_gradient_tableau(guide=F) + ggtitle("p from original data vs permuted data") +
  geom_text(data=sig.assoc.p, aes(feat, label=p.emp, vjust=vjust), x=Inf, y=Inf, hjust=1.05)

# plot hist of beta against background
burden.res.perm.emax %>% filter(feat %in% sig.assoc.p$feat) %>%
  ggplot(aes(est.max)) + geom_histogram(bins=50, aes(fill=..count..)) + theme_bw() + facet_wrap(~feat, scales="free") +
  geom_segment(data=sig.assoc.p %>% mutate("x"=abs(est)), mapping=aes(x=x, y=0, xend=x, yend=5), arrow=arrow(length=unit(0.03,"npc"), ends="first"), color="red") + 
  ggthemes::scale_fill_gradient_tableau(guide=F) + ggtitle("beta from original data vs permuted data")

# ----- 2. permutation-derived cell-bin-specific threshold -----
# min p per permutation
sig.thresh = burden.res.perm %>% filter(feat %in% potent.assoc$feat & var %in% genes.to.test$gene.2) %>% 
  group_by(perm.seed) %>% mutate(p.min = min(p)) %>% 
  dplyr::select(perm.seed,p.min) %>% unique() %>% arrange(p.min) %>% 
  head(ceiling(length(unique(burden.res.perm$perm.seed))*0.05)) %>% tail(1) %>% pull(p.min)

# significant associations
sig.assoc = burden.res.bulk %>% filter(p < sig.thresh)

# histogram of p value in bins
burden.res.perm %>% ggplot(aes(p)) + geom_bar() + scale_x_binned(n.breaks=1/0.05) + theme_bw()

# plot histrogram of p-values for permuted data
burden.res.perm %>% group_by(seed) %>% mutate(p.min = min(p)) %>% dplyr::select(seed,p.min) %>% unique() %>% arrange(p.min) %>% 
  ggplot(aes(-log10(p.min))) + geom_histogram(bins=50, aes(fill=..count..)) + theme_bw() +
  geom_segment(data=potent.assoc %>% mutate("x"=-log10(p)), mapping=aes(x=x, y=0, xend=x, yend=2), color="red") + 
  ggthemes::scale_fill_gradient_tableau(guide=F) + ggtitle("observed p-values against background distribution", subtitle=round) + 
  geom_vline(xintercept=-log10(sig.thresh))

# # ---- number of donors with variants per associations -----
# donors.variant = gene.vars %>% summarise_at(grep("Meta|Sample", colnames(.), value=T, invert=T), ~nnzero(.x)) %>%
#   colSums() %>% data.frame("n.vars" = .) %>% tibble::rownames_to_column("gene") %>%
#   inner_join(sigResults %>% mutate("gene" = sub("_","-",var)), by="gene")
# 
# donors.variant %>% ggplot(aes(p)) + geom_bar(aes(fill=factor(n.vars))) + scale_x_binned(n.breaks=1/0.05) + theme_bw() +
#   labs(title="permuted data")
# 
# donors.variant %>% ggplot(aes(factor(n.vars),-log10(p))) + geom_boxplot(width=0.3, outlier.size=0.25) + theme_bw() + labs(x="donors with variants") +
#   ggtitle(round, subtitle="all associations") + guides(color=F)
# # -----------

# plots
p1 = donors.variant %>% ggplot(aes(n.vars)) + geom_histogram(bins=30, aes(fill=..count..)) + theme_bw() + labs(x="donors with variants") +
  scale_x_continuous(breaks = seq(min(donors.variant$n.vars), max(donors.variant$n.vars), by = 4)) + 
  ggtitle(round) + guides(fill=F)
p2 = donors.variant %>% ggplot(aes(n.vars,-log10(p))) + geom_point(aes(color = -log10(p))) + theme_bw() + labs(x="donors with variants") +
  scale_x_continuous(breaks = seq(min(donors.variant$n.vars), max(donors.variant$n.vars), by = 4)) + 
  ggtitle(round) + guides(color=F)
p3 = burden.res %>% ggplot(aes(p)) + geom_histogram(bins=30, aes(fill=..count..)) + theme_bw() + 
  ggtitle(round, subtitle="all associations") + guides(color=F) + guides(fill=F)
gridExtra::grid.arrange(p3, p1, p2, nrow = 1)

###
donors.variant %>% ggplot(aes(factor(n.vars),-log10(p))) + geom_boxplot() + theme_bw() + labs(x="donors with variants") +
  scale_x_continuous(breaks = seq(min(donors.variant$n.vars), max(donors.variant$n.vars), by = 4)) + 
  ggtitle(round) + guides(color=F)
###
donors.variant %>% group_by(n.vars) %>% count() %>% dplyr::rename("total"="n") %>% 
  inner_join(donors.variant %>% group_by(n.vars) %>% filter(p < 0.05) %>% count(), by="n.vars") %>% 
  ggplot(aes(factor(n.vars), n/total)) + geom_col() + 
  labs(x="number of donors with variants", y="% of associations with p < 0.05")
###
donors.variant %>% ggplot(aes(factor(n.vars),est)) + geom_boxplot() + theme_bw() + labs(x="donors with variants") +
  scale_x_continuous(breaks = seq(min(donors.variant$n.vars), max(donors.variant$n.vars), by = 4)) + 
  ggtitle(round) + guides(color=F)

}
# plot associations, qq, microscopy images
if(F){
# ----------------- Plot associations ---------------
# take significant associations
# burden.res.bulk %>% arrange(p) %>% filter(p < 10^-6 & feat == "Cytoplasm_Granularity_3_RNA")
frame.to.plot = sig.assoc.p %>% dplyr::select(feat,var,est,p) %>%
  mutate(p = format(p, digits=1, sc=T), variable = sub("_", "-", sub("_inter","",var)), 
         assoc = ifelse(grepl("inter",var), "inter", "gene"))

# plot each feature individually
sapply(unique(frame.to.plot$feat), function(x) plot_association(feature = x, features.values = data.sampled.bulk.gaussian, 
              gene.values = gene.vars, frame.to.plot = frame.to.plot, interaction.feat = "Metadata_Sex", 
              model = reg.model, round=round, level=level, plot=T))

# plot all features together
plt.list = NULL
for (feat in unique(frame.to.plot$feat)){
  plt.list[[feat]] = plot_association(feature = feat, features.values = data.sampled.bulk.gaussian, gene.values = gene.vars, 
                                      frame.to.plot = frame.to.plot, interaction.feat = "Metadata_Disease_Category", 
                                      model = reg.model, round = round, level=level, plot=F)  
}
cowplot::plot_grid(plotlist=plt.list, ncol=4)

# --------------- Get ID of cell images -------------
# meta data of images
# images.meta = read.csv("~/Documents/cmQTL/images/cmqtl.3/sample_images_allPlates.csv") %>% dplyr::select(-one_of("URL"))
images.meta = read.csv("images/sample_images.csv") %>% dplyr::select(-one_of("URL"))
images.meta.sites = read.csv("bekar.cheeze/plate.metadata.csv") %>% unique() %>% 
  inner_join(images.meta, by = c("Metadata_Plate", "Metadata_Well"))

feature = "Cells_RadialDistribution_RadialCV_Mito_1of4"
gene = "PRLR"

images.meta.sites %>% filter(Metadata_Plate=="BR00106708" & grepl("r06c04f03p01",filename))

# get image id for baseline model
features.values %>% dplyr::select(c(Metadata_line_ID, !!as.name(feature))) %>%
  inner_join(gene.vars %>% dplyr::select(one_of("Metadata_Project.Alias",gene)), by = c("Metadata_line_ID"="Metadata_Project.Alias")) %>%
  filter(!!as.name(gene) != 0) %>% mutate(avg = median(!!as.name(feature)), dist.avg = abs(avg - !!as.name(feature))) %>% 
  arrange(dist.avg) %>% head(20) %>%
  inner_join(images.meta.sites, by = c("Metadata_line_ID")) %>% 
  mutate(image.name = paste(Metadata_Plate, sub("-.*","",filename), sep=":")) %>% 
  dplyr::select(Metadata_line_ID,Metadata_Plate,Metadata_Well,!!as.name(feature),!!as.name(gene),avg,dist.avg,image.name) %>% 
  unique() %>% dplyr::select(Metadata_line_ID,image.name) %>% unique()

# ----- temp -----
# compare seeded cells across cell lines with and without variants
features.values.cellCount = features.values %>% dplyr::select(c(Metadata_line_ID,!!as.name(feature)),Cells_Neighbors_NumberOfNeighbors_Adjacent) %>%
  inner_join(gene.vars %>% dplyr::select(one_of("Metadata_Project.Alias",gene)), by = c("Metadata_line_ID"="Metadata_Project.Alias")) %>%
  inner_join(cell.count.summ.line, by="Metadata_line_ID")

cellCount.withV = features.values.cellCount %>% filter(!!as.name(gene) != 0) %>% pull(Metadata_CellCount_avg)
cellCount.withoutV = features.values.cellCount %>% filter(!!as.name(gene) == 0) %>% pull(Metadata_CellCount_avg)
test = wilcox.test(cellCount.withV, cellCount.withoutV)

features.values.cellCount %>% ggplot(aes(as.factor(!!as.name(gene)),Metadata_CellCount_avg)) + geom_boxplot(width=0.5) + theme_bw() +
  annotate(geom="text", x=Inf, y=Inf, label = round(test$p.value,3), hjust=1.2, vjust=2) +
  geom_jitter(width=0.2)


gene="PRLR"
gene.cellCount = list.files(path="metadata", pattern="cell.cat.*.summ", full.names=T) %>% purrr::map_df(read.delim) %>%
  filter(Metadata_line_ID %in% donor.meta$Metadata_Project.Alias & cell.cat == "isolate") %>%
  inner_join(gene.vars %>% dplyr::select(one_of("Metadata_Project.Alias",gene)), by=c("Metadata_line_ID"="Metadata_Project.Alias"))

test = wilcox.test(gene.cellCount %>% filter(!!as.name(gene) != 0) %>% pull(n), 
                   gene.cellCount %>% filter(!!as.name(gene) == 0) %>% pull(n))

gene.cellCount %>%  ggplot(aes(as.factor(!!as.name(gene)),n)) + geom_boxplot(width=0.5) + theme_bw() + geom_jitter(width=0.2) +
  annotate(geom="text", x=Inf, y=Inf, label = paste("p =",round(test$p.value,3)), hjust=1.2, vjust=2) +
  ggtitle(gene)
# ----- temp -----

# --------- QQ plots for genes and features ---------
library("pryr")
burden.res = burden.res.bulk %>% arrange(p)
associated.feat = burden.res %>% filter(p.adj < 0.05) %>% pull(feat) %>% unique()
associated.genes = burden.res %>% filter(p.adj < 0.05) %>% pull(var) %>% unique()

# lambda and qq for associated genes
pdf(paste0("plots/qq.genes.",round,"_",reg.model,".pdf"), width = 2.5, height = 7)
par(mfrow=c(4,1))
for(g in associated.genes){
  pvals = burden.res %>% filter(var == !!g & data == "orig") %>% pull(p)
  a %<a-% qqman::qq(pvector = pvals, main = paste0(g, "\nlambda = ", round(GenABEL::estlambda(pvals, method="median")$estimate,2)))
  a; rm(a)
  print(g)
}
dev.off()

# qq for associated features
pdf(paste0("plots/qq.feats.",round,"_",reg.model,".pdf"), width = 10, height = 2.5)
par(mfrow=c(1,3))
for(f in associated.feat){
  # raw p_values
  ps.orig = burden.res %>% filter(feat == f) %>% pull(p)
  lambda = round(GenABEL::estlambda(ps.orig, method="median")$estimate,2)
  a %<a-% qqman::qq(pvector = ps.orig, 
                    main = paste(f,"\nin",ifelse(round=="intcol","Colony","Isolated"),"cells, lambda =",lambda))
  a; rm(a)
  
  # # permuted p_values
  # if(burden.res %>% filter(data == "perm") %>% nrow() > 0){
  #   ps.perm = burden.res %>% filter(feat == f & data == "perm") %>% pull(p)
  #   b %<a-% qqman::qq(pvector = ps.perm, main = paste0(f, "\n(p_values permute) lambda = ", round(GenABEL::estlambda(ps.perm)$estimate,2)))
  #   b; rm(b)
  # }
  print(f)
}
dev.off()


# -------- Plot summary of associated genes ---------
# beta coef vs p value of genes (5.5x4)
pdf(paste("plots/genes.feats",round,level,reg.model,"pdf",sep="."), width = 5.5, height = 4)
burden.res[sample(rownames(burden.res)), ] %>% mutate(label = ifelse(p.adj < 0.05 & est != 0, var, "")) %>% filter(p < 0.05 & data == "orig") %>%
  ggplot(aes(est, -log10(p))) + geom_point(stroke = 1.5, size = 1, shape = 1, aes(color = feat.cat)) + theme_bw() + 
  theme(panel.grid.major.x = element_blank()) + ggrepel::geom_text_repel(aes(label = label)) +
  ggtitle(paste(round, ", variants effect",variants), paste("features tested =", nrow(features.uncor), ", genes tested =",nrow(genes.to.test))) + 
  xlab("beta coeff. of gene") + scale_color_manual(values = feat.cat.cols) + 
  geom_hline(yintercept = -log10(0.05/(nrow(features.uncor)*nrow(genes.to.test))), linetype = "dashed", color = "blue") + 
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 13))
dev.off()

# frequency of associated genes per feature category
burden.res %>% filter(p.adj < 0.05 & data == "orig") %>% group_by(var) %>% mutate(freq = n()) %>% filter(freq > 0) %>%
  ggplot(aes(forcats::fct_infreq(var))) + geom_bar(stat="count", width = 0.5, aes(fill = feat.cat)) + theme(panel.grid = element_blank()) + 
  ggtitle(round) + coord_flip() + scale_fill_manual(values = feat.cat.cols) + ylab("number of associated features") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 13)) + xlab("")

}
# circos, manhattan plots, pathways enrichment
if(F){
library("biomaRt")
library("BioCircos")
  
# # all associations
# gwasResults = rbind(potent.assoc.intermediate.bulk %>% mutate("cell.group" = "intermediate"), 
#                     potent.assoc.isolate.bulk %>% mutate("cell.group" = "isolate"), 
#                     potent.assoc.colony.bulk %>% mutate("cell.group" = "colony"))
# sigResults = rbind(sig.assoc.p.intermediate.bulk %>% mutate("cell.group" = "intermediate"), 
#                    sig.assoc.p.isolate.bulk %>% mutate("cell.group" = "isolate"), 
#                    sig.assoc.p.colony.bulk %>% mutate("cell.group" = "colony"))
  
gwasResults = rbind(potent.assoc.isolate.bulk.rareCM %>% mutate("cell.group" = "isolate"), 
                    potent.assoc.intcol.bulk.rareCM %>% mutate("cell.group" = "intcol"))
sigResults = rbind(potent.assoc.isolate.bulk.rareCM %>% filter( p < 0.05/(nrow(genes.to.test)*length(features.to.test)) ) %>% mutate("cell.group" = "isolate"),
                   potent.assoc.intcol.bulk.rareCM %>% filter( p < 0.05/(nrow(genes.to.test)*length(features.to.test)) ) %>% mutate("cell.group" = "intcol"))

# venn diagram of overlap berween different cell groups
VennDiagram::venn.diagram( x=list(sigResults %>% filter(cell.group == "intcol") %>% pull(var), sigResults %>% filter(cell.group != "intcol") %>% pull(var)), 
             category.names=c("colony","isolate"), filename="plots/venn.png", resolution=400, 
             height=2, width=2, units="in", col=c("#440154ff", '#21908dff'), 
             cat.cex = 0.8, cat.dist = c(0,0), cat.default.pos = "text", cat.just=list(c(0.5,2),c(0.5,2)),
             fill=c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)), 
             main="associated genes" )

# load results
for(round in c("isolate","intermediate")){
  all.res.files = list.files("extreme.values/burden.out", pattern = paste("rv",round,"bulk",phenotype,sum.or.present,"*.baseline.rareCM.tab",sep="."), full.names = T)
  print(paste(length(all.res.files), "files for", round))
  assign(x = paste0(round,".1.res"), value = all.res.files %>% purrr::map_df(read.delim))
}

for(round in c("anycells","isolate","colony","intermediate")){
  #all.res.files = list.files("extreme.values/burden.out", pattern = paste("rv",round,"well",phenotype,sum.or.present,"*.baseline.rareCM.tab",sep="."), full.names = T)
  all.res.files = list.files("extreme.values/burden.out", pattern = paste("rv",round,"bulk",phenotype,sum.or.present,"*.baseline.rareCM.db.time.tab",sep="."), full.names = T)
  print(paste(length(all.res.files), "files for", round))
  assign(x = paste0(round,".2.res"), value = all.res.files %>% purrr::map_df(read.delim))
}

# get all gene symbols from human genome
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
all.genes = getBM(attributes = c("hgnc_symbol","chromosome_name"), mart=mart, uniqueRows = T) %>% filter(chromosome_name %in% 1:22) %>% 
  dplyr::pull(hgnc_symbol)

# ------------------ manhattan plot ----------------
###
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
attributes = c("ensembl_gene_id","start_position","end_position","hgnc_symbol","chromosome_name")

# get positon of genes
genes.pos = getBM(attributes=attributes, filters=c("hgnc_symbol"), values=list(hgnc_symbol=unique(gwasResults$var)), mart=mart, uniqueRows=T)

# merge associations and gene positions
gwasResults.pos = genes.pos %>% filter(chromosome_name %in% 1:22) %>% mutate("chrom" = as.numeric(chromosome_name)) %>% 
  inner_join(gwasResults, by = c("hgnc_symbol"="var")) %>% 
  dplyr::select(-one_of("ensembl_gene_id","data","variants","std.err","p.adj","chromosome_name")) %>%
  mutate(feat.cat = sub(".*_(.*?)_.*","\\1", feat), "pos" = (start_position+end_position)/2)

###
gwasResults.pos.format = gwasResults.pos %>% group_by(chrom) %>% summarise(chr_len=max(pos)) %>% 
  arrange(chrom) %>% mutate(tot=cumsum(chr_len)-chr_len) %>%
  left_join(gwasResults.pos, by="chrom") %>% arrange(chrom,pos) %>% mutate(pos_cum=pos+tot) %>% 
  left_join(sigResults %>% mutate(sig="y") %>% dplyr::select(feat,var,sig,cell.group), by=c("feat","hgnc_symbol"="var","cell.group")) %>%
  mutate("sig" = ifelse(p < 10^-6, feat.cat, ifelse((chrom/2)%%1==0, "grey", "black"))) %>%
  mutate("alpha" = ifelse(sig %in% c("grey","black"), 0.75, 1)) %>%
  mutate("cell.group" = ifelse(!is.na(cell.group), cell.group, "non_sig")) %>%
  mutate("size" = ifelse(p < 10^-6, 1, 0.5)) %>%
  mutate(cell.group=plyr::mapvalues(cell.group, c("isolate","intcol"), c("Isolated cells","Cells in colony")))
  
# axis ticks position
axisdf = gwasResults.pos.format %>% group_by(chrom) %>% summarize( center=(min(pos_cum)+max(pos_cum))/2 )

# plot manhattan
label.data = gwasResults.pos.format %>% filter(p < 0.05/(nrow(genes.to.test)*length(features.to.test))) %>%
  group_by(pos_cum) %>% slice_max(order_by=-log10(p), n=1)

# plot
man.p = gwasResults.pos.format %>% arrange(chrom,pos_cum) %>%
  ggplot(aes(pos_cum,-log10(p))) + geom_point(aes(alpha=alpha, color=sig, shape=cell.group, size=factor(size))) +
  labs(x="Chromosome", y="-Log10(p)") + guides(alpha=F, size=F, shape=guide_legend(title="Cell type")) + 
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center, expand = c(0,0)) + scale_y_continuous(expand = c(0.025,0)) +
  theme_bw() + scale_size_manual(values=c(0.5,1)) + theme(panel.grid.major=element_blank()) +
  scale_color_manual(values=c(feat.cat.cols,"grey"="#d9d9d9","black"="#525252"), name="Trait category", 
                     breaks=names(feat.cat.cols), labels=names(feat.cat.cols)) +
  geom_hline(yintercept=6, size=0.35, color="grey", linetype="dashed") + 
  geom_hline(yintercept=-log10(0.05/(nrow(genes.to.test)*length(features.to.test))), size=0.25, linetype="dashed") +
  theme(legend.key.size=unit(0.5,"cm"), legend.text=element_text(size=7), legend.title=element_text(size=8)) +
  ggrepel::geom_text_repel(data=label.data, aes(pos_cum,-log10(p), label=hgnc_symbol), size=3)

# save plot
#png(filename="plots/rare.associations.manhattan.png", width=8, height=2.75, units='in', res=100)
pdf(file="plots/rare.associations.manhattan.pdf", width=8, height=2.75)
print(man.p)
dev.off()

###
# sigResults %>% group_by(feat,var) %>% mutate(n = n()) %>% arrange(-n)

# ------------------- circos plot ------------------
# functin to make track for Biocircos
make_track <- function(rv.res, chr.lengths, track.y){
  
  ###
  mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  attributes = c("ensembl_gene_id","start_position","end_position","hgnc_symbol","chromosome_name")
  ###
  genes.pos = getBM(attributes=attributes, filters=c("hgnc_symbol"), values=list(hgnc_symbol=unique(rv.res$var)), mart=mart, uniqueRows = T)
  
  ###
  gene.coords = genes.pos %>% filter(chromosome_name %in% 1:22) %>% mutate("chrom" = paste0("chr",chromosome_name)) %>% 
    inner_join(chr.lengths, by = "chrom") %>% inner_join(rv.res, by = c("hgnc_symbol"="var")) %>% mutate(feat.cat = sub(".*_(.*?)_.*","\\1", feat)) %>%
    dplyr::select(one_of("feat.cat","hgnc_symbol","start_position","end_position","chrom","length","p","p.adj"))
  
  # add one track for each chromosome
  tracks = BioCircosTracklist()
  for(chr in unique(gene.coords$chrom)){
    # define histogram/bars to be displayed
    barValues = gene.coords %>% filter(chrom == chr) %>% pull("p") %>% log10()*-1
    barColor = gene.coords %>% filter(chrom == chr) %>% mutate("color" = ifelse(p.adj < 0.1, feat.cat.cols[feat.cat], "#bababa")) %>% pull("color")
    
    # start and end positions
    st = gene.coords %>% filter(chrom == chr) %>% pull(start_position)/10^6 - 0.1
    en = gene.coords %>% filter(chrom == chr) %>% pull(end_position)/10^6
    
    # add a bar track
    # tracks_bar = tracks_bar + BioCircosBarTrack(trackname = paste0("bars",chr), chromosome = chr, starts = st, ends = en, values = barValues, 
    #                                             color = barColor, range = c(0,12), maxRadius = 1.4, minRadius = 1.2)
    # add a snp track
    tracks = tracks + BioCircosSNPTrack(trackname = paste0("snps_",chr), chromosomes = chr, positions = en, values = barValues, size = 2,
                                        minRadius=track.y, maxRadius=track.y+0.25, colors = barColor, range = c(-log10(min(gene.coords$p)),0) )
  }
  return(tracks)
}

# define a custom genome with chr names and their length
chr.lengths = GenomicFeatures::getChromInfoFromUCSC("hg38") %>% head(22)
lengthChr = setNames(chr.lengths$length/10^6, chr.lengths$chrom)

# make tracks
tr1 = make_track(rv.res = colony.1.res %>% filter(p < 10^-3 & data == "orig"), chr.lengths = chr.lengths, track.y = 0.725)
tr2 = make_track(rv.res = colony.2.res %>% filter(p < 10^-3 & data == "orig"), chr.lengths = chr.lengths, track.y = 0.725-0.225)
# tr2 = make_track(rv.res = intermediate.1.res %>% filter(p < 10^-3 & data == "orig"), chr.lengths = chr.lengths, track.y = 0.725-0.225)
# tr3 = make_track(rv.res = colony.1.res %>% filter(p < 10^-3 & data == "orig"), chr.lengths = chr.lengths, track.y = 0.725-(0.225*2))
# tr4 = make_track(rv.res = isolate.1.res %>% filter(p < 10^-3 & data == "orig"), chr.lengths = chr.lengths, track.y = 0.725-(0.225*3))


# add background to tracks
tr1 = tr1 + BioCircosBackgroundTrack("any cells", fillColors = "#FFFFFF", minRadius = 0.80, maxRadius = 0.99, borderColors = "#bababa")
tr2 = tr2 + BioCircosBackgroundTrack("isolate", fillColors = "#FFFFFF", minRadius = 0.55, maxRadius = 0.79, borderColors = "#bababa")
tr3 = tr3 + BioCircosBackgroundTrack("colony", fillColors = "#FFFFFF", minRadius = 0.35, maxRadius = 0.54, borderColors = "#bababa")
tr4 = tr4 + BioCircosBackgroundTrack("colony", fillColors = "#FFFFFF", minRadius = 0.05, maxRadius = 0.34, borderColors = "#bababa")

# plot circos
tracks = tr1 + tr2 + tr3 + tr4
BioCircos(tracklist = tracks, genome = as.list(lengthChr), genomeTicksDisplay = F, genomeLabelDy = 10, chrPad = 0.01, 
          genomeFillColor = "Spectral", zoom = F, displayGenomeBorder = F, SNPMouseOverDisplay = F, genomeLabelOrientation = 90)


# --------------- compare 2 results ----------------
# 3x7
# function to compare and plot two results
compare_two_res <- function(res1, name1, res2, name2){
  
  # merge all results from multiple results
  res.both = res1 %>% dplyr::select(feat,var,p,est,std.err) %>%
    inner_join(res2 %>% dplyr::select(feat,var,p,est,std.err), by=c("feat","var"), suffix=c(".1",".2")) %>%
    mutate("r1" = est.1/std.err.1, "r2" = est.2/std.err.2) %>% 
    mutate("z1" = (est.1 - mean(est.1))/sd(est.1), "z2" = (est.2 - mean(est.2))/sd(est.2))
  print(paste("number of common associations",nrow(res.both)))
  
  # correlation between p and estimates
  r.p = round(as.numeric(cor.test(res.both$p.1, res.both$p.2)$estimate),2)
  r.b = round(as.numeric(cor.test(res.both$r1, res.both$r2)$estimate),2)
  r.z = round(as.numeric(cor.test(res.both$z1, res.both$z2)$estimate),2)
  
  ###
  title = paste(name1, "vs", name2)
  # frame.plot = res.both %>% filter(p.1 < 0.00001 | p.2 < 0.00001) #%>% mutate("sig" = ifelse(p.adj.1 < 0.1 | p.adj.2 < 0.1, "yes", "no"))
  frame.plot = res.both %>% mutate("sig" = ifelse(p.1 < 0.05/(nrow(genes.to.test)*length(features.to.test)) | p.2 < 0.05/(nrow(genes.to.test)*length(features.to.test)), "yes", "no"))
  
  # 1. p values
  p1 = ggplot(frame.plot, aes(-log10(p.1),-log10(p.2))) + geom_point(alpha = 0.25, aes(color=sig)) +
    theme_bw() + ggthemes::scale_color_tableau() +
    xlim( min(-log10(c(frame.plot$p.1, frame.plot$p.2))), max(-log10(c(frame.plot$p.1, frame.plot$p.2))) ) +
    ylim( min(-log10(c(frame.plot$p.1, frame.plot$p.2))), max(-log10(c(frame.plot$p.1, frame.plot$p.2))) )
    labs(x=paste("-log10(p) for",name1),y=paste("-log10(p) for",name2)) +
    annotate(geom="text", x = Inf, y = -Inf, label = paste("r =",r.p), hjust = 1.2, vjust = -1.25, size = 5, parse = F)
  
  # 2. beta by std error
  p2 = ggplot(frame.plot, aes(r1,r2)) + geom_point(size = 0.5, alpha = 0.25, aes(color=sig)) +
    theme_bw() + ggthemes::scale_color_tableau() +
    labs(x=paste("beta for",name1),y=paste("beta for",name2)) +
    annotate(geom="text", x = Inf, y = -Inf, label = paste("r =",r.b), hjust = 1.2, vjust = -1.25, size = 5, parse = F)
  
  # 3. z(beta by std error)
  p3 = ggplot(frame.plot, aes(z1,z2)) + geom_point(alpha = 0.25, aes(color=sig)) +
    theme_bw() + ggthemes::scale_color_tableau() + 
    xlim( min(c(frame.plot$z1, frame.plot$z2)), max(c(frame.plot$z1, frame.plot$z2)) ) +
    ylim( min(c(frame.plot$z1, frame.plot$z2)), max(c(frame.plot$z1, frame.plot$z2)) ) +
    labs(x=paste("z(beta) for",name1),y=paste("z(beta) for",name2)) + 
    annotate(geom="text", x = Inf, y = -Inf, label = paste("r =",r.z), hjust = 1.2, vjust = -1.25, size = 5, parse = F)
  
  # save plot
  png(filename="plots/compare.two.res.png", height=5, width=3.5, units="in", res=100)
  gridExtra::grid.arrange(p1,p3, ncol=1, top=title)
  dev.off()
  
}

# inputs and call to function
res1 = potent.assoc.intcol.bulk.rareCM #%>% filter(var %in% "WASF2" & p < 0.05)
name1 = "colony_all_rare"
res2 = potent.assoc.intcol.bulk.rareCM.Gnomad #%>% filter(var %in% "WASF2" & p < 0.05)
name2 = "colony_GnomadRare"

# 3x11
compare_two_res(res1=res1, name1=name1, res2=res2, name2=name2)

# --------------- pathway enrichment ---------------
# function to calculate pathway enrichment
calculate_enrichment <- function(de.genes, pathways, min.genes.inPath, all.genes.names){
  # iterate over each pathway, and check if DE genes show any enrichment
  permt.res = NULL
  ftest.res = NULL
  pathways.tested = 0
  
  for(path.name in names(pathways)){
    
    genes.inPath = intersect(de.genes, pathways[[path.name]])
    
    if(length(genes.inPath) >= min.genes.inPath){
      
      print(paste("checking pathway", path.name))
      pathways.tested = pathways.tested + 1
      
      # ----- Permut test -----
      rand.len = NULL
      for(i in 1:100000){
        rand.len = c(rand.len, length(intersect(sample(x = all.genes.names, size = length(de.genes)), pathways[[path.name]])))
      }
      p.val = (length(which(rand.len >= length(genes.inPath)))/100000) #* length(pathways)
      
      if(p.val < 0.05) permt.res = rbind(permt.res, data.frame("pathway"=path.name, "in.path.genes"=length(genes.inPath), "pval.perm"=round(p.val,7)))
      rm(p.val)
      
      # ----- Fisher test -----
      a = length(genes.inPath)
      b = length(de.genes) - a
      c = length(intersect(setdiff(all.genes.names, de.genes), pathways[[path.name]]))
      d = length(setdiff(all.genes.names, de.genes)) - c
      cont.mat = matrix(c(a,b,c,d), nrow = 2)
      colnames(cont.mat) = c("de.genes", "all.genes")
      rownames(cont.mat) = c("in.path", "not.in.path")
      ftest = fisher.test(cont.mat)
      p.val = ftest$p.value #* length(pathways)
      
      if(p.val < 0.05) ftest.res = rbind(ftest.res, data.frame("pathway"=path.name, "in.path.genes"=length(genes.inPath), "or"=round(ftest$estimate,2), "pval.ftest"=round(p.val,7)))
      rm(p.val)
      
    } # if condition on genes.inPath
    
  } # path.name loop
  
  rownames(ftest.res) = NULL
  return(list("permut" = permt.res, "ftest" = ftest.res, "pathways.tested" = pathways.tested))
  
}

# genes to test for enrichment
associated.genes = rbind(burden.res.bulk.iso %>% filter(p < 10^-6), burden.res.bulk.col %>% filter(p < 10^-6)) %>% 
  pull(var) %>% unique()

# get all gene symbols from human genome
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
all.genes = getBM(attributes = c("hgnc_symbol","chromosome_name"), mart=mart, uniqueRows = T) %>% 
  filter(chromosome_name %in% 1:22) %>% dplyr::pull(hgnc_symbol)

# load pathway file
pathway.files = list.files(path = "gsea/pathways", pattern = "c[5].*.v7.2.symbols.gmt", full.names = T)

# calculate enrichment
path.res = NULL
for(path.file in pathway.files){
  print(paste("on pathway", path.file))
  if(grepl("c2",path.file)) min.genes=3 else min.genes=3
  # enrichment in pathways
  enrich = calculate_enrichment(de.genes = associated.genes, pathways = fgsea::gmtPathways(path.file), 
                                min.genes.inPath = min.genes, all.genes.names = setdiff(all.genes,associated.genes))
  # filter enrichment res
  if(enrich$pathways.tested != 0 & !is.null(enrich$permut)){
    path.res = rbind( path.res, inner_join(x = enrich$permut, y = enrich$ftest, by = c("pathway","in.path.genes")) %>% 
                       mutate("p.adj" = p.adjust(pval.perm,"bon",enrich$pathways.tested), 
                              "path.cat" = sub(".*/(.*).v[0-9].*.symbols.gmt","\\1",path.file)) )
  }
}

# pre-process results and plot
path.res.to.plot = path.res %>% filter(!grepl("UNKNOWN", pathway)) %>% 
  mutate( pathway = Hmisc::capitalize(tolower(pathway)), "p.adj.2" = p.adjust(p.adj,"bon",length(unique(path.res$path.cat))) )

# dots (preferred)
p = path.res.to.plot %>% filter(p.adj.2 < 0.05 & !grepl("regulation_of_biosynthetic_process", pathway)) %>% 
  arrange(path.cat,-or) %>% mutate("pathway"=factor(pathway, levels=unique(pathway))) %>%
  ggplot(aes(pathway,or)) + geom_point(aes(size = factor(in.path.genes), color = path.cat)) + 
  coord_flip() + theme_bw() + labs(y="enrichment (or)", x="pathways", size="genes in pathway", color="pathway category") +
  ggthemes::scale_color_tableau() + scale_size_manual(values=c(2,3,4)) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size=8))
pdf(file="plots/pathway.enrich.pdf", width=7, height=2.5)
print(p)
dev.off()

# bars
path.res.to.plot %>% ggplot(aes(pathway, or)) + geom_col(width = 0.35, aes(fill = path.cat)) + 
  labs(y = "enrichment (odds ratio)", x = "pathways") + ggthemes::scale_fill_tableau() + coord_flip() + theme_bw() + 
  theme(panel.grid.minor = element_blank())

# --------------- disease enrichment ---------------
# ebi gwas catalog
gwas.cat = readr::read_delim("gsea/gwas_catalog_v1.0-associations_e100_r2020-08-13.tsv", delim = "\t")

# function to get gene-associated genes
get_associated_diseases <- function(gene.set, gwas.cat){
  associated.diseases = NULL
  for(gene in gene.set){
    associated.diseases = c(associated.diseases, gwas.cat %>% dplyr::select(one_of("DISEASE/TRAIT","REPORTED GENE(S)","MAPPED_GENE","OR or BETA","P-VALUE")) %>% 
                              filter((grepl(gene,`REPORTED GENE(S)`) | grepl(paste0("^",gene,"$"),MAPPED_GENE)) & !is.na(`OR or BETA`) & !is.na(`REPORTED GENE(S)`) & `P-VALUE` < 10^-8) %>%
                              pull("DISEASE/TRAIT") %>% unique())
  }
  return(associated.diseases %>% unique() %>% length())
}

# genes to test for enrichment
associated.genes = sigResults %>% pull(var) %>% unique()
n.assoc.diseases = get_associated_diseases(gene.set=associated.genes, gwas.cat=gwas.cat)

# random set of genes
all.genes = setdiff(genes.to.test$gene.2, associated.genes)
n.assoc.diseases.rand = NULL
for(i in 1:1000){
  print(i)
  n.assoc.diseases.rand = c(n.assoc.diseases.rand, 
                            get_associated_diseases(gene.set=sample(all.genes,length(associated.genes)), gwas.cat=gwas.cat))
}

# p value of enrichment
length(which(n.assoc.diseases.rand > n.assoc.diseases))/length(n.assoc.diseases.rand)


disease.table = table(associated.diseases) %>% as.data.frame(stringsAsFactors = F) %>% dplyr::rename("disease" = "associated.diseases")
disease.table = disease.table %>% mutate("cat" = ifelse(grepl("cup-disc|macular|optic|myopia|corneal|Cortical",disease,ignore.case=T),"optic", disease))
disease.table = disease.table %>% mutate("cat" = ifelse(grepl("count|of white|red cells",disease,ignore.case=T),"cell count", cat))
disease.table = disease.table %>% mutate("cat" = ifelse(grepl("Asthma",disease,ignore.case=T),"asthma", cat))
disease.table = disease.table %>% mutate("cat" = ifelse(grepl("Anthro|body|cholesterol|omega-6|Peripheral artery",disease,ignore.case=T),"cardio", cat))
disease.table = disease.table %>% mutate("cat" = ifelse(grepl("Parkinson",disease,ignore.case=T),"parkinson", cat))
disease.table = disease.table %>% mutate("cat" = ifelse(grepl("Parkinson|Insomnia|Morning|Cognitive|Schizo|Bipolar|hurt|neuroticism|mood",disease,ignore.case=T),"pysc", cat))
disease.table = disease.table %>% mutate("cat" = ifelse(grepl("level",disease,ignore.case=T),"protein level", cat))
disease.table = disease.table %>% mutate("cat" = ifelse(grepl("Headache|Migraine",disease,ignore.case=T),"migraine", cat))
disease.table = disease.table %>% mutate("cat" = ifelse(grepl("bowel|crohn",disease,ignore.case=T),"ibd", cat))

# plot
disease.table %>% group_by(cat) %>% mutate("count" = n()) %>% arrange(count) %>% filter(count >= 2) %>% dplyr::select(cat,count) %>% unique() %>%
  ggplot(aes(cat,count)) + geom_col(width = 0.35, aes(fill = cat)) + theme_bw() + ggthemes::scale_fill_tableau(guide=F) + 
  labs(x = "disease category", y = "number of linked genes")
  
}
# baaki cheeza
if(F){
# # ----- temp -----
# cut.res = NULL
# for(cut in seq(0.5,1,0.1)){
#   cv = findCorrelation(x = cor(ftrs), cutoff = cut, names = T, exact = F)
#   cut.res = rbind(cut.res, data.frame(cut, "nfeats" = length(setdiff(colnames(ftrs), cv))))
# }
# cut.res %>% ggplot(aes(factor(cut), nfeats)) + geom_col(width=0.35) + theme_bw() +
#   labs(x="correlation cutoff", y = "number of features")

ftrs.cor = data.sampled.bulk.gaussian %>% dplyr::select(one_of( replicate_correlation_values %>% filter(median >= 0.5) %>%
                                                                  pull(variable) )) %>% cor(.)
feature = "Cells_Intensity_IntegratedIntensity_ER"
ftrs.cor %>% as.data.frame() %>% filter(!!as.name(feature) > 0.9) %>% dplyr::select(one_of(feature)) %>%
  arrange(!!as.name(feature))

ftrs.cor["Cells_AreaShape_Zernike_0_0",] %>% data.frame("r"=.) %>% tibble::rownames_to_column("feat") %>% 
  filter(grepl("Cells_AreaShape",feat)) %>% arrange(r)
# ----- temp -----
  
# agreement between typical regression and on residues
res = regress_with_residues(feature = feature, genes.to.test = genes.to.test, covars = covars, data = feature.gene.var, inter.var = inter.var)

frame = res %>% mutate(feat = feature, cells = round) %>%
  inner_join(burden.res.bulk %>% filter(feat == feature), by=c("feat","var"), suffix=c("_resid","_orig")) %>%
  mutate("z_orig" = (est_orig - mean(est_orig))/sd(est_orig), "z_resid" = (est_resid - mean(est_resid))/sd(est_resid))

# correlation between p and estimates
r.p = round(as.numeric(cor.test(frame$p_orig, frame$p_resid)$estimate),2)
r.b = round(as.numeric(cor.test(frame$est_orig, frame$est_resid)$estimate),2)
r.z = round(as.numeric(cor.test(frame$z_orig, frame$z_resid)$estimate),2)

pe = frame %>% ggplot(aes(est_orig,est_resid)) + geom_point(alpha = 0.8, color = "grey") + ggtitle(feature) + 
  geom_abline(slope=1, intercept=0) + theme_bw() + 
  annotate(geom="text", x = Inf, y = -Inf, label = paste("r =",r.b), hjust = 1.2, vjust = -1.25, size = 5, parse = F)

pz = frame %>% ggplot(aes(z_orig,z_resid)) + geom_point(alpha = 0.8, color = "grey") + ggtitle(feature) + 
  geom_abline(slope=1, intercept=0) + theme_bw() + 
  annotate(geom="text", x = Inf, y = -Inf, label = paste("r =",r.z), hjust = 1.2, vjust = -1.25, size = 5, parse = F)

pp = frame %>% ggplot(aes(-log10(p_orig),-log10(p_resid))) + geom_point(alpha = 0.8, color = "grey") + ggtitle(feature) + 
  geom_abline(slope=1, intercept=0) + theme_bw() + 
  annotate(geom="text", x = Inf, y = -Inf, label = paste("r =",r.p), hjust = 1.2, vjust = -1.25, size = 5, parse = F)

gridExtra::grid.arrange(pp, pe, pz, nrow=1)


# summary of variant's impact
donor.gt %>% select(c(CHROM, POS, IMPACT)) %>% unique() %>% count(IMPACT) %>%
  ggplot(aes(IMPACT, n)) + geom_bar(width = 0.35, stat = "identity", aes(fill = IMPACT)) + theme_bw() +
  geom_text(aes(label=n), vjust=-0.25, size = 5) + xlab("variant impact") + ylab("count") + ggtitle("included variants") +
  scale_fill_manual(values = my.cols[-c(1:2)], guide = F) + scale_y_continuous(expand = expansion(mult = c(0,.1))) +
  theme(axis.text.x = element_text(size = 13), axis.title = element_text(size = 13), panel.grid.major = element_blank())

# summary of variant's effect, 7x6
donor.gt %>% select(c(CHROM, POS, EFFECT)) %>% unique() %>% group_by(EFFECT) %>% mutate(n = n()) %>% filter(n > 100) %>%
  ggplot(aes(forcats::fct_infreq(EFFECT))) + geom_bar(width = 0.35, stat = "count", aes(fill = EFFECT)) + theme_bw() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.25) + scale_y_continuous(expand = expansion(mult = c(0,.1))) +
  xlab("variant category") + ggtitle("high and moderate impact variants") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12), axis.title = element_text(size = 13), panel.grid.major = element_blank()) +
  scale_fill_manual(values = my.cols, guide = F)


# ----------------------------------------
# cross validation across all types of cells
reg.model = "baseline"
burden.bulk.res.all = NULL
burden.well.res.all = NULL
for(round in c("anycells","isolate")){
  print(round)
  # bulk level
  all.bulk.files = list.files("extreme.values/burden.out", pattern = paste("rv",round,"bulk",phenotype,sum.or.present,"*",reg.model,"rareGnomad.rareCM","tab",sep="."), full.names = T)
  burden.bulk.res.all = rbind(burden.bulk.res.all, all.bulk.files %>% purrr::map_df(read.delim))
  print("bulk done")
  # well level
  all.well.files = list.files("extreme.values/burden.out", pattern = paste("rv",round,"well",phenotype,sum.or.present,"*",reg.model,"rareGnomad.rareCM","tab",sep="."), full.names = T)
  burden.well.res.all = rbind(burden.well.res.all, all.well.files %>% purrr::map_df(read.delim))
  print("well done")
}

# --------------- 1. correlation among p values between cell types at SAME level ---------------
# tmp = burden.bulk.res.all %>% filter(p.adj < 0.05 & data == "orig") %>% mutate(variable = sub("_.*","",var)) %>% select(feat,variable) %>% unique() %>%
#   {filter(burden.bulk.res.all %>% mutate(variable = sub("_.*","",var)), variable %in% .$variable & feat %in% .$feat)} %>% 
#   select(variable,feat,est,std.err,p,cells) %>% mutate(beta.se = est/std.err)
# # frame = reshape2::dcast(data = tmp, variable + feat ~ cells, value.var = "p") #%>% select(-one_of("isolate"))
# 
# # pearson correlation test
# ctest = cor.test(frame[complete.cases(frame), "anycells"], frame[complete.cases(frame), "isolate"])
# 
# # plot correlation
# frame[complete.cases(frame), ] %>% 
#   ggplot(aes(anycells, isolate)) + geom_point(size = 2, alpha = 0.85) + #geom_abline(slope = ctest$estimate) +
#   annotate(geom = "text", x = Inf, y = Inf, label = paste("r =", round(ctest$estimate,2)), hjust = 1.1, vjust = 2, size = 4, parse = F) +
#   annotate(geom = "text", x = Inf, y = Inf, label = paste("P =", round(ctest$p.value,2)), hjust = 1.1, vjust = 3.5, size = 4, parse = F) + 
#   theme_bw() + ggthemes::scale_color_tableau() + ggtitle("beta/std. error")
#   
# # plot summary of associated genes
# burden.all.res %>% filter(p.adj < 0.05 & data == "orig") %>% group_by(variable) %>% mutate(freq = n()) %>%
#   ggplot(aes(forcats::fct_infreq(variable))) + geom_bar(stat="count", width = 0.5, aes(fill = cells)) + 
#   theme_bw() + theme(panel.grid.major.x = element_blank()) + 
#   coord_flip() + scale_fill_manual(values = my.cols) + ylab("number of associated features") +
#   theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13)) + xlab("variant burden in genes")

# --------------- 2. correlation among p values between cell types at DIFFERENT level ---------------
# get those things at well level which are significant at donor level
tmp1 = burden.bulk.res.all %>% filter(p.adj < 0.05 & data == "orig") %>% mutate(level = "bulk") %>% select(-one_of("variants","data")) %>%
  inner_join(burden.well.res.all %>% mutate(level = "well") %>% select(-one_of("variants","data")), by = c("var","feat","cells"), suffix = c(".bulk",".well"))
# other way around
tmp2 = burden.well.res.all %>% filter(p.adj < 0.05 & data == "orig") %>% mutate(level = "well") %>% select(-one_of("variants","data")) %>%
  inner_join(burden.bulk.res.all %>% mutate(level = "bulk") %>% select(-one_of("variants","data")), by = c("var","feat","cells"), suffix = c(".well",".bulk"))

# combine data
frame.plot = rbind(tmp1,tmp2) %>% unique()
dim(frame.plot)

# 1. p values
p1 = ggplot(frame.plot[sample(rownames(frame.plot)),], aes(-log10(p.bulk), -log10(p.well))) + geom_point(size = 2, alpha = 0.85, aes(color = cells)) + #ggrepel::geom_text_repel(aes(label = feat)) +
  theme_bw() + ggthemes::scale_color_tableau() #+ ggtitle("donor and well level associations")
# 2. beta by std error
p2 = ggplot(frame.plot[sample(rownames(frame.plot)),], aes(est.bulk/std.err.bulk, est.well/std.err.well)) + geom_point(size = 2, alpha = 0.85, aes(color = cells)) + #ggrepel::geom_text_repel(aes(label = feat)) +
  theme_bw() + ggthemes::scale_color_tableau() #+ ggtitle("donor and well level associations")

# save plot
pdf("plots/tmp.pdf", width = 6, height = 2)
gridExtra::grid.arrange(p1, p2, nrow = 1)
dev.off()

}
#------------------------------------------------------------
#------------------- cna on single cells --------------------
#------------------------------------------------------------
if(F){
# recode metadata at donor level for cna
grep("Meta", colnames(features.values), value=T)
features.values %>% dplyr::select(Metadata_line_ID,Metadata_Plate,Metadata_Sex,Metadata_onEdge) %>%
  mutate(Metadata_Sex = recode(Metadata_Sex, "Male"=0, "Female"=1)) %>%
  mutate(Metadata_Plate = as.numeric(as.factor(Metadata_Plate))) %>%
  saveRDS(file="cna/meta.bulk.nonRedun.rds")
print("rds of non redundant metadata saved")
  
#data.bulk = readRDS("~/Documents/cmQTL/single.cells/intcol.sampled.rds")
data.bulk = readRDS("rds.objects/singleCell/sampled.data/intcol.sampled.rds")
print("sampled single cells loaded")
print(dim(data.bulk))

# combine donor information
data.bulk = merge(x=data.bulk, y=donor.meta, by.x="Metadata_line_ID", by.y="Metadata_Project.Alias", sort=F)
print("donor info added to data with all features")
print(dim(data.bulk))

# gaussianize each bulk feature across all plates
data.sampled.bulk.gaussian = data.bulk %>% mutate_at(.vars = grep("Meta|Sample", colnames(data.bulk), value=T, ignore.case=T, invert=T), .funs = RNOmni::rankNorm)
print("features gaussianized")
print(dim(data.sampled.bulk.gaussian))
rm(data.bulk)

# remove certain columns
cols.to.remove = grep("sample_|Sample_|_Children_|_Number_|_Parent_|_Center_|_Location_|_Count_|Granularity_1[4-6]|_Orientation|Euler", colnames(data.sampled.bulk.gaussian), value = T)
data.sampled.bulk.gaussian = data.sampled.bulk.gaussian %>% dplyr::select(-one_of(cols.to.remove))
print("certain columns removed")
print(dim(data.sampled.bulk.gaussian))
print(head(data.sampled.bulk.gaussian$Cells_Neighbors_NumberOfNeighbors_Adjacent))

# remove features with near zero variance
ftrs = data.sampled.bulk.gaussian[, grep("Meta", colnames(data.sampled.bulk.gaussian), value = T, invert = T)]
nzv = nearZeroVar(x = ftrs, saveMetrics = F, freqCut = 90/10, allowParallel = T, names = T)
data.sampled.bulk.gaussian = data.sampled.bulk.gaussian %>% dplyr::select(-one_of(nzv))
print("near zero variance features removed")
print(dim(data.sampled.bulk.gaussian))
rm(ftrs)

# data.sampled.bulk.gaussian = readRDS("~/Documents/cmQTL/single.cells/intcol.sampled.gaussian.rds")
# save traits of non-redundant cells
all.meta.cols = grep("Meta", colnames(data.sampled.bulk.gaussian), value=T)
features.values %>% dplyr::select(Metadata_line_ID, Metadata_Plate) %>%
  inner_join(data.sampled.bulk.gaussian, by=c("Metadata_line_ID","Metadata_Plate")) %>%
  dplyr::select(-one_of(all.meta.cols)) %>% saveRDS("cna/cells_traits_nonRedun.rds")
print("rds of non redundant and gaussianized data saved")

# save metadata of non-redundant cells 
features.values %>% dplyr::select(Metadata_line_ID, Metadata_Plate) %>%
  inner_join(data.sampled.bulk.gaussian, by=c("Metadata_line_ID","Metadata_Plate")) %>%
  dplyr::select(Metadata_line_ID,Metadata_Plate,Metadata_Sex) %>% dim
  saveRDS("cna/meta.cells.nonRedun.rds")
print("rds of metadata of non redundant and gaussianized data saved")

# saveRDS(object=data.sampled.bulk.gaussian, file="rds.objects/singleCell/sampled.data/intcol.sampled.gaussian.rds")
# print("rds of gaussianized data saved")
# print(dim(data.sampled.bulk.gaussian))

# data.sampled.bulk.gaussian = readRDS("rds.objects/singleCell/sampled.data/intcol.sampled.gaussian.rds")
# 
# ###
# meta.cols = c("Metadata_line_ID","Metadata_Plate","Metadata_Well","Metadata_Sex","Cells_Neighbors_NumberOfNeighbors_Adjacent","Metadata_Disease_Status")
# data.sampled.bulk.gaussian %>% 
#   dplyr::select(one_of(meta.cols)) %>% saveRDS("cna/meta.rds")
# 
# ###
# all.meta.cols = grep("Meta", colnames(data.sampled.bulk.gaussian), value=T)
# data.sampled.bulk.gaussian %>% 
#   dplyr::select(-one_of(all.meta.cols)) %>% saveRDS("cna/cells_traits.rds")

# PCA on gaussianised traits
# pca.res = data.sampled.bulk.gaussian %>% mutate("Meta_id" = paste(Metadata_Plate,Metadata_Well,Metadata_line_ID,"cell",row_number(),sep = "_")) %>%
#   tibble::remove_rownames() %>% tibble::column_to_rownames("Meta_id") %>%
#   dplyr::select(one_of(grep("Meta|Sample", colnames(data.sampled.bulk.gaussian), invert=T, value=T))) %>% prcomp()
# print("pca done")
# print(dim(pca.res$x))
# 
# pca.res$x %>% as.data.frame() %>% dplyr::select(one_of(paste0("PC",1:10))) %>% tibble::rownames_to_column("Meta_id") %>%
#   readr::write_delim(path="eigenvals.singlecells.tab", delim="\t", col_names=T)
# print("rds of pca values saved")
# 
# pca.res = read.delim("eigenvals.singlecells.tab")
# pca.res %>% mutate(Metadata_line_ID = sub(".*_[A-Z].*_(.*)_cell.*", "\\1", Meta_id),
#                    Metadata_Plate = sub("_.*", "", Meta_id),
#                    Metadata_Well = sub(".*_([A-Z][0-9][0-9])_.*", "\\1", Meta_id)) %>% 
#   ggplot(aes(PC1,PC2)) + geom_point(aes(color=Metadata_Plate), size=0.25) + guides(color=F)

}
#------------------------------------------------------------
#----------------- common and rare variants -----------------
#------------------------------------------------------------
# combine gwas results from samira
if(F){
cores = 50
print("--> combining gwas results...")

pheno.names = read.table("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/SecondRound_JatinFeatures/intermediate/intermediate.phenotypeNames.txt", header=F, col.names=c("idx","name"))
all.res.files = list.files("/data/srlab2/sasgari/Projects/V2F/Analysis/GWAS/SecondRound_JatinFeatures/intermediate", 
                        pattern="chr.*pheno.*.LR.fastGWA.gz", full.names=T)
print(length(all.res.files))
print(head(all.res.files,2))
block = split(1:length(all.res.files), cut_interval(1:length(all.res.files), n=cores))
gwas.files = all.res.files[block[[seed]]]
print(paste("working on",length(gwas.files)))

gwas.res = NULL
for(file in gwas.files){
  print(file)
  if(file.info(file)$size > 500){
    feat = pheno.names %>% filter(idx == as.numeric(sub(".*chr.*_pheno([0-9]{1,}).LR.*", "\\1", file))) %>% pull(name) %>% as.character()
    gwas.res = rbind(gwas.res, read.delim(file) %>% filter(P < 10^-5) %>% mutate("feat" = feat))
  }
}
readr::write_delim(x=gwas.res, file=paste("gwas/intermediate",seed,"tab", sep="."), delim="\t")
print("gwas results saved")

###
gwas.files = list.files("gwas", pattern=".tab", full.names=T)
print(length(gwas.files))
print(head(gwas.files,2))
gwas.res = gwas.files %>% purrr::map_df(read.delim)

###
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
attributes = c("ensembl_gene_id","start_position","end_position","hgnc_symbol","chromosome_name")
genes.pos = getBM(attributes=attributes, filters=c("hgnc_symbol"), values=list(hgnc_symbol=unique(sig.assoc.p$var)), 
                  mart=mart, uniqueRows = T) %>% filter(chromosome_name %in% 1:22)
sig.assoc.p = sig.assoc.p %>% inner_join(genes.pos, by=c("var"="hgnc_symbol"))

###
sig.assoc.p.snp = NULL
for(i in 1:nrow(sig.assoc.p)){
  # s = sig.assoc.p %>% filter(chromosome_name == gwas.res.th.p[i,"CHR"] & start_blc < gwas.res.th.p[i,"POS"] & end_blc > gwas.res.th.p[i,"POS"])
  s = gwas.res %>% filter(CHR == as.numeric(sig.assoc.p[i,"chromosome_name"]) & feat == sig.assoc.p[i,"feat"] &
                                 POS >= (sig.assoc.p[i,"start_position"]-1000000) & POS <= (sig.assoc.p[i,"end_position"]+1000000))
  if(nrow(s) > 0) sig.assoc.p.snp = rbind(sig.assoc.p.snp, cbind(sig.assoc.p[i,] %>% dplyr::select(feat,var,est,p), s %>% dplyr::select(SNP,POS,BETA,P)))
  if(nrow(s) > 0) print(i)
}

# get LD information
# tab = LDlinkR::LDpair(var1="chr20:29297594", var2="chr7:29297594", token="df65c945c99c")
# tab = LDlinkR::LDmatrix(c("rs3", "rs4", "rs148890987"), "YRI", "r2", token="df65c945c99c")
ld.info = LDlinkR::LDproxy(snp="rs148890987", pop = c("CEU","AFR","AMR","ASW"), r2d="d", token="df65c945c99c")

df.hg19 = ld.info %>% filter(Dprime > 0.9) %>% mutate(chr=sub("\\:.*","",Coord), start=sub(".*\\:","",Coord), end=sub(".*\\:","",Coord)) %>%
  dplyr::select(chr,start,end)
df.hg38 = DeepBlueR::deepblue_liftover(regions=GenomicRanges::makeGRangesFromDataFrame(df=df.hg19), 
                             source="hg19", target="hg38")

}
#------------------------------------------------------------
# association feature value and disease status
if(F){
# merge feature values with pcs info
data = features.values %>% inner_join(pcs.terra, by = c("Metadata_line_ID" = "Metadata_Project.Alias"))
print("data prepared")
print(dim(data))

# select covariates - the cell cycle is Nuclei_Intensity_IntegratedIntensity_DNA
covars = c("Metadata_Plate", "Metadata_iPSC_Origin", "Metadata_Sex", "Metadata_Disease_Status", paste0("PC",1:5))
if(level == "well") covars = union(covars, c("Metadata_Well", "Metadata_line_ID"))
if(round != "isolate") covars = union("Cells_Neighbors_NumberOfNeighbors_Adjacent", covars)

res = NULL
for(feature in features.to.test){
  if(level == "well"){
    effects = paste(c(sub("(Meta.*)","(1|\\1)", grep("Disease_Status|Age", covars, value=T, invert=T)), grep("Age", covars, value=T)), collapse = " + ")
    baseline = as.formula(paste(feature,"~", "Metadata_Disease_Status", "+", effects))
    lm.b = lmer(formula = baseline, data = data, control = lmerControl(optimizer = "bobyqa"), REML=F)
    ###
    pval = car::Anova(lm.b) %>% as.data.frame() %>% tibble::rownames_to_column("variable") %>%
      filter(variable == "Metadata_Disease_Status") %>% pull("Pr(>Chisq)")
    ###
    tmp = coefficients(summary(lm.b)) %>% as.data.frame() %>% tibble::rownames_to_column("variable") %>%
      filter(variable == "Metadata_donor_Affected.StatusYes") %>% mutate("feature" = !!feature, "p" = pval) %>%
      select(one_of("variable","Estimate","p","feature"))
  }
  if(level == "bulk"){
    baseline = as.formula(paste(feature,"~", "Metadata_Disease_Status", "+", paste(grep("Disease_Status", covars, value=T, invert=T), collapse = "+")))
    lm.b = lm(formula = baseline, data = data)
    tmp = coefficients(summary(lm.b)) %>% as.data.frame() %>% tibble::rownames_to_column("variable") %>%
      filter(grepl("Disease_Status",variable)) %>% mutate("feature" = !!feature) %>% dplyr::rename("p" = "Pr(>|t|)") %>% filter(p < 0.05) %>%
      dplyr::select(one_of("variable","Estimate","p","feature"))
  }
  res = rbind(res, tmp)
  print(feature)
}

res %>% mutate(p.adj = p.adjust(p,"bon",length(features.to.test))) %>% arrange(p) %>% head()

}
#------------------------------------------------------------
# remove batch effect from features and get extreme values
if(F){
print(paste("looking at extreme feature values for", round))
  
# regress out plate, image-batch, well etc. covariates
data.model = data.sampled.bulk.gaussian[!is.na(data.sampled.bulk.gaussian$Metadata_donor_Age), ]
print("donors who age is unknown removed from raw data")
print(dim(data.model))

# split features for parallel processing
features.to.test = grep("Meta|_Neighbors_", colnames(data.model), invert = T, value = T)
block = split(1:length(features.to.test), cut_interval(1:length(features.to.test), n = cores))
print(paste("testing",length(features.to.test[block[[seed]]])))

# variance in feature accounted by fixed and random effects
extreme.res = NULL
residuals.res = NULL
for(feature in features.to.test){
  
  # prepare a mixed model
  print(feature)
  # rand.effects = "(1|Metadata_Plate) + (1|Metadata_PlateWell) + (1|Metadata_Well) + (1|Metadata_image_batch) + (1|Metadata_donor_Sex) + (1|Metadata_donor_Source) + (1|Metadata_donor_Race) + (1|Metadata_donor_Affected.Status) + (1|Metadata_line_ID)"
  rand.effects = "(1|Metadata_Plate) + (1|Metadata_image_batch)"
  fixed.effects = c("Cells_Neighbors_NumberOfNeighbors_Adjacent")
  mix.effc.formula = as.formula(paste(feature, "~" , fixed.effects, "+", rand.effects))
  
  # calculate mixed effect models and save residuals
  model.effc = lmer(formula = mix.effc.formula, data = data.model, control = lmerControl(optimizer = "bobyqa"), REML = F)
  frame = data.frame("feature" = residuals(model.effc), data.model[, grep("Meta", colnames(data.model), value = T)])
  colnames(frame)[1] = feature

  # store residuals for burden test
  if(is.null(residuals.res)) residuals.res = frame else residuals.res = merge(x = residuals.res, y = frame, by = grep("Meta", colnames(data.model), value = T))
  print("residuals saved")
  
  # which cell lines has extreme value for selected features
  tmp = frame %>% mutate( zscore = (!!as.name(feature) - mean(!!as.name(feature)))/sd(!!as.name(feature)) ) %>% dplyr::filter(abs(zscore) > 3 & abs(zscore) == max(abs(zscore))) %>% 
    dplyr::select(c(Metadata_line_ID, Metadata_Plate, Metadata_donor_Description, Metadata_donor_Sex, Metadata_donor_Source, Metadata_donor_Affected.Status, Metadata_donor_Age))
  if(nrow(tmp) > 0) extreme.res = rbind(extreme.res, data.frame(tmp, feature))
  
}

# save donor info with extreme feature values
saveRDS(object = extreme.res, file = paste0("rds.objects/singleCell/others/extremes/extreme.res.", seed, ".rds"))
print(paste("rds of extreme res saved for", seed))
print(dim(extreme.res))

# get residuas for donor with wgs
residuals.res.donor = merge(x = residuals.res, y = donor.meta[, c("Metadata_donor_Project.Alias", "Sample_ID")], by.x = "Metadata_line_ID", by.y = "Metadata_donor_Project.Alias")
residuals.res.donor.wgs = merge(x = residuals.res.donor, y = wgs.good, by = "Sample_ID")
print(dim(residuals.res.donor.wgs))

# map pcs+terra to wgs sample id
pcs.terra.wgs = residuals.res.donor.wgs %>% select(c(Metadata_line_ID, Sample)) %>% inner_join(pcs.terra, by = c("Metadata_line_ID" = "Project.Alias"))

# save covariates in ped format
covar = data.frame("fid" = pcs.terra.wgs$Sample, "iid" = pcs.terra.wgs$Sample, "fatid" = 0, "matid" = 0, "sex" = 0, pcs.terra.wgs[, paste0("PC",1:5)])
print(dim(covar))
write.table(x = covar, file = "extreme.values/covar.ped", append = F, quote = F, sep = " ", row.names = F, col.names = T)

# take those sample with covars
residuals.res.donor.wgs = residuals.res.donor.wgs %>% filter(Metadata_line_ID %in% pcs.terra.wgs$Metadata_line_ID)

# save residuals as phenotype in ped format
tfam = data.frame("fid" = residuals.res.donor.wgs$Sample, "iid" = residuals.res.donor.wgs$Sample, "fatid" = 0, "matid" = 0, "sex" = 0, residuals.res.donor.wgs[, grep("Meta", colnames(residuals.res), value = T, invert = T)])
print(dim(tfam))
write.table(x = tfam, file = "extreme.values/pheno.ped", append = F, quote = F, sep = " ", row.names = F, col.names = T)


# terra = terra[complete.cases(terra), ]
# terra = terra[!duplicated(terra$Project.Alias), ]
# a = merge(x = terra[, c("IID", "Project.Alias")], y = pcs, by.x = "IID", by.y = "V1")
# terra %>% filter(Project.Alias %in% setdiff(residuals.res.donor.wgs$Metadata_line_ID, a$Project.Alias))
# b = merge(x = a, y = residuals.res.donor.wgs[, c("Metadata_line_ID", "Sample")], by.x = "Project.Alias", by.y = "Metadata_line_ID", all = F)

}
# enrichment of extreme values
if(F){
# bring in raw sampled data
data.model = data.sampled[!is.na(data.sampled$Metadata_donor_Age), ]
print("donors who age is unknown removed from raw data")
print(dim(data.model))
# data.model$Metadata_donor_disease = sub("(^CONTROL).*", "\\1", data.model$Metadata_donor_Description)

# format metadata
### c("ALZHEIMER|AMD|ASD|ID","EPILEPSY|CEREBRAL","HCV","BLINDING|GLAUCOMA","DIABETIC|DR|DIABETES","LIVER","CARDIO")
donor.meta = donor.meta %>% filter(Metadata_donor_Project.Alias %in% data.model$Metadata_line_ID)
donor.meta$disease = sub("(^CONTROL).*", "\\1", donor.meta$Metadata_donor_Description)
donor.meta = donor.meta %>% mutate(disease_groups = ifelse( grepl("ALZHEIMER|AMD|ASD|ID", disease), "physc", ifelse( grepl("DIABETIC|DR|DIABETES", disease), "sugar", "others")) )
print(dim(donor.meta))

# load extreme value stats
extreme.res.files = list.files(path = "rds.objects/singleCell/others/extremes", pattern = "*.rds", full.names = T)
print(paste("number of files", length(extreme.res.files)))
extreme.res = NULL
for(file in extreme.res.files){
  tmp = readRDS(file)
  extreme.res = rbind(extreme.res, tmp)
}

# pre-process extreme value stats
extreme.res$disease = sub("(^CONTROL).*", "\\1", extreme.res$Metadata_donor_Description)
extreme.res = extreme.res %>% mutate(disease_groups = ifelse( grepl("ALZHEIMER|AMD|ASD|ID", disease), "physc", ifelse( grepl("DIABETIC|DR|DIABETES", disease), "sugar", "others")) )
print(paste("extreme res: number of cell lines", length(unique(extreme.res$Metadata_line_ID))))
print(paste("extreme res: number of features", length(unique(extreme.res$feature))))

# plot metadata in hist bins
# overall
meta = "Metadata_donor_Sex"
extreme.res %>% group_by(Metadata_line_ID) %>% mutate(count = n()) %>% dplyr::select(c(count, Metadata_line_ID, !!as.name(meta))) %>% distinct() %>%
  ggplot(aes(count)) + geom_bar(width = 0.8, aes(fill = !!as.name(meta))) + scale_x_binned(n.breaks = 10) + scale_fill_viridis_d(option = "C") + theme_bw() +
  xlab("#extreme values")
# zoomed in
extreme.res %>% group_by(Metadata_line_ID) %>% mutate(count = n()) %>% filter(count > 50) %>% group_by(Metadata_line_ID) %>% mutate(count = n()) %>% dplyr::select(c(count, Metadata_line_ID, !!as.name(meta))) %>% distinct() %>%
  ggplot(aes(count)) + geom_bar(width = 0.8, aes(fill = !!as.name(meta))) + scale_x_binned(n.breaks = 10) + scale_fill_viridis_d(option = "C") + theme_bw() +
  xlab("#extreme values")

extreme.ids = extreme.res %>% group_by(Metadata_line_ID) %>% mutate(count = n()) %>% filter(count > 50) %>% .$Metadata_line_ID %>% unique()

# tfam = data.frame("fid" = wgs.good$Sample, "iid" = wgs.good$Sample, "fatid" = 0, "matid" = 0, "sex" = 0)
# wgs.good %>% filter(Sample_ID %in% donor.meta$Sample_ID) %>% dim()

# donor.meta %>% filter(Sample_ID %in% wgs.good$Sample_ID) %>% .$Metadata_donor_Project.Alias %>%
#   {filter(residuals.res, Metadata_line_ID %in% .) %>% head()} 

# tfam.pheno = donor.meta %>% filter(Metadata_donor_Project.Alias %in% extreme.ids) %>% .$Sample_ID %>%
#   {filter(wgs.good, Sample_ID %in% .) %>% .$Sample} %>% {mutate(tfam, phenotype = ifelse(Sample %in% ., "2", "1"))}

# fid iid fatid matid sex y1 y2 y3 y4
# P1 P1 0 0 0 1.7642934435605 -0.733862638327895 -0.980843608339726 2
# P2 P2 0 0 0 0.457111744989746 0.623297281416372 -2.24266162284447 1

  # Family ID ('FID')
  # Within-family ID ('IID'; cannot be '0')
  # Within-family ID of father ('0' if father isn't in dataset)
  # Within-family ID of mother ('0' if mother isn't in dataset)
  # Sex code ('1' = male, '2' = female, '0' = unknown)
  # Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
  # Phenotype column should be coded as {0,1,2} meaning {Missing, Control, Case}
  
  
  
  # summary of cell lines having extreme values
  a = extreme.res %>% group_by(Metadata_line_ID) %>% mutate(count = n()) %>% filter(count > 50) %>% .$Metadata_donor_Source %>% table()
  b = extreme.res %>% group_by(Metadata_line_ID) %>% mutate(count = n()) %>% filter(count < 50) %>% .$Metadata_donor_Source %>% table()
  rbind(a,b)
  fisher.test(rbind(a,b))
  
  
  a = extreme.res %>% group_by(Metadata_line_ID) %>% mutate(count = n()) %>% filter(count > 50) %>% .$Metadata_line_ID %>% unique() %>%
    {filter(donor.meta, Metadata_donor_Project.Alias %in% .) %>% .$Metadata_donor_Source %>% table()}
  b = table(donor.meta$Metadata_donor_Source)
  
  # plot metadata in bars
  # frame = extreme.res %>% group_by(Metadata_line_ID) %>% mutate(count = n()) %>% filter(count > 10) %>% dplyr::select(c(Metadata_line_ID,count, Metadata_donor_Source)) %>% dplyr::distinct() %>% data.frame()
  # frame$Metadata_line_ID = factor(frame$Metadata_line_ID, frame[order(frame$count), "Metadata_line_ID"])
  # ggplot(frame, aes(Metadata_line_ID, count)) + geom_col(width = 0.5, aes(fill = Metadata_donor_Source)) + theme_bw() + theme(axis.text.x = element_blank()) +
  #   scale_fill_manual(values = col.arr, guide = F)
  
  
  # 
  # ###
  # extreme.res.good = extreme.res %>% filter(Metadata_line_ID != out.lines)
  # donor.meta.good = donor.meta %>% filter(Metadata_donor_Project.Alias != out.lines)
  # for(meta in c("Metadata_donor_Source","Metadata_donor_Affected.Status","Metadata_donor_Sex")){
  # 
  #   c = donor.meta.good[grep(pattern, donor.meta.good$disease, invert = T), c("Metadata_donor_Project.Alias")]
  #   d = extreme.res.good %>% filter(Metadata_line_ID %in% c) %>% .$Metadata_line_ID %>% unique()
  #   
  #   mat = rbind(table(extreme.res.good[, meta]), table(donor.meta.good[, meta]))
  #   ftest = fisher.test(mat)
  #   ftest$p.value
  #   print(data.frame(meta, or = round(ftest$estimate,3), p = round(ftest$p.value,3)))
  # }
  # 
  # ###
  # extreme.res.good = extreme.res %>% group_by(Metadata_line_ID) %>% mutate(count = n()) %>% filter(count > 10)
  # donor.meta.good = donor.meta %>% filter(Metadata_donor_Project.Alias %in% extreme.res.good$Metadata_line_ID)
  # 
  # pval.arr = NULL
  # for(pattern in c("ALZHEIMER|AMD|ASD|ID","EPILEPSY|CEREBRAL","HCV","BLINDING|GLAUCOMA","DIABETIC|DR|DIABETES","LIVER","CARDIO")){
  #   
  #   a = donor.meta.good[grep(pattern, donor.meta.good$disease), c("Metadata_donor_Project.Alias")]
  #   b = extreme.res.good %>% filter(Metadata_line_ID %in% a) %>% .$Metadata_line_ID %>% unique()
  #   
  #   c = donor.meta.good[grep(pattern, donor.meta.good$disease, invert = T), c("Metadata_donor_Project.Alias")]
  #   d = extreme.res.good %>% filter(Metadata_line_ID %in% c) %>% .$Metadata_line_ID %>% unique()
  # 
  #   mat = matrix( c(length(a), length(b), length(c), length(d)), nrow = 2)
  #   ftest = fisher.test(mat)
  #   pval.arr = c(pval.arr, ftest$p.value)
  #   if(ftest$p.value < 0.05) print(data.frame(pattern, or = round(ftest$estimate,3), p = round(ftest$p.value,3)))
  #   
  # }
  
  # a = sub(".*;[[:space:]]", "", donor.meta$Metadata_donor_Description)
  # b = sub("[|][[:space:]](.*),[[:space:]].*", "\\1", a)
  # c = sub("[[:space:]][|].*", "", sub("(CONTROL).*", "\\1", b))
  # frame.total = data.frame(table(sub("[,][[:space:]].*", "", c)))
  
  
  
  # extreme.res$disease = sub("(^CONTROL).*", "\\1", extreme.res$Metadata_donor_Description)
  # a = sub("(^CONTROL).*", "\\1", unique(extreme.res$Metadata_donor_Description)) 
  # b = sub(".*?;[[:space:]](.*)([|].*;[[:space:]](.*))*", "\\1 \\3", a)
  # sub("[,][[:space:]].*", "", sub("[[:space:]][|]", "", a))
  # 
  # a = sub(".*;[[:space:]](.*)[|].*;[[:space:]](.*)", "\\1 \\2", )
  # b = sub("[|][[:space:]](.*),[[:space:]].*", "\\1", a)
  # c = sub("[[:space:]][|].*", "", sub("(CONTROL).*", "\\1", b))
  # extreme.res$disease = sub("[,][[:space:]].*", "", c)
  # 
  # a = sub(".*;[[:space:]]", "", extreme.res$Metadata_donor_Description)
  # b = sub("[|][[:space:]](.*),[[:space:]].*", "\\1", a)
  # c = sub("[[:space:]][|].*", "", sub("(CONTROL).*", "\\1", b))
  # frame = data.frame(table(sub("[,][[:space:]].*", "", c)))
  
  # merge(frame, frame.total, by = "Var1", suffixes = c("_extremes","_donors")) %>% ggplot(aes(Freq_donors, Freq_extremes)) + geom_point() + theme_bw() +
  #   geom_abline(intercept = 0 , slope = 1)
  # 
  # for (var in unique(m$Var1)){
  #   r1 = m[m$Var1 == var, -1]
  #   r2 = colSums(m[m$Var1 != var, -1])
  #   pval = fisher.test(rbind(r1, r2))$p.value
  #   if(pval < 0.05) print(c(var, pval))
  # }
  
  # b = unique(data.sampled[grep(pattern, data.sampled$Metadata_donor_disease), "Metadata_line_ID"])
  # d = unique(data.sampled[grep(pattern, data.sampled$Metadata_donor_disease, invert = T), "Metadata_line_ID"])
  
}
# summary of rare variants in 1kg and cmqtl
if(F){
#---------------------
# freq of rare 1kg variants in cmQTL
freq.cm.in.kg = read.delim("wgs/SCBB_CIRM_WGScallset_2020_v1.bi.kg.freq.afreq")
# plot
ggplot(freq.cm.in.kg, aes(ALT_FREQS)) + geom_histogram(bins = 100) + theme_bw() + ggtitle("freq of rare 1kg variants\nin cmqtl") +
  annotate(geom = "text", x=0.75, y=Inf, size=5, vjust = 2, label = nrow(freq.cm.in.kg)) 
# zoomed plot
freq.cm.in.kg %>% filter(ALT_FREQS < 0.1) %>% ggplot(aes(ALT_FREQS)) + geom_histogram(bins = 30) + theme_bw() + 
  geom_vline(xintercept = 0.01, color = "red")

#---------------------
# variants rare in 1kg and cmqtl
donor.gt = read.delim("wgs/rareKG.rareCM.high.tab", header = F, stringsAsFactors = F) %>% slice(-1)
colnames(donor.gt) = read.delim("bekar.cheeze/wgs_cols.txt")$Cols
hist(as.numeric(donor.gt$AF), breaks = seq(0,1,0.01), main = "rare variants in cmqtl")

#---------------------
# rare.kg = read.delim("wgs/kg.high.tab", header = F, stringsAsFactors = F) %>% slice(-1)
# rare.maf1 = read.delim("wgs/maf1.high.tab", header = F, stringsAsFactors = F) %>% slice(-1)
# # combine 1kg variants with maf1 in cmQTL
# rare.kg.maf1 = rare.kg[,c("V3","V9")] %>% inner_join(rare.maf1[,c("V3","V9")], by = "V3", suffix = c("_kg","_maf1"))
# rare.kg.maf1 %>% dim()

# # combine 1kg variants with maf1 in cmQTL
# n.out.1 = rare.kg.maf1 %>% filter(V9_kg > 0.01 & V9_maf1 < 0.01) %>% nrow()
# n.out.2 = rare.kg.maf1 %>% filter(V9_kg < 0.01 & V9_maf1 > 0.01) %>% nrow()
# ggplot(rare.kg.maf1, aes(as.numeric(V9_kg), as.numeric(V9_maf1))) + geom_point(size=3) + theme_bw() +
#   annotate(geom = "text", x=0.75, y=0.75, size=5, label = paste0(nrow(rare.kg.maf1), "\n", n.out.1)) + 
#   xlab("freq rare variants in 1kg") + ylab("freq rare variants in cmqtl") + stat_binhex() + scale_fill_viridis_c(option = "D") +
#   theme(axis.title = element_text(size = 15), axis.text = element_text(size = 13)) + ggtitle("freq of cmqtl\nrare variants in 1kg")

}
#------------------------------------------------------------------------------------------------------------------
print("there were following warnings")
print(warnings())
print("rscript done")
#------------------------------------------------------------------------------------------------------------------
if(F){

  # outliers having a lot of variants and outliers in PC space
  # donor.meta = donor.meta %>% filter(!(Metadata_ID %in% c("CW10054","CW10162","CW10100","CW10101","CW10163","CW10095","CW20252")))
  
  # if(level != "single") all.data.files = list.files(path = paste0("rds.objects/",level, "/archive"), pattern = paste0(round,".*donor.rds"), full.names = T)
  
  # # get uncorrelated features for burden testing
  # ftrs = data.sampled.bulk.gaussian %>% as.data.frame() %>% dplyr::select(one_of( replicate_correlation_values %>% filter(median >= 0.5) %>% pull(variable) ))
  # cv = findCorrelation(x = cor(ftrs), cutoff = 0.9, names = T, exact = F)
  # features.uncor = data.frame(feature = setdiff(colnames(ftrs), cv))
  # print(paste("uncorrelated features taken for downsteam testing", nrow(features.uncor)))
  
  # #-----------------------
  # #--------- tmp ---------
  # # comparison of gaussianisation at donor level with well level then mean
  # feature = c("Nuclei_Granularity_1_DNA","Cytoplasm_AreaShape_Zernike_9_1","Cells_AreaShape_Area")[2]
  # 
  # # gaussianisation at donor level
  # a = data.sampled.bulk.gaussian.bulk %>% dplyr::select(Metadata_line_ID,!!as.name(feature))
  # a = a[!duplicated(a$Metadata_line_ID), ]
  # 
  # # gaussianisation at well level
  # b = data.sampled.bulk.gaussian.from.well %>% as.data.frame() %>%  dplyr::select(Metadata_line_ID,!!as.name(feature))
  # b = b[!duplicated(b$Metadata_line_ID), ]
  # 
  # # comparison
  # d = a %>% inner_join(b, by="Metadata_line_ID", suffix = c(".donor", ".well"))
  # cr = cor.test(d[,paste0(feature,".donor")], d[,paste0(feature,".well")])
  # plot(d[,paste0(feature,".donor")], d[,paste0(feature,".well")], main=paste(feature,"\n r = ",round(cr$estimate,2)),
  #      xlab="donor level gaussianized", ylab="well level gaussianized, then median")
  # #-----------------------
  # #--------- tmp ---------
  
  
  # effects = paste(c(sub("(Meta.*)","(1|\\1)", grep(paste0("Age|CellCount|",inter.var), covars, value = T, invert = T)), grep("Age|CellCount", covars, value = T)), collapse = " + ")
  # baseline = as.formula(paste(feature,"~", gene,"+", inter.var, "+", effects))
  # alternative = as.formula(paste(feature,"~", gene,"*", inter.var, "+", effects))
  
  # distribution of feature values (n=224)
  features.values %>% dplyr::select(one_of(features.to.test)) %>% reshape2::melt() %>% 
    ggplot(aes(value)) + geom_histogram(bins=50) + facet_wrap(~variable, ncol=5)
  
  # effective number of test or traits
  ftrs = data.sampled.bulk.gaussian %>% dplyr::select(one_of(features.to.test))
  effective.traits = poolr::meff(R=cor(ftrs), method="gao")
  burden.res.bulk %>% filter(p < 0.00000156/effective.traits)
  
  d = readr::read_delim("~/Google Drive/cmQTL_Jatin_Samira/intermediate.trait.gene.p.mat", delim="\t")
  
  # combine p value for all genes per trait
  e = d %>% tibble::column_to_rownames("feat") %>% t() %>% as.data.frame() %>% 
    summarise_at(.vars=colnames(.), .funs=ACAT::ACAT)
  e %>% t() %>% as.data.frame() %>% filter(V1 < 0.05/224)
  
  # combine p value for all traits per gene
  f = d %>% summarise_at(.vars=setdiff(colnames(d),"feat"), .funs=ACAT::ACAT)
  f %>% t() %>% as.data.frame() %>% filter(V1 < 0.05/9033)
  
  res = NULL
  for(gene in setdiff(colnames(d),"feat")){
    res = rbind( res, d %>% dplyr::select(feat,!!as.name(gene)) %>% 
                   mutate("p2" = p.adjust(!!as.name(gene), method="BH", n=9033*224)) %>% 
                   arrange(p2) %>% filter(p2 < 0.05) %>% dplyr::select(feat,p2) %>%
                   mutate("gene" = gene) )
  }
  
  
  # b = data.sampled.bulk.gaussian %>% group_by(Metadata_line_ID,Metadata_Plate) %>%
  #   summarise_at(.vars=grep("Meta", colnames(data.sampled.bulk.gaussian), value=T, invert=T, ignore.case=T), .funs=mean) %>% unique()
  # 
  # d = data.sampled.bulk.gaussian %>% group_by(Metadata_line_ID) %>% dplyr::select(-Metadata_Well) %>%
  #   mutate_at(.vars=grep("Meta", colnames(data.sampled.bulk.gaussian), value=T, invert=T, ignore.case=T), .funs=mean) %>% unique()
  
  
  # # ---- temp ----
  # frame.to.plot = res %>% arrange(p) %>% head(1) %>% dplyr::select(feat,var,est,p) %>%
  #   mutate(p = format(p, digits=1, sc=T), variable = sub("_", "-", sub("_inter","",var)), assoc = ifelse(grepl("inter",var), "inter", "gene"))
  # # ---- temp ----
  
  # # -------------------- interactions -----------------
  # # interaction
  # all.res.files = list.files("extreme.values/burden.out", pattern = paste("rv",round,level,phenotype,sum.or.present,"*.interaction",variant.source,inter.var,"tab",sep="."), full.names = T)
  # print(length(all.res.files))
  # burden.res.inter <- all.res.files %>% purrr::map_df(read.delim)
  # 
  # # significant interactions
  # frame.to.plot = burden.res.inter %>% filter(data == "orig" & p.adj < 0.1) %>% dplyr::select(feat,var,est,p) %>% 
  #   mutate(p = format(p, digits=1, sc=T), variable = sub("_", "-", sub("_inter","",var)), assoc = ifelse(grepl("inter",var), "inter", "gene"))
  
  
  # ----- 1. select associations based on FDR -----
  # critical value
  burden.res.bulk.fdr = burden.res.bulk %>% filter(var %in% genes.to.test$gene.2) %>% add_tally() %>% arrange(p) %>% 
    mutate("i" = row_number(), "cv" = (i/n)*0.2, "p_cv" = p < cv)
  
  # potential associations at fdr
  potent.assoc = burden.res.bulk.fdr %>% filter(p_cv == T) %>% tail(1) %>%
    {filter(burden.res.bulk.fdr, row_number() %in% 1:.$i)}
  print(paste("number of associations at fdr 0.2 is", nrow(potent.assoc)))
  
  # # select significant associations using feature-gene-specific background
  # # --- temp ---
  # cores = 10
  # print("--> combining permutation results...")
  # all.res.files = list.files(paste0("/data/srlab/jatin/cmqtl/permut/",round), pattern=paste(level,phenotype,sum.or.present,gene.freq,"*.baseline",variant.source,"tab",sep="."), full.names=T)
  # print(length(all.res.files))
  # print(head(all.res.files,2))
  # 
  # block = split(1:length(all.res.files), cut_interval(1:length(all.res.files), n=cores))
  # files = all.res.files[block[[seed]]]
  # print(paste("working on",length(files)))
  # 
  # # filter permutation results
  # for(file in files){
  #   print(file)
  #   readr::read_delim(file, delim="\t") %>% filter(feat %in% potent.assoc$feat & var %in% potent.assoc$var) %>%
  #     dplyr::select(feat,var,est,p,perm.seed) %>% unique() %>%
  #     readr::write_delim(sub(".*[/](.*).tab","extreme.values/permut.out/\\1.02.fdr",file), delim="\t")
  # }
  # print("quitting after all...")
  # quit()
  # # --- temp ---
  
  # # or load permutation associations
  # setwd("~/Documents/cmQTL/permut")
  # burden.res.perm = readRDS("burden_perm_intermediate_0.02.res")
  # setwd("/Volumes/ja286/cmQTL")
  
  # load permutation results
  all.res.files = list.files("extreme.values/permut.out", pattern=paste0(round,".",level,".*.",gene.freq,".*.baseline.",variant.source,".02.fdr"), full.names = T)
  length(all.res.files)
  head(all.res.files,2)
  burden.res.perm <- all.res.files %>% purrr::map_df(read.delim)
  
  # plot
  plt.list = NULL
  for(i in 1:nrow(potent.assoc)){
    plt.list[[i]] = burden.res.perm %>% filter(feat == potent.assoc[i,"feat"] & var == potent.assoc[i,"var"]) %>% 
      ggplot(aes(abs(est))) + geom_histogram(bins=50) + ggtitle(paste(potent.assoc[i,"feat"], potent.assoc[i,"var"], sep=":")) + 
      geom_segment(data=potent.assoc[i,], mapping=aes(x=abs(est), y=0, xend=abs(est), yend=5), arrow=arrow(length=unit(0.03,"npc"), ends="first"), color="red")
  }
  cowplot::plot_grid(plotlist=plt.list, nrow=5)
  
  
  
  
  
  burden.res.perm.pmin %>% arrange(p.min) %>% head(20)
  burden.res.perm.pmin %>% group_by(perm.seed) %>% summarize_at("p.min", min) %>%
    arrange(p.min) %>% head(50) %>% tail(1)
  
  # # effects.map = c("Metadata_Plate" = "Plate", "Metadata_Well" = "Well",
  # #                 "Metadata_line_ID" = "Cell line", "Metadata_donor_Race" = "Donor race", "Metadata_donor_Source" = "iPSC source",
  # #                 "Metadata_donor_Affected.Status" = "Donor disease", "Metadata_PlateWell" = "PlateWell",
  # #                 "Metadata_donor_Sex" = "Donor sex", "Metadata_donor_Ethnicity" = "Donor ethnicity", "Metadata_donor_Age" = "Donor age",
  # #                 "Metadata_image_batch" = "Imaging batch", "Cell neighbors" = "Cell neighbors")
  # # effect.cols = my.cols[1:length(as.character(effects.map))]
  # # names(effect.cols) = as.character(effects.map)
  
  # over-represenation of a feature class in associated features
  assoc.feat.summ = table(sub(".*_(.*?)_.*","\\1", c(sig.assoc.p$feat))) %>% as.data.frame()
  allfeat.summ = table(sub(".*_(.*?)_.*","\\1", features.to.test)) %>% as.data.frame()
  
  class="AreaShape"
  a = assoc.feat.summ %>% filter(Var1 == !!class) %>% pull(Freq)
  b = assoc.feat.summ %>% filter(Var1 != !!class) %>% pull(Freq) %>% sum()
  c = allfeat.summ %>% filter(Var1 == !!class) %>% pull(Freq)
  d = allfeat.summ %>% filter(Var1 != !!class) %>% pull(Freq) %>% sum()
  
  cont.mat = matrix(c(a,b,c,d), 2, byrow=T, dimnames=list(c("associated","not_associated"),c("specific_class","other_classes")))
  fisher.test(cont.mat)
  
  
  cell.count.summ = list.files("metadata", pattern="cell.category.*.summ", full.names=T) %>% purrr::map_df(read.delim) %>%
    filter(Metadata_line_ID %in% donor.meta$Metadata_Project.Alias)

  cell.count.summ %>% ggplot(aes(n)) + geom_histogram(bins=50, aes(fill=cell.cat), color=NA) + facet_wrap(~cell.cat, scales="free") +
    ggthemes::scale_fill_tableau() + theme_bw() + scale_x_continuous(breaks=scales::extended_breaks(n=7)) + 
    labs(x="cell count for a cell line") + 
    geom_segment(data=cell.count.summ %>% filter(Metadata_line_ID %in% lines), aes(x=n, xend=n), y=0, yend=2)

  # ## --- temp - group by acat ---
  # rand = readRDS("tmp/intermediate.rand.5files.rds")
  # rand %>% filter(var %in% genes.to.test$gene.2) %>% group_by(var) %>% summarise("p.acat" = ACAT::ACAT(p)) %>% arrange(p.acat) %>%
  #   ggplot(aes(p.acat)) + geom_histogram(aes(fill=..count..), stat="count") + ggtitle("all features combined per gene") + 
  #   theme_bw() + scale_x_binned(n.breaks=1/0.05) + labs(x="combined p from acat")
  # 
  # potent.assoc %>% filter(var == "PRLR") %>% arrange(p) %>% head(2)
  # ## --- temp - group by acat ---
  

  # # ----------------
  # # ----- temp -----
  # if(level == "well") features.values = features.values %>% inner_join(cell.count.summ.well, by=c("Metadata_Plate","Metadata_Well")) %>%
  #   mutate_at(.vars = "Metadata_CellCount_avg", .funs = RNOmni::rankNorm)
  # if(level == "bulk") features.values = features.values %>% inner_join(cell.count.summ.line, by="Metadata_line_ID") %>%
  #   mutate_at(.vars = "Metadata_CellCount_avg", .funs = RNOmni::rankNorm)
  # print("cell count appended and gaussianized")
  # # ----- temp -----
  # # ----------------
  
  
  # # significant associations from baseline model
  # all.res.files = list.files("extreme.values/burden.out", pattern = paste("rv",round,level,phenotype,sum.or.present,gene.freq,".*.baseline",variant.source,"tab",sep="."), full.names = T)
  # length(all.res.files); head(all.res.files,2)
  # burden.res.bulk <- all.res.files %>% purrr::map_df(read.delim) %>% filter(data == "orig")
  # # control for fdr
  # burden.res.bulk.fdr = burden.res.bulk %>% filter(var %in% genes.to.test$gene.2) %>% add_tally() %>% arrange(p) %>% 
  #   mutate("i" = row_number(), "cv" = (i/n)*0.2, "p_cv" = p < cv)
  # # potential associations
  # features.to.test = burden.res.bulk.fdr %>% filter(p_cv == T) %>% tail(1) %>%
  #   {filter(burden.res.bulk.fdr, row_number() %in% 1:.$i)} %>% pull(feat) %>% unique()
  
  
  # 25k single cells sampled from all 7 plates
  data.bulk = readRDS("tmp/anycells.sampled.25ksub.rds")
  
  ###
  data.cell.count = data.sampled.bulk.gaussian %>% inner_join(cell.count.summ.well, by=c("Metadata_Plate","Metadata_Well")) %>%
    dplyr::select(Metadata_line_ID,Metadata_CellCount_avg)
  data.cell.count %>% ggplot(aes(factor(Metadata_line_ID), Metadata_CellCount_avg)) + geom_point(size=0.5) + 
    theme(axis.text.x=element_blank())
  ###
  
  # ---- temp -----
  # distribution of features
  data.bulk %>% dplyr::select(Metadata_line_ID,Cells_AreaShape_Compactness,Nuclei_AreaShape_Perimeter) %>%
    reshape2::melt("Metadata_line_ID") %>%
    ggplot(aes(factor(Metadata_line_ID),value)) + geom_boxplot(width=0.5, outlier.size=0.1) + 
    facet_wrap(~variable, ncol=1, scales="free_y") + geom_jitter(size=0.1) + guides(fill=F)
  
  # distribution of feature variance, %>% group_by(Metadata_line_ID) 
  vars.sc = data.bulk.sc %>% summarise_at(.vars=features.to.test, .funs=var) %>% #tibble::column_to_rownames("Metadata_line_ID") %>% 
    reshape2::melt()
  vars.well = data.bulk.well %>% summarise_at(.vars=features.to.test, .funs=var) %>% #tibble::column_to_rownames("Metadata_line_ID") %>% 
    reshape2::melt()
  vars.sc %>% inner_join(vars.well, by="variable", suffix = c(".sc",".well")) %>% filter(value.sc < 0.01 & value.well < 10) %>%
    ggplot(aes(value.sc, value.well)) + geom_point()
  
  vars %>% tibble::column_to_rownames("Metadata_line_ID") %>%
  reshape2::melt() %>% ggplot(aes(variable, value)) + geom_boxplot(width=0.2) + 
    theme(axis.text.x=element_blank())
  
  vars %>% ggplot(aes(Cells_AreaShape_Zernike_1_1)) + geom_density()
  # ---- temp -----
  
  # # beta estimate
  # sig.assoc.b = rbind(sig.assoc.b, burden.res.perm.emax %>% filter(feat == !!feature) %>% arrange(-est.max) %>%
  #                       group_by(feat) %>% add_tally() %>% head(floor(nrow(.)*0.05)) %>% tail(1) %>%
  #                       {filter(potent.assoc, feat == !!feature & abs(est) > .$est.max)})
  
  # # min p per feature per permutation
  # burden.res.perm.pmin = burden.res.perm %>% filter(feat %in% potent.assoc$feat & var %in% genes.to.test$gene.2) %>%
  #   group_by(feat,perm.seed) %>% mutate(p.min = min(p)) %>% dplyr::select(feat,perm.seed,p.min) %>% unique()
  # 
  # # max est per feature per permutation
  # burden.res.perm.emax = burden.res.perm %>% filter(feat %in% potent.assoc$feat & var %in% genes.to.test$gene.2) %>%
  #   group_by(feat,perm.seed) %>% mutate(est.max = max(abs(est))) %>% dplyr::select(feat,perm.seed,est.max) %>% unique()
  
  ##################
  
  
  cell.count.summ.neigh = data.bulk %>% dplyr::select(Metadata_line_ID,Cells_Neighbors_NumberOfNeighbors_Adjacent) %>% 
    group_by(Metadata_line_ID) %>% mutate("Cells_Neighbors_avg" = mean(Cells_Neighbors_NumberOfNeighbors_Adjacent)) %>% 
    dplyr::select(Metadata_line_ID,Cells_Neighbors_avg) %>% unique() %>%
    inner_join(cell.count.summ, by="Metadata_line_ID")
  
  feature.gene.var %>% inner_join(cell.count.summ, by="Metadata_line_ID") %>% dim()
  
  lines = cell.count.summ.neigh %>% arrange(Count_Cells_avg) %>% head(10) %>% pull(Metadata_line_ID)
  data.bulk %>% mutate("flag" = ifelse(Metadata_line_ID %in% lines,"low","high")) %>%
    ggplot(aes(Nuclei_AreaShape_Area, Nuclei_Granularity_1_DNA)) + geom_point(aes(color=flag))
  
  cor.test(cell.count.summ.neigh$Cells_Neighbors_avg, cell.count.summ.neigh$Count_Cells_avg)
  
  cell.count.summ.neigh %>% ggplot(aes(Count_Cells_avg,Cells_Neighbors_avg)) + geom_point() + theme_bw()
  ##################

  mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  attributes = c("ensembl_gene_id","start_position","end_position","hgnc_symbol","chromosome_name","gene_biotype")
  
  genes.pos = getBM(attributes=attributes, filters=c("hgnc_symbol"), values=list(hgnc_symbol=unique(grep("Meta|Sample|-|&",colnames(gene.vars), invert=T, value=T))), mart=mart, uniqueRows = T)
  
  dat = gene.vars %>% summarise_at(grep("Meta|Sample|-|&", colnames(.), value=T, invert=T), ~nnzero(.x)) %>% 
    colSums() %>% data.frame("n.vars" = .) %>% tibble::rownames_to_column("hgnc_symbol")
  
  genes.pos %>% mutate("len" = end_position-start_position) %>% filter(chromosome_name %in% 1:22 & gene_biotype %in% "protein_coding") %>% 
    inner_join(dat, by="hgnc_symbol") %>% filter(len < 250000 & n.vars < 100) %>%
    ggplot(aes(len, n.vars)) + geom_point(size = 0.5)
  
  # permutate
  res.perm = NULL
  # if(nrow(res %>% filter(p.adj < 0.05)) > 0){
  #   genes.to.test.perm = res %>% filter(p.adj < 0.05 & data == "orig") %>% dplyr::select("var") %>% mutate(gene.2 = sub("_", "-", sub("_inter","",var))) %>% dplyr::select("gene.2")
  #   print(paste("permuting for sig associated gene", nrow(genes.to.test.perm)))
  #   # permute interaction term or feature value
  #   for(i in 1:10000){
  #     if(reg.model == "baseline") feature.gene.var[, feature] = sample(feature.gene.var[, feature])
  #     if(reg.model == "interaction") feature.gene.var[, inter.var] = sample(feature.gene.var[, inter.var])
  #     res.perm = rbind(res.perm, do.call("rbind", apply(genes.to.test.perm, 1, function(x){ regress_feature_allele(feature = feature, gene = x, covars = covars, data = feature.gene.var, model = reg.model, level = level, inter.var = inter.var) })))
  #   }
  #   res.perm = res.perm %>% mutate(feat = feature, cells = round, variants = sum.or.present, 
  #                                  p.adj = p.adjust(p, "bon", n = nrow(genes.to.test)*length(features.to.test)), data = "perm")
  #   print("permutation done")
  # }
  
  # # ----- temp - save for cell health by greg -----
  # write.table(x = data.sampled.bulk.gaussian, file = paste0("bekar.cheeze/",round,".well.forCellHealth.tab"), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
  # # ----- temp - save for cell health by greg -----
  
  
  # 
  # frame.to.plot = res.perm %>% filter(feat == feature & var %in% c("GFM2","ANHX","CABYR")) %>% dplyr::select(feat,var,est,p) %>%
  #   mutate(p = format(p, digits=1, sc=T), variable = sub("_", "-", sub("_inter","",var)), assoc = ifelse(grepl("inter",var), "inter", "gene"))
  # # ---- temp ----
  
  # filter inter associations which are significant for baseline model
  # frame.to.plot = burden.res.bulk %>% filter(data == "orig" & p.adj < 0.05) %>% select(one_of("feat","var")) %>% dplyr::rename("var.gene"="var") %>%
  #   inner_join(burden.res.inter %>% filter(data == "orig" & p < 0.05) %>% mutate(variable = sub("_inter","",var)), by = c("feat","var.gene"="variable")) %>%
  #   dplyr::select(feat,var,est,p) %>% mutate(p = format(p, digits=1, sc=T), variable = sub("_", "-", sub("_inter","",var)), assoc = ifelse(grepl("inter",var), "inter", "gene"))
  
  
  d = read.table("tmp.tab", header=T)
  res = NULL
  for(maf in c(0.01, 0.1, 0.2)){
    res = rbind(res, data.frame(maf, "n"= d %>% filter(AF > maf & AF_Gn < 0.01) %>% nrow()))
  }
  d.rare = d %>% filter(AF < 0.01 | AF_Gn < 0.01)
  d.rare[sample(1:nrow(d.rare), 1000000), ] %>% ggplot(aes(AF,AF_Gn)) + geom_point(size = 0.5) + 
    geom_vline(data=res, aes(xintercept=maf), color = "blue") + geom_text(data=res, aes(x=maf+0.05, y=0.2, label = n)) + 
    geom_hline(yintercept=0.01, color = "blue") + theme_bw()
  
  
  burden.res.bulk %>% filter(p.adj < 0.1 & data == "orig") %>% dplyr::select(one_of("feat","p")) %>% dplyr::rename("p.orig" = "p") %>%
    inner_join(burden.res.perm.pmin, by = "feat") %>% group_by(feat) %>% mutate(n = length(which(p.min < p.orig))) %>%
    dplyr::select(feat,n) %>% unique()
  
  # # permutate
  # res.perm = NULL
  # for(i in 1:10){
  #   print(paste("permutation round",i))
  #   
  #   myop <- function(x){
  #     return(sample(x))
  #   }
  #   # 1. permute feature
  #   #feature.gene.var = feature.gene.var %>% mutate_at(.vars=c(feature), .funs=myop)
  #   # 2. permute genes
  #   feature.gene.var = feature.gene.var %>% mutate_at(.vars=c(genes.to.test$gene.2), .funs=myop)
  #   
  #   # burden test
  #   res.perm = rbind(res.perm, do.call("rbind", apply(genes.to.test, 1, function(x){ regress_feature_allele(feature = feature, gene = x, covars = covars, data = feature.gene.var, model = reg.model, level = level, inter.var = inter.var) })))
  # }
  # res.perm = res.perm %>% mutate(feat = feature, cells = round, variants = sum.or.present, 
  #                                p.adj = p.adjust(p, "bon", n = nrow(genes.to.test)*length(features.to.test)), data = "perm")
  # print("permutation done")
  
  
  # --------------
  # p value from permutation
  feats.assoc = burden.res.bulk %>% filter(p.adj < 0.1 & data == "orig") %>% dplyr::select(one_of("feat","var","p")) %>% dplyr::rename("p.orig" = "p") %>%
    inner_join(burden.res.perm %>% filter(data == "perm"), by = "feat") %>% group_by(feat) %>% mutate(n = length(which(p < p.orig))) %>%
    dplyr::select(feat,p.orig,n) %>% unique()
  
  
  ftrs = gene.vars %>% dplyr::select(one_of(grep("Meta|Sample",colnames(gene.vars), value=T, invert=T)))
  nzv = nearZeroVar(x = ftrs, saveMetrics = F, freqCut = 97/3, allowParallel = T, names = T)
  ftrs = ftrs %>% dplyr::select(-one_of(nzv))
  dim(ftrs)
  
  cv = findCorrelation(x = cor(ftrs), cutoff = 0.9, names = T, exact = F)
  length(cv)
  
  
  
  donor.meta %>% group_by(Metadata_Disease_Status,Metadata_Sex) %>% count() %>% pull(n) %>% matrix(.,2) %>% fisher.test()
  donor.meta %>% group_by(Metadata_Disease_Status,Metadata_Sex) %>% count()
  
  
  # ###
  # feats.assoc = burden.res.inter %>% filter(p.adj < 0.05 & data == "orig") %>% dplyr::select(one_of("feat","var","p")) %>% dplyr::rename("p.orig" = "p") %>%
  #   inner_join(burden.res.inter %>% filter(data == "perm"), by = c("feat","var")) %>% group_by(feat,var) %>% mutate(n = length(which(p < p.orig))) %>%
  #   dplyr::select(feat, var, n) %>% unique() %>% mutate(p=n/10000) %>% filter(p < 0.01)
  # 
  # xint.data = burden.res.inter %>% filter(p.adj < 0.05 & data == "orig" & feat %in% feats.assoc$feat & var %in% feats.assoc$var) %>% 
  #   dplyr::select(feat,var,p) %>% mutate(xint = -log10(p))
  # 
  # burden.res.inter %>% filter(data == "perm" & feat %in% feats.assoc$feat & var %in% feats.assoc$var) %>%
  #   ggplot(aes(-log10(p))) + geom_histogram(bins = 50) + facet_wrap(~feat+var, ncol = 4) + theme_bw() + 
  #   geom_vline(data = xint.data, color = "blue", aes(xintercept = xint))
  # ###
  
  # # feature to plot
  # feature = "Cytoplasm_Texture_InfoMeas1_DNA_20_02"
  # 
  # # gene asssociated to feature
  # pval.frame = res %>% arrange(p) %>% head(1) %>% dplyr::select(var,est,p) %>%
  #   mutate(p = format(p, digits=1, scientific=T), variable = sub("_", "-", sub("_inter","",var))) %>% mutate(assoc = ifelse(grepl("inter",var), "inter", "gene"))
  # 
  # # prepare data for plot
  # data = features.values %>% dplyr::select(c(Metadata_line_ID, Metadata_Disease_Category, !!as.name(inter.var), !!as.name(feature))) %>%
  #   inner_join(gene.vars %>% dplyr::select(-one_of(inter.var)), by = c("Metadata_line_ID"="Metadata_Project.Alias")) %>%
  #   filter(Metadata_line_ID %in% pcs.terra$Metadata_Project.Alias) %>%
  #   dplyr::select(c(!!as.name(feature), unique(pval.frame$variable), Metadata_line_ID, !!as.name(inter.var), Metadata_Disease_Category)) %>%
  #   reshape2::melt(id = c(feature,"Metadata_line_ID",inter.var,"Metadata_Disease_Category"))
  
  # # color for plot
  # disease.cols = union(my.cols[-c(4:3)],viridis::plasma(7))[-c(10:11)]
  # cond.cols = my.cols[c(4:3)]
  # 
  # # plot the feature according to variant in associated gene
  # if(nrow(pval.frame) == 1) wid = 3.5 else wid = 2.8*nrow(pval.frame)
  # pdf("plots/tmp.pdf", width = wid, height = 4.5)
  # p = ggplot(data, aes(x = factor(value), y = !!as.name(feature))) +
  #   geom_point(position = position_jitterdodge(), inherit.aes = T, alpha = 0.75, pch = 21, size = 2, stroke=0.5, aes(group = !!as.name(inter.var), fill = Metadata_Disease_Category))
  #   if(reg.model == "interaction") p = p + geom_boxplot(outlier.size = 0, alpha = 0.5, aes(color = !!as.name(inter.var)))
  #   if(reg.model == "baseline") p = p + geom_boxplot(outlier.size = 0, alpha = 0.5) 
  #   p + ggthemes::scale_fill_tableau(name = "disease") + scale_color_manual(values = cond.cols) +
  #   facet_wrap(~variable, scales = "free", nrow = 1) + theme_bw() +
  #   theme(panel.grid.major.x = element_blank(), strip.background = element_blank()) +
  #   theme(axis.title = element_text(size=12), axis.text = element_text(size=12), strip.text = element_text(size=12), title = element_text(size=12)) +
  #   xlab("number of variants") + ylab("gaussianised feature value") + ggtitle(feature, subtitle = round) +
  #   geom_text(data = pval.frame, mapping = aes(x = Inf, y = Inf, label = paste("beta =", round(est,2))), hjust = 1.1, vjust = 1.5, size = 4, parse = F) +
  #   geom_text(data = pval.frame, mapping = aes(x = Inf, y = Inf, label = paste("P =", p)), hjust = 1.1, vjust = 3.5, size = 4, parse = F) +
  #   geom_text(data = pval.frame, mapping = aes(x = Inf, y = Inf, label = assoc), hjust = 1.1, vjust = 5.5, size = 4, parse = F) +
  #   guides(fill = guide_legend(ncol = 1))
  # dev.off()
  # # geom_smooth(method = "lm", se=F, aes(color=Metadata_donor_Affected.Status), aes(group=1))
  
  # gene.vars = wgs.map.big.meta %>% mutate("wgs_id" = ifelse(test = grepl("CW20020",to_sample_id), "CW20020-P12-DH-12-4-17", toupper(gsub("\\_|\\.","-",collaborator_sample_id)))) %>% 
  #   dplyr::select(one_of("meta_id_final","wgs_id")) %>% unique() %>% 
  #   inner_join(donor.gt.re.sum %>% as.data.frame() %>% tibble::rownames_to_column("wgs_id") %>% mutate("wgs_id" = sub("^X","",toupper(gsub("\\_|\\.","-",wgs_id)))), by = "wgs_id") %>%
  #   inner_join(donor.meta %>% dplyr::select(c(Metadata_ID,Metadata_Project.Alias,!!as.name(inter.var))), by = c("meta_id_final"="Metadata_ID")) %>%
  #   dplyr::select(-one_of("meta_id_final","wgs_id"))
  # print("summarized variant count per gene per donor prepared")
  # print(dim(gene.vars))
  
  all.files = list.files("profiles", '.csv', full.names = T)
  profiles = all.files %>% purrr::map_df(read.csv)
  length(unique(profiles$Metadata_line_ID))
  
  sort(setdiff(profiles$Metadata_line_ID, donor.meta.wgs$Metadata_donor_Project.Alias))
  
  dim(donor.meta)
  write.csv(x = donor.meta, file = "metadata/metaData_donor_final_28july2020.csv", row.names = F)
  
  write.csv(x = donor.meta.wgs[, c("Metadata_donor_Project.Alias","Sample","Sample_ID")], 
            file = "metadata/metaData_donor_wgs.csv", row.names = F)
  

  # # genes to test for variant burden
  # genes.to.test = gene.vars %>% dplyr::select(-one_of(grep("Meta|Sample", colnames(gene.vars), value = T))) %>% 
  #   summarise_all(.funs = ~nnzero(.x), na.rm = T) %>% t() %>% data.frame("fr" = .) %>% tibble::rownames_to_column("gene") %>% 
  #   filter(fr > ceiling(nrow(gene.vars)*0.01) & fr < ceiling(nrow(gene.vars)*0.1), !grepl('\\.[0-9]|^RP[0-9]|^RPL|&|-[0-9]', gene)) %>% 
  #   mutate(gene.2 = sub("-", "_", gene)) %>% dplyr::select(gene.2)
  # print(paste("gene to test", nrow(genes.to.test)))
  
  
  ##################### feature vs disease status ########################
  if(F){
    print("inside feature vs disease")
    disease.to.feature = NULL
    for(feature in features.to.test){
      
      if(level == "bulk") covars = c("Metadata_donor_Affected.Status", "Metadata_Plate", "Metadata_image_batch", "Metadata_donor_Sex", "Metadata_donor_Age", paste0("PC",1:5))
      if(level == "well") covars = c("Metadata_donor_Affected.Status", "Metadata_Plate", "Metadata_image_batch", "Metadata_donor_Sex", "Metadata_donor_Age", paste0("PC",1:5), "Metadata_Well")
      if(round != "isolate") covars = union("Cells_Neighbors_NumberOfNeighbors_Adjacent", covars)
      
      feature.gene.var = features.values %>% dplyr::select(c(Metadata_line_ID, grep("PC[0-9]",covars,invert=T,value=T), !!as.name(feature))) %>%
        inner_join(pcs.terra[, c("Project.Alias", paste0("PC",1:5))], by = c("Metadata_line_ID" = "Project.Alias"))
      
      baseline = as.formula(paste(feature,"~", "Metadata_donor_Affected.Status", "+", paste(grep("Affected", covars, value=T, invert=T), collapse = "+")))
      lm.b = glm(formula = baseline, data = feature.gene.var)
      summ = coefficients(summary(lm.b))
      
      res = data.frame("var" = feature, "est" = summ["Metadata_donor_Affected.StatusYes", "Estimate"], "std.err" = summ["Metadata_donor_Affected.StatusYes", "Std. Error"], "p" = summ["Metadata_donor_Affected.StatusYes", "Pr(>|t|)"], row.names = NULL)
      disease.to.feature = rbind(disease.to.feature, res)
      print(paste("regression done for", feature))
      
    }
    
    # prepare data for plot
    pval.frame = disease.to.feature %>% mutate(p.adj = p.adjust(p, "bon", length(features.to.test))) %>% filter(p.adj < 0.05)
    pval.frame = disease.to.feature %>% filter(p < 0.05) %>% arrange(p) %>% head(1)
    feature = "Cytoplasm_Texture_AngularSecondMoment_RNA_20_01"
    
    features.values %>% dplyr::select(c(Metadata_line_ID, Metadata_donor_disease.name, Metadata_donor_Affected.Status, !!as.name(feature))) %>%
      ggplot(aes(Metadata_donor_Affected.Status, !!as.name(feature))) + geom_point(stroke=0.25, pch = 21, alpha = 0.5, size = 3, position = position_jitter(width = 0.25)) +
      geom_boxplot(alpha = 0.6, width = 0.6, outlier.shape = NA)
  }
  ########################################################################
  
  
  # for(findex in images.selected$FileIndex){
  #   images.selected.tiff = raster(paste0("images/cmqtl/", sub(":.*","",findex), "/", sub(".*:","",findex)))
  #   plot(images.selected.tiff, main = sub(":.*","",findex))
  # }
  
  
  data %>% filter(Metadata_donor_Affected.Status == "No" & value == 0) %>% arrange(!!as.name(feature)) %>%
    filter(row_number()==1 | row_number()==n()) %>% dplyr::select(Metadata_line_ID,!!feature,value) %>%
    {inner_join(x = data.sampled.bulk.gaussian, y = ., by = c("Metadata_line_ID",feature))} %>% 
    dplyr::select(Metadata_line_ID,!!feature,value,Metadata_Plate) %>%
    {inner_join(x = images.meta.onesite, y = ., by = c("Metadata_Plate","Metadata_line_ID"))} %>%
    mutate(FileIndex = paste(Metadata_Plate, sub("-.*","",FileName_OrigDNA), sep=":"))
  
  data %>% filter(Metadata_donor_Affected.Status == "No" & value == 0) %>% filter(!!as.name(feature) > -0.075 & !!as.name(feature) < 0.075) %>%
    dplyr::select(Metadata_line_ID,!!feature,value) %>%
    {inner_join(x = data.sampled.bulk.gaussian, y = ., by = c("Metadata_line_ID",feature))} %>% 
    dplyr::select(Metadata_line_ID,!!feature,value,Metadata_Plate) %>%
    {inner_join(x = images.meta.onesite, y = ., by = c("Metadata_Plate","Metadata_line_ID"))} %>%
    mutate(FileIndex = paste(Metadata_Plate, sub("-.*","",FileName_OrigDNA), sep=":"))
  
  data %>% filter(Metadata_donor_Affected.Status == "No" & value == 0) %>% arrange(!!as.name(feature)) %>%
    filter( (!!as.name(feature) > -0.01 & !!as.name(feature) < 0.01) | (row_number()==1 | row_number()==n()) ) %>%
    dplyr::select(Metadata_line_ID,!!feature,value) %>%
    {inner_join(x = data.sampled.bulk.gaussian, y = ., by = c("Metadata_line_ID",feature))} %>% 
    dplyr::select(Metadata_line_ID,!!feature,value,Metadata_Plate) %>%
    {inner_join(x = images.meta.onesite, y = ., by = c("Metadata_Plate","Metadata_line_ID"))} %>%
    mutate(FileIndex = paste(Metadata_Plate, sub("-.*","",FileName_OrigDNA), sep=":"))
  
  # ------------- Plots -------------
  if(F){
    # v3 height 5.5/7.5 x 4.5  
    pval.frame = res %>% filter(p.adj < 0.05) %>% arrange(p) %>% dplyr::select(var, p) %>%
      mutate(p = format(p, digits=1, scientific=T), variable = sub("_.*","",var)) %>% mutate(assoc = ifelse(grepl("inter",var), "inter", "gene"))
    disease.cols = union(my.cols[-c(4:3)],viridis::plasma(7))[-c(10:11)]
    cond.cols = my.cols[c(4:3)]
    
    feature.gene.var %>% dplyr::select(c(!!as.name(feature), unique(pval.frame$variable), Metadata_line_ID, Metadata_donor_Affected.Status)) %>%
      {inner_join(x = ., y = donor.meta[, c("Metadata_donor_Project.Alias","Metadata_donor_disease.name")], by = c("Metadata_line_ID"="Metadata_donor_Project.Alias"))} %>%
      reshape2::melt(id = c(feature,"Metadata_line_ID","Metadata_donor_Affected.Status","Metadata_donor_disease.name")) %>%
      mutate(xticks = paste(Metadata_donor_Affected.Status, value, sep = "_var")) %>%
      ggplot(aes(xticks, !!as.name(feature))) +
      geom_boxplot(width = 0.6, outlier.shape = NA, aes(color = Metadata_donor_Affected.Status)) +
      geom_point(stroke=0.1, pch = 21, alpha = 0.5, size = 3, position = position_jitter(width = 0.25), aes(fill = Metadata_donor_disease.name)) +
      facet_wrap(~variable, scales = "free") + theme_bw() +
      theme(panel.grid.major.x = element_blank(), strip.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      theme(axis.title = element_text(size=11), axis.text = element_text(size=12), strip.text = element_text(size=12), title = element_text(size=12)) +
      scale_fill_manual(values = disease.cols, name = "disease") + scale_color_manual(values = cond.cols, guide = F) +
      xlab("number of variants") + ylab("gaussianised feature value") + ggtitle(feature, subtitle = round) +
      geom_text(data = pval.frame, mapping = aes(x = Inf, y = Inf, label = paste("P =", p)), hjust = 1.1, vjust = 2, size = 4, parse = F)
    # geom_text(data = pval.frame, mapping = aes(x = Inf, y = Inf, label = assoc), hjust = 1.1, vjust = 3.5, size = 4, parse = F)
    
    # geom_smooth(method = "lm", se=F, aes(color=Metadata_donor_Affected.Status), aes(group=1))
    # feature.gene.var %>% filter(Metadata_donor_Affected.Status == "Yes" & BVES != 0) %>% pull(Metadata_line_ID) %>% 
    #   {filter(donor.meta, Metadata_donor_Project.Alias %in% .) %>% select(Metadata_donor_Project.Alias, Metadata_donor_disease.name)}
    
  }
  
  
  data.bulk = NULL
  for(file in list.files(path = paste0("rds.objects/",level), pattern = paste0(round,".*donor.rds"), full.names = T)){
    print(file)
    data.bulk.plate = readRDS(file)
    if(!is.null(data.bulk)){
      common.cols = intersect(colnames(data.bulk), colnames(data.bulk.plate))
      data.bulk = rbind(data.bulk[, common.cols], data.bulk.plate[, common.cols])
    } else data.bulk = data.bulk.plate
  }
  print("bulk data loaded from all plates")
  print(dim(data.bulk))
  
  
  
  #################################################
  # correlation per feature among replicates
  profile_files <- list.files("profiles", pattern = "[1-9|mt]_augmented.csv", full.names = T)
  profiles <- profile_files %>% purrr::map_df(readr::read_csv)
  replicate_correlation_values <- cytominer::replicate_correlation(profiles, names(profiles) %>% str_subset("Cells_|Cytoplasm_|Nuclei_"), 
                                                                   strata = "Metadata_line_ID", replicates = 8, split_by = "Metadata_Plate", cores = 4)
  
  replicate_correlation_values %>% readr::write_csv("bekar.cheeze/replicate_correlation_values.csv")
  
  # plot correlation for specific featues
  replicate_correlation_values = readr::read_csv("bekar.cheeze/replicate_correlation_values.csv")
  replicate_correlation_values %>% filter(str_detect(variable, "Cells_AreaShape_Zernike")) %>% 
    tidyr::separate("variable", c("x1", "x2", "x3", "n", "m")) %>% 
    ggplot(aes(n, m, size = median, label = sprintf("%.2f", median))) + geom_label()
  #################################################
  
  
  a = replicate_correlation_values %>% filter(grepl("AGP",variable)) %>% pull(median)
  b = replicate_correlation_values %>% filter(grepl("Edge",variable)) %>% pull(median)
  summary(b)
  
  # # get top genes to plot
  # gene.to.plot = res %>% filter(p.adj < 0.05 & !grepl("disease|inter",var)) %>% arrange(p.adj) %>% pull(var) %>% unique()
  # # get top interactions to plot
  # gene.to.plot = res %>% filter(p < 0.0001 & grepl("inter",var)) %>% mutate(gene = sub("_.*","",var)) %>%
  #   pull(gene) %>% unique()
  
  # pval.frame = res %>% filter(p.adj < 0.05 & grepl("inter",var)) %>% dplyr::select(var, p) %>% #rename(c(variable=var)) %>%
  #   mutate(p = format(p, digits=1, scientific=T), variable = sub("_.*","",var)) %>% mutate(assoc = ifelse(grepl("inter",var), "inter", "gene"))
  
  # ############################
  # # gene.to.plot = res %>% filter(p.adj < 0.05) %>% arrange(p.adj) #%>% mutate(gene = sub("_.*","",var)) %>% pull(gene) %>% unique()
  
  # 
  # # qqplot for all tested genes 4x3.5
  # if(reg.model == "interaction"){
  #   ps = burden.res %>% filter(grepl("inter",var) & data == "orig") %>% pull(model.imp)
  # } else if(reg.model == "baseline"){
  #   ps = burden.res %>% filter(!grepl("inter",var) & data == "orig") %>% pull(p)
  # }
  # qqman::qq(pvector = ps, main = paste0(feature, paste("\nlambda = ",round(GenABEL::estlambda(ps, method="regression")$estimate,2))))
  # ############################
  
  # pval.frame = res %>% filter(var %in% gene.to.plot) %>% dplyr::select(var, p) %>% #rename(c(variable=var)) %>%
  #   mutate(p = format(p, digits=1, scientific=T), variable = sub("_.*","",var))
  
  # # v1 feature value vs. variant counts in top genes 3/5 x 4
  # feature.gene.var %>% dplyr::select(c(!!as.name(feature), gene.to.plot, "Metadata_donor_Affected.Status")) %>%
  #   reshape2::melt(id = c(feature, "Metadata_donor_Affected.Status")) %>%
  #   ggplot(aes(factor(value), !!as.name(feature))) + geom_boxplot(width = 0.4, outlier.size = 0, aes(color = factor(value))) +
  #   geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, aes(color = Metadata_donor_Affected.Status)) + theme_bw() + facet_wrap(~variable, scales = "free") +
  #   theme(panel.grid.major.x = element_blank(), strip.background = element_blank()) +
  #   theme(axis.title = element_text(size=13), axis.text = element_text(size=13), strip.text = element_text(size=13), title = element_text(size=13)) +
  #   scale_color_manual(values = c(my.cols[1:2],my.cols[4:3])) + xlab("number of variants") + ggtitle(round) +
  #   geom_text(data = pval.frame, mapping = aes(x = Inf, y = Inf, label = paste("P =", p)), hjust = 1.05, vjust = 2, size = 4, parse = F) +
  #   stat_summary(fun=median, geom="line", color="blue", aes(group=1))
  
  # # v2 height 3.5
  # feature.gene.var %>% dplyr::select(c(!!as.name(feature), pval.frame$variable, "Metadata_donor_Affected.Status")) %>%
  #   reshape2::melt(id = c(feature, "Metadata_donor_Affected.Status")) %>%
  #   ggplot(aes(factor(value), !!as.name(feature))) + 
  #   geom_boxplot(width = 0.4, outlier.size = 0, aes(color = Metadata_donor_Affected.Status)) +
  #   geom_point(pch = 21, alpha = 0.5, size = 1.5, position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2), aes(color = Metadata_donor_Affected.Status)) +
  #   facet_wrap(~variable, scales = "free") + theme_bw() +
  #   theme(panel.grid.major.x = element_blank(), strip.background = element_blank()) +
  #   theme(axis.title = element_text(size=11), axis.text = element_text(size=11), strip.text = element_text(size=11), title = element_text(size=11)) +
  #   scale_color_manual(values = c(my.cols[4:3]), guide = F) + xlab("number of variants") + ylab("gaussianised feature value") + ggtitle(feature, subtitle = round) +
  #   geom_text(data = pval.frame, mapping = aes(x = Inf, y = Inf, label = paste("P =", p)), hjust = 1.1, vjust = 2, size = 4, parse = F) +
  #   geom_text(data = pval.frame, mapping = aes(x = Inf, y = Inf, label = assoc), hjust = 1.1, vjust = 3.5, size = 4, parse = F)
  
  
  
  
  # load bulk gaussianized data and get features to test
  if(T){
    # combine summarized from plates
    data.bulk = NULL
    for(file in list.files(path = "rds.objects/bulk", pattern = paste0(round,".*donor.rds"), full.names = T)){
      print(file)
      data.bulk.plate = readRDS(file)
      if(!is.null(data.bulk)){
        common.cols = intersect(colnames(data.bulk), colnames(data.bulk.plate))
        data.bulk = rbind(data.bulk[, common.cols], data.bulk.plate[, common.cols])
      } else data.bulk = data.bulk.plate
    }
    print("bulk data loaded from all plates")
    print(dim(data.bulk))
    
    # gaussianize each bulk feature across all plates
    myop <- function(x){
      return( RNOmni::rankNorm(u = x))
    }
    data.sampled.bulk.gaussian = data.bulk %>% mutate_at(.vars = grep("Meta", colnames(data.bulk), value = T, invert = T), .funs = myop)
    print("features at bulk level gaussianized")
    print(dim(data.sampled.bulk.gaussian))
    
    # remove certain columns
    cols.to.remove = grep("Sample_|_Children_|_Number_|_Parent_|_Center_|_Location_|_Count_|Granularity_1[4-6]", colnames(data.sampled.bulk.gaussian), value = T)
    data.sampled.bulk.gaussian = data.sampled.bulk.gaussian %>% dplyr::select(-one_of(cols.to.remove))
    print("certain columns removed")
    print(dim(data.sampled.bulk.gaussian))
    print(head(data.sampled.bulk.gaussian$Cells_Neighbors_NumberOfNeighbors_Adjacent))
    
    # remove features with near zero variance
    ftrs = data.sampled.bulk.gaussian[, grep("Meta", colnames(data.sampled.bulk.gaussian), value = T, invert = T)]
    nzv = nearZeroVar(x = ftrs, saveMetrics = F, freqCut = 90/10, allowParallel = T, names = T)
    data.sampled.bulk.gaussian = data.sampled.bulk.gaussian %>% dplyr::select(-one_of(nzv))
    print("near zero variance features removed")
    print(dim(data.sampled.bulk.gaussian))
    
    # take uncorrelated features for burden testing
    ftrs = data.sampled.bulk.gaussian[, grep("Meta", colnames(data.sampled.bulk.gaussian), value = T, invert = T)]
    cv = findCorrelation(x = cor(ftrs), cutoff = 0.9, names = T, exact = F)
    features.uncor = data.frame(feature = setdiff(colnames(ftrs), cv))
    print(paste("uncorrelated features taken for burden testing", nrow(features.uncor)))
    
    # temp #
    feat.cat = unique(sub(".*_(.*?)_.*", "\\1", grep("Meta",colnames(data.sampled.bulk.gaussian),value = T,invert = T)))
    feat.cat.cols = my.cols[1:length(feat.cat)]
    names(feat.cat.cols) = feat.cat
    
    # combine donor information
    data.sampled.bulk.gaussian = merge(x = data.sampled.bulk.gaussian, y = donor.meta, by.x = "Metadata_line_ID", by.y = "Metadata_donor_Project.Alias", sort = F)
    print("donor info added to data with all features")
    print(dim(data.sampled.bulk.gaussian))
    
    # add image batch information
    data.sampled.bulk.gaussian = data.sampled.bulk.gaussian %>% mutate(Metadata_image_batch = ifelse(Metadata_Plate %in% c("BR00106708","BR00106709"), "a", ifelse(Metadata_Plate %in% c("BR00107339","BR00107338"), "b", "c")) )
    print("metadata image batch added")
    print(dim(data.sampled.bulk.gaussian))
    
    # remove cell lines which are on >1 plates
    data.sampled.bulk.gaussian = data.sampled.bulk.gaussian[!duplicated(data.sampled.bulk.gaussian$Metadata_line_ID), ]
    print("duplicated lines removed")
    print(dim(data.sampled.bulk.gaussian))
    
  }
  
  
  
  # # load bulk level gaussianied data # VERIFY #
  # d = read.delim(paste0("rds.objects/bulk/INTFeatures_Indlevel_mean_SingleCell.txt"), sep = " ")
  # print(paste("bulk gaussianised data loaded for", round))
  # print(dim(data.sampled.bulk.gaussian))
  
  
  
  # # genes to test for variant burden
  # # 1. based on total variant count across donor
  # # genes.to.test = feature.gene.var %>% select(-c(!!as.name(feature), covars)) %>% colSums() %>% data.frame("sum" = .) %>% 
  # #   tibble::rownames_to_column("gene") %>% filter(sum > 0, !grepl('\\.[0-9]|^RP[0-9]|^RPL|&', gene)) %>% select(gene)
  # 
  # # 2. based on how many donor have variants
  # genes.to.test = feature.gene.var %>% select(-c(!!as.name(feature), covars)) %>% summarise_all(.funs = ~nnzero(.x), na.rm = T) %>% 
  #   t() %>% data.frame("fr" = .) %>% tibble::rownames_to_column("gene") %>% 
  #   filter(fr > nrow(feature.gene.var)*0.01, !grepl('\\.[0-9]|^RP[0-9]|^RPL|&', gene)) %>% select(gene)
  # print(paste("gene to test", nrow(genes.to.test)))
  
  # load bulk level gaussianied data # VERIFY #
  round = "Isolate"
  data.sampled.bulk.gaussian = read.delim(paste0("rds.objects/bulk/INTFeatures_Indlevel_mean_", round, ".txt"), sep = " ")
  print(paste("bulk gaussianised data loaded for", round))
  print(dim(data.sampled.bulk.gaussian))
  
  # take uncorrelated features
  res = NULL
  for(cor.cut in seq(0.5,0.9,0.1)){
    ftrs = data.sampled.bulk.gaussian[, grep("Meta", colnames(data.sampled.bulk.gaussian), value = T, invert = T)]
    cv = findCorrelation(x = cor(ftrs), cutoff = cor.cut, names = T, exact = F)
    res = rbind(res, data.frame("cor" = cor.cut, "nfeat" = nrow(data.frame(feature = setdiff(colnames(ftrs), cv)))))
  }
  res = rbind(res, data.frame("cor" = 1, "nfeat" = ncol(data.sampled.bulk.gaussian)))
  ggplot(res, aes(factor(cor), nfeat)) + geom_bar(stat = "identity", width = 0.35, fill = "blue") + theme_bw() + theme(panel.grid.major.x = element_blank()) +
    geom_text(aes(label=nfeat), vjust=-0.5) + xlab("inter-feature correlation") + ylab("number of features") + ggtitle(round) +
    scale_y_continuous(expand = expansion(mult = c(0,.1)))
  

  
  # ####################################
  # setwd("~/Documents/cmQTL/codes/")
  #system('python ~/Documents/cmQTL/codes/temp2.py cmqtlpl1.5-31-2019-mt:r04c09f05p01', wait=F)
  use_python("/usr/local/bin/python3")
  #reticulate::source_python("~/Documents/cmQTL/codes/temp2.py")
  reticulate::py_run_file("~/Documents/cmQTL/codes/temp2.py")
  # data.sampled = readRDS("r2.more.cleaned.rds")
  # print("sampled processed data loaded")
  # print(dim(data.sampled))
  # setwd("/Volumes/ja286/cmQTL")
  # ####################################

    
  # load single cell data
  # data.sampled = readRDS(paste0("rds.objects/singleCell/sampled.data/cleaned.data/", round, ".more.cleaned.rds"))
  # print(paste("cleaned raw data from wallace loaded for", round))
  # print(dim(data.sampled))
  
  # convert single to donor level
  if(F){
    print(paste("summarizing data at donor level for", round))
    
    # make bulk data
    data.sampled.bulk = as.data.frame(data.sampled %>% group_by(Metadata_Plate, Metadata_line_ID) %>% summarise_at(grep("Meta", colnames(data.sampled), invert = T, value = T), mean))
    print(dim(data.sampled.bulk))
    print(paste("number of cell lines in bulk is", length(unique(data.sampled.bulk$Metadata_line_ID))))
    
    # gaussianize each bulk feature across all plates
    data.gaussian = apply(X = data.sampled.bulk[, grep("Meta", colnames(data.sampled.bulk), value = T, invert = T)], MARGIN = 2, FUN = function(x) {RNOmni::rankNorm(u = x)})
    data.sampled.bulk.gaussian = cbind(data.sampled.bulk[, grep("Meta", colnames(data.sampled.bulk), value = T)], data.gaussian)
    print("features at bulk level gaussised")
    print(dim(data.sampled.bulk.gaussian))
    
    # combine donor information
    data.sampled.bulk.gaussian = merge(x = data.sampled.bulk.gaussian, y = donor.meta, by.x = "Metadata_line_ID", by.y = "Metadata_donor_Project.Alias", sort = F)
    print("donor info added to data with all features")
    print(dim(data.sampled.bulk.gaussian))
    
    # add image batch information
    data.sampled.bulk.gaussian = data.sampled.bulk.gaussian %>% mutate(Metadata_image_batch = ifelse(Metadata_Plate %in% c("BR00106708","BR00106709"), "a", ifelse(Metadata_Plate %in% c("BR00107339","BR00107338"), "b", "c")) )
    print("metadata image batch added")
    print(dim(data.sampled.bulk.gaussian))
    
    # remove cell lines which are on >1 plates
    data.sampled.bulk.gaussian = data.sampled.bulk.gaussian[!duplicated(data.sampled.bulk.gaussian$Metadata_line_ID), ]
    print("duplicated lines removed")
    print(dim(data.sampled.bulk.gaussian))
    
    # on plate edge or not
    # data.model$Metadata_onEdge = sub("[0-9]+", "", data.model$Metadata_Well) %in% c("A","B","O","P") + sub("[A-Z]+", "", data.model$Metadata_Well) %in% c("01","24")
    
  }
  
  # files = list.files(path = "extreme.values/burden.out", pattern = "CMC.assoc", full.names = T)
  # burden.res = NULL
  # for(file in files){
  #   print(sub(".*\\/(.*).CMC.assoc", "\\1", file))
  #   tmp = read.delim(file) %>% select(Gene, Pvalue) %>% mutate(feature = sub(".*\\/(.*).CMC.assoc", "\\1", file))
  #   burden.res = rbind(burden.res, tmp)
  # }
  # print(dim(burden.res))
  # burden.res %>% group_by(feature) %>% mutate(Pvalue.adj = p.adjust(p = Pvalue, method = "bon")) %>% filter(Pvalue.adj < 0.05) %>% mutate()
  #   arrange(Pvalue.adj) %>%
  #   ggplot(aes())
  # 
  # saveRDS(object = burden.res, file = "burden1.rds")
  

d = read.delim("wgs/snpsif.rare.tab")
d %>% filter(ANN.0..EFFECT != "" | ANN.0..IMPACT != "") %>% .$ANN.0..IMPACT %>% table()
  

seq = NULL
for(id in donor.meta$Metadata_donor_ID){
  hits = grep(id, as.character(wgs$V1), value = T)
  if(length(hits) > 0) seq = rbind(seq, data.frame(id, hits))
}

commons = intersect(donor.meta$Metadata_donor_BSP.Collab.Sample.ID, wgs$V1)

pat = paste(sub("SCBB-", "", setdiff(donor.meta$Metadata_donor_BSP.Collab.Sample.ID, wgs$V1)), collapse = "|")
pat = paste(donor.meta[!donor.meta$Metadata_donor_BSP.Collab.Sample.ID %in% commons, "Metadata_donor_BSP.Collab.Participant.ID"], collapse = "|")
pat = paste(donor.meta[!donor.meta$Metadata_donor_BSP.Collab.Sample.ID %in% commons, "Metadata_donor_BSP.Collab.Sample.ID"], collapse = "|")

commons = union(commons, grep(pat, wgs$V1, value = T))

length(intersect(wgs$V1, commons))

length(intersect(donor.meta$Metadata_donor_BSP.Collab.Participant.ID, wgs$V1))


id = donor.meta$Metadata_donor_BSP.Collab.Sample.ID[1]


seq = NULL
for(id in donor.meta$Metadata_donor_BSP.Collab.Sample.ID){
  a = sub("SCBB-([0-9]+).*", "\\1", id)
  pat = sub("_.*", "", a)
  hits = grep(pat, as.character(wgs$V1), value = T)
  if(length(hits) == 1) seq = rbind(seq, data.frame(id, hits)) else if(length(hits) >1) seq = rbind(seq, data.frame(id, "hits" = hits[1]))
}
dim(seq)
write.table(x = seq$hits, file = "wgs.good.samples.txt", append = F, quote = F, sep = "\n", row.names = F, col.names = F)



# # 1. from terra data sample section
# samples = read.delim("bekar.cheeze/sample_july2020.tsv")
# rownames(samples) = sub("_.*", "", sub("SCBB[-_]([0-9]+)[-_]*.*", "\\1", samples$collaborator_sample_id))

# 2. from 1st row of vcf file
#wgs = read.table("bekar.cheeze/WGS_AllSamples.txt") %>% mutate("id" = sub("_.*", "", sub("SCBB[-_]([0-9]+)[-_]*.*", "\\1", V1)))
wgs = 
length(intersect(wgs$id, donor.meta$Metadata_donor_BSP.Collab.Sample.ID))

# 3. donor meta data
rownames(donor.meta) = sub("_.*", "", sub("SCBB[-_]([0-9]+)[-_]*.*", "\\1", donor.meta$Metadata_donor_BSP.Collab.Sample.ID))

#merge(x = donor.meta[, "Metadata_donor_BSP.Collab.Sample.ID", drop=F], y = samples[, "collaborator_sample_id", drop=F], by = 0)
good.samples = intersect(rownames(samples), rownames(donor.meta))
length(good.samples)
write.table(x = intersect(samples[good.samples, "collaborator_sample_id"], wgs$V1), file = "wgs.good.samples.txt", append = F, quote = F, sep = "\n", row.names = F, col.names = F)

setdiff(rownames(samples), rownames(donor.meta))

# donor.meta[setdiff(rownames(donor.meta), rownames(samples)), ]
# nowgs.confirmed = c("CW20020", "CW60331", "CW60365", "CW60040", "CW60289")
length(intersect(samples$collaborator_sample_id, wgs$V1))

#CW60094FF1_P13_MT_4_6_18 CW60094FF1

# note 1. collaborator_sample_id is name of sample in wgs file
# note 2. collaborator_participant_id is name of sample in metadata file
wgs = read.table("bekar.cheeze/tmp.txt") %>% mutate("id" = sub("_.*", "", sub("SCBB[-_]([0-9]+)[-_]*.*", "\\1", V1)))
dm = donor.meta %>% mutate("wgs.map.id" = gsub("\\.","_",Metadata_donor_BSP.Collab.Sample.ID))
setdiff(wgs$V1, dm$wgs.map.id)

head(donor.meta)
donor.meta %>% filter(grepl("1749",Sample_ID))
donor.meta %>% filter(grepl("1749",Metadata_donor_ID))
donor.meta %>% filter(grepl("1749",Metadata_donor_BSP.Collab.Participant.ID))
donor.meta %>% filter(grepl("1749",Metadata_donor_BSP.Collab.Sample.ID))

}
#------------------------------------------------------------------------------------------------------------------
