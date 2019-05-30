## cmQTL ## 

"Hypothesis: colony-specific-effects and isolated-specific-effects might be driving the patterns
we see in the optimal plate density/time analysis. 

The point of this analysis is to see if our notion of 'isolated' cell versus 'colony' cell
has any ability to distinguish the features that characterize each of the cells. 
if these two categories do not distinguish cells in the feature space, maybe 
there is no colony-specific-effect versus isolated-specific-effect that will confound our 
decisions regarding optimal plating time and density"

#notes on data: 
#BR00103267 is 6hrs, with densities: 5000 15000 10000 20000
#BR00103268 is 24hrs, with densities: 1000 2000 3000 4000 
#took colony cells from each of the wells and meta analyzed features
#took isolated cells from each of the wells and meta analyzed features

data_files <- c("BR00103267","BR00103268")
timepoints <- c("6hrs","24hrs")
library(umap)
for (i in 1:length(data_files)){
  dat_colony <- read.table("data/",paste0(data_files[i],"_colony_normalized_variable_selected.csv"), header = T, stringsAsFactors = F, sep = ",")
  dat_isolated <- read.table("data/",paste0(data_files[i],"_isolated_normalized_variable_selected.csv"), header = T, stringsAsFactors = F, sep = ",")
  dim(dat_colony) 
  dim(dat_isolated) 
  
  table(dat_isolated$Metadata_line_ID) #cell line
  d <- data.frame(V1 = rbind(dat_isolated[,9:ncol(dat_isolated)], dat_colony[,9:ncol(dat_colony)])) #include only cell morphology features (not ID columns)
  label_vec <- c(rep(rgb(0.2,0.5,0.2,0.5), nrow(dat_isolated)), rep(rgb(0.6,0.2,0.2,0.5),nrow(dat_colony))) #add ID for isolated versus colony
  
  dat_prc <- prcomp(na.omit(t(d)),center = TRUE, scale = TRUE) #NaN values cause error w/o na.omit()
  dim(dat_prc$rotation)
  dat_prc_loadings <- dat_prc$rotation
  
  #umap on 30 PCs (where elbow of sd plot is)
  dat_prc_umap <- umap(dat_prc_loadings[,1:30])
  
  pdf(paste0("UMAP_",timepoints[i],"_perDens_IsolatedvColony_acrossCellLines.pdf"), width = 10, height = 7)
  par(mfrow = c(2,4))
  #color by single / colony 
  densities <- c(dat_isolated$Metadata_plating_density, dat_colony$Metadata_plating_density)
  densities_uniq <- sort(unique(densities),decreasing = F)
  for (i in 1:length(densities_uniq)){
    w_plot <- which(densities == densities_uniq[i])
    plot(dat_prc_umap$layout[w_plot,1], dat_prc_umap$layout[w_plot,2], xlab= "UMAP 1", ylab= "UMAP 2", col = label_vec[w_plot], pch = 19, main=paste0("density=",densities_uniq[i]))
    legend("topright", col = c(rgb(0.2,0.5,0.2,0.5),rgb(0.6,0.2,0.2,0.5)), legend = c("isolated","colony"), bty = "n", pch = 19)
  }
  
  #color by cell type 
  celllines <- c(dat_isolated$Metadata_line_ID, dat_colony$Metadata_line_ID)
  celllines_unique <- sort(unique(celllines),decreasing = F)
  mycols <- match(celllines, celllines_unique)
  mycols_key <- c(rgb(0.2,0.5,0.2,0.5), rgb(0.5,0.2,0.2,0.5), rgb(0.2,0.2,0.5,0.5), rgb(0.6,0.3,0.6,0.5), rgb(0.3,0.6,0.6,0.5), rgb(0.6,0.6,0.3,0.5))
  mycols <- mycols_key[mycols]
  for (i in 1:length(densities_uniq)){
    w_plot <- which(densities == densities_uniq[i])
    plot(dat_prc_umap$layout[w_plot,1], dat_prc_umap$layout[w_plot,2], xlab= "UMAP 1", ylab= "UMAP 2", col = mycols[w_plot], pch = 19, main=paste0("density=",densities_uniq[i]))
    legend("topright", col = mycols_key, legend = c("A","B","C","D","E","F"), bty = "n", pch = 19)
  }
  dev.off()
}

##########################################################

#### investigate super tight cell cluster #### 

#what characterizes the super tight cluster in 24hrs
w <- which(dat_isolated_prc_umap$layout[w_plot,2] > 6)
s <- sapply(1:nrow(d), function(x) length(which(abs(d[x,]) == "Inf")))
boxplot(s[w],s[-w]) #more Inf values in tight cluster than not. 
wilcox.test(s[w],s[-w]) #p = 0.045

s <- sapply(1:nrow(d), function(x) length(which(abs(d[x,]) == "NaN")))
boxplot(s[w],s[-w]) #comparable number of NAs, but tight cluster has more (not significant) 
wilcox.test(s[w],s[-w]) #p = 0.15

rs_w <- rowSums(d[w,], na.rm = T)
rv_w <- rowVars(d[w,], na.rm = T)
rs_notw <- rowSums(d[-w,], na.rm = T)
rv_notw <- rowVars(d[-w,], na.rm = T)
boxplot(rs_w,rs_notw, na.rm = T)
t.test(rs_w[abs(rs_w)!="Inf"],rs_notw[abs(rs_notw)!="Inf"]) #not significant
boxplot(rv_w[abs(rv_w)!="Inf"],rv_notw[abs(rv_notw)!="Inf"])
t.test(rv_w[abs(rv_w)!="Inf"],rv_notw[abs(rv_notw)!="Inf"]) #not significant

##########################################################

########### combine time point / density / cell type into one UMAP ##########

dat_colony_6hrs <- read.table("data/BR00103267_colony_normalized_variable_selected.csv", header = T, stringsAsFactors = F, sep = ",")
dat_isolated_6hrs <- read.table("data/BR00103267_isolated_normalized_variable_selected.csv", header = T, stringsAsFactors = F, sep = ",")
dat_colony_24hrs <- read.table("data/BR00103268_colony_normalized_variable_selected.csv", header = T, stringsAsFactors = F, sep = ",")
dat_isolated_24hrs <- read.table("data/BR00103268_isolated_normalized_variable_selected.csv", header = T, stringsAsFactors = F, sep = ",")

dat_all <- rbind(dat_colony_6hrs,dat_isolated_6hrs,dat_colony_24hrs,dat_isolated_24hrs)
dat_prc <- prcomp(na.omit(t(dat_all[,9:ncol(dat_all)])),center = TRUE, scale = TRUE)
dim(dat_prc$rotation) #1460 x 584
dat_prc_loadings <- dat_prc$rotation
dat_prc_umap <- umap(dat_prc_loadings[,1:30])

pdf("UMAP_alldata.pdf", width = 10, height = 10)
par(mfrow = c(2,2))

#color by isolated/colony
label_vec <- c(rep(rgb(0.2,0.5,0.2,0.5), nrow(dat_colony_6hrs)), rep(rgb(0.6,0.2,0.2,0.5),nrow(dat_isolated_6hrs)),rep(rgb(0.2,0.5,0.2,0.5), nrow(dat_colony_24hrs)), rep(rgb(0.6,0.2,0.2,0.5),nrow(dat_isolated_24hrs)))
plot(dat_prc_umap$layout[,1], dat_prc_umap$layout[,2], xlab= "UMAP 1", ylab= "UMAP 2", col = label_vec, pch = 19, main = "Isolated/Colony")
legend("topright", col = c(rgb(0.2,0.5,0.2,0.5),rgb(0.6,0.2,0.2,0.5)), legend = c("colony","isolated"), bty = "n", pch = 19)

#color by time point
label_vec <- c(rep(rgb(0.2,0.5,0.2,0.5), nrow(dat_colony_6hrs) + nrow(dat_isolated_6hrs)), rep(rgb(0.6,0.2,0.2,0.5),nrow(dat_isolated_24hrs) + nrow(dat_colony_24hrs)))
plot(dat_prc_umap$layout[,1], dat_prc_umap$layout[,2], xlab= "UMAP 1", ylab= "UMAP 2", col = label_vec, pch = 19, main = "Time Point")
legend("topright", col = c(rgb(0.2,0.5,0.2,0.5),rgb(0.6,0.2,0.2,0.5)), legend = c("6 hours","24 hours"), bty = "n", pch = 19)

#color by density 
densities <- dat_all$Metadata_plating_density
densities_uniq <- sort(unique(densities),decreasing = F)
m <- match(densities,densities_uniq)
mycols_key <- c(rgb(0.2,0.5,0.2,0.5), rgb(0.5,0.2,0.2,0.5), rgb(0.2,0.2,0.5,0.5), rgb(0.6,0.3,0.6,0.5),rgb(0.3,0.6,0.6,0.5), rgb(0.6,0.6,0.3,0.5), rgb(0.5,0.5,0.5,0.5), rgb(0.1,0.3,0.5,0.5))
mycols <- mycols_key[m]
plot(dat_prc_umap$layout[,1], dat_prc_umap$layout[,2], xlab= "UMAP 1", ylab= "UMAP 2", col = mycols, pch = 19, main = "Densities")
legend("bottomleft", col = mycols_key, legend = densities_uniq, bty = "n", pch = 19)

#color by cell line 
celllines <- dat_all$Metadata_line_ID
celllines_unique <- sort(unique(celllines),decreasing = F)
mycols <- match(celllines, celllines_unique)
mycols_key <- c(rgb(0.2,0.5,0.2,0.5), rgb(0.5,0.2,0.2,0.5), rgb(0.2,0.2,0.5,0.5), rgb(0.6,0.3,0.6,0.5), rgb(0.3,0.6,0.6,0.5), rgb(0.6,0.6,0.3,0.5))
mycols <- mycols_key[mycols]
plot(dat_prc_umap$layout[,1], dat_prc_umap$layout[,2], xlab= "UMAP 1", ylab= "UMAP 2", col = mycols, pch = 19, main = "Cell Line")
legend("bottomleft", col = mycols_key, legend = c("A","B","C","D","E","F"), bty = "n", pch = 19)

dev.off()

