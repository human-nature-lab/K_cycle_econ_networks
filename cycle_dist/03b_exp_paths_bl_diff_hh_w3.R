library(tidyverse)
library(intergraph)
library(igraph)

#-------------------------- Load data --------------------------------
setwd("/scratch")

# Create a data path for saving cleaned / augmented data sets
data_path <- "/scratch"

# Load node / edge list for borrowing ties
load(paste(c(data_path, "/borrow_lend_diff_hh_undir_edge_node_list_w1w3_*.rda"), collapse=""))

# Load generated network attributes to access kcycle prevalence per village
load(paste(c(data_path, "/net_attr_borrow_lend_undir_ds2_w1w3_*.rda"), collapse=""))

#------------------------------____--------------------------

#--------------------- Path density ------------------------
# Calculate the path density using the degree distribution of the observed graphs
# (individual borrow and lend ties) to generate random graphs for estimating the expected cycles.

##Do 50 random villages per village
all_vill_w3_vec <- c(1:176)
vill_w3_vec <- all_vill_w3_vec[!all_vill_w3_vec %in% c(155, 156)]

#vil_all_bl_diff_hh_kp2<-array(NA,dim=c(length(vill_w3_vec),50))
vil_all_bl_diff_hh_kp3<-array(NA,dim=c(length(vill_w3_vec),50))
vil_all_bl_diff_hh_kp4<-array(NA,dim=c(length(vill_w3_vec),50))
vil_all_bl_diff_hh_kp5<-array(NA,dim=c(length(vill_w3_vec),50))

library(sna, quietly = TRUE)

for(j in 1:length(vill_w3_vec)){
  vill_code <- vill_w3_vec[j]
  edgelist_subset <- subset(borrow_lend_diff_hh_undir_edgelist_w3, village_code_w3==vill_code) # & relationship=="trust_borrow_money"
  com <- graph_from_data_frame(edgelist_subset, directed=FALSE)
  com <- simplify(com)
  outd <- igraph::degree(com, mode = "out")
  ind <- igraph::degree(com, mode = "in")
  
  for(kk in 1:50) {
    sim <- sample_degseq(outd, ind, method = "simple.no.multiple")
    sim <- asNetwork(sim)
    # sim_adj = get.adjacency(sim,sparse=FALSE)
    
    kpath_count=kpath.census(sim, maxlen = 5, mode = "graph", tabulate.by.vertex = TRUE,
                             path.comembership = "none", dyadic.tabulation = "none")
    # detach("package:sna", unload = TRUE)
    
    #vil_all_bl_diff_hh_kp2[j,kk] <- kpath_count$path.count[2,1]
    vil_all_bl_diff_hh_kp3[j,kk] <- kpath_count$path.count[3,1]
    vil_all_bl_diff_hh_kp4[j,kk] <- kpath_count$path.count[4,1]
    vil_all_bl_diff_hh_kp5[j,kk] <- kpath_count$path.count[5,1]
  }
  print(j)
}

### Path density
# Data frame of the number of paths per village
vil_all <- kstats_borrow_lend_diff_hh_undir_w3_df %>%
  select(kpath3, kpath4, kpath5) %>%
  data.frame

path_density<-array(NA,dim=c(length(vill_w3_vec),4))
colnames(path_density)<-c("kpath3","kpath4","kpath5")

for(j in 1:length(vill_w3_vec)){
  tempc<-mean(vil_all_bl_diff_hh_kp2[j,])
  tp<-ifelse(tempc>1,vil_all[j,1]/tempc,vil_all[j,1])
  path_density[j,1]<-ifelse(tp>0,log(tp),tp)
  
  tempc<-mean(vil_all_bl_diff_hh_kp3[j,])
  tp<-ifelse(tempc>1,vil_all[j,2]/tempc,vil_all[j,2])
  path_density[j,2]<-ifelse(tp>0,log(tp),tp)
  
  tempc<-mean(vil_all_bl_diff_hh_kp4[j,])
  tp<-ifelse(tempc>1,vil_all[j,3]/tempc,vil_all[j,3])
  path_density[j,3]<-ifelse(tp>0,log(tp),tp)
  
  tempc<-mean(vil_all_bl_diff_hh_kp5[j,])
  tp<-ifelse(tempc>1,vil_all[j,4]/tempc,vil_all[j,4])
  path_density[j,4]<-ifelse(tp>0,log(tp),tp)
}

borrow_lend_diff_hh_undir_pathdens_w3_sds <- as.data.frame(path_density)

#------------------------ saving data ---------------------------
# Saving cleaned data to data folder; saved files are ready for analysis.
save(borrow_lend_diff_hh_undir_pathdens_w3_sds, file = "borrow_lend_diff_hh_undir_pathdens_w3_*.rda")
