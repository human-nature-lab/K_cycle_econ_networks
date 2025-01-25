# Trying different methods for calculating the expected number of cycles (denominator in cyclic 
# density calculation)

library(tidyverse)
library(dplyr)
library(igraph)
library(tidygraph)
library(ggraph)
library(intergraph)

library(ggplot2)

#-------------------------- Load data --------------------------------
setwd("/home/stl22/palmer_scratch")

# Create a data path for saving cleaned / augmented data sets
data_path <- "/home/stl22/palmer_scratch"

# Load node / edge list for borrowing ties
load(paste(c(data_path, "/borrow_lend_diff_hh_undir_edge_node_list_w1w3_20240728.rda"), collapse=""))

# Load generated network attributes to access kcycle prevalence per village
load(paste(c(data_path, "/net_attr_borrow_lend_undir_ds2_w1w3_20240728.rda"), collapse=""))

#-------------- method 2: sample_degseq() -----------------
# Example code from Yanick:
# (where "com" is the network)
# outd <- igraph::degree(com, mode = "out")
# ind <- igraph::degree(com, mode = "in")
#
# sim <- sample_degseq(outd, ind, method = "simple.no.multiple")
# sim <- asNetwork(sim)

# Create 100 random graphs that are generated using the same degree distribution
# each village network observed. Then calculate the observed number of cycles in
# each random graph generated and cyclic density for observed villages.

#------------ indiv borrow lend network
##Do 100 random villages per village
vc <- c(1:176)

vil_indiv_w1_bl_diff_hh_kcy2<-array(NA,dim=c(length(vc),100))
vil_indiv_w1_bl_diff_hh_kcy3<-array(NA,dim=c(length(vc),100))
vil_indiv_w1_bl_diff_hh_kcy4<-array(NA,dim=c(length(vc),100))
vil_indiv_w1_bl_diff_hh_kcy5<-array(NA,dim=c(length(vc),100))

library(sna, quietly = TRUE)

for(j in 1:176){
  edgelist_subset <- subset(borrow_lend_diff_hh_undir_edgelist_w1, village_code_w1==j)
  g <- graph_from_data_frame(edgelist_subset, directed=FALSE)
  g <- simplify(g)
  outd <- igraph::degree(g, mode = "out")
  ind <- igraph::degree(g, mode = "in")
  
  for(kk in 1:100) {
    sim <- sample_degseq(outd, ind, method = "simple.no.multiple")
    sim <- asNetwork(sim)
    # sim_adj = get.adjacency(sim,sparse=FALSE)
    
    
    apot.g=kcycle.census(sim, maxlen=5, mode = "graph", tabulate.by.vertex =TRUE,
                         cycle.comembership = c("none"))
    # detach("package:sna", unload = TRUE)
    
    vil_indiv_w1_bl_diff_hh_kcy2[j,kk] <- apot.g$cycle.count[1,1]
    vil_indiv_w1_bl_diff_hh_kcy3[j,kk] <- apot.g$cycle.count[2,1]
    vil_indiv_w1_bl_diff_hh_kcy4[j,kk] <- apot.g$cycle.count[3,1]
    vil_indiv_w1_bl_diff_hh_kcy5[j,kk] <- apot.g$cycle.count[4,1]
  }
  print(j)
}

### Cyclic density
# Data frame of the number of cycles per village
vil_all <- kstats_borrow_lend_diff_hh_undir_w1_df %>%
  select(kcycle2, kcycle3, kcycle4, kcycle5) #%>%
# data.frame()

cyc_density<-array(NA,dim=c(176,4))
colnames(cyc_density)<-c("kcycle2","kcycle3","kcycle4","kcycle5")

for(j in 1:176){
  tempc<-mean(vil_indiv_w1_bl_diff_hh_kcy2[j,])
  tp<-ifelse(tempc>1,vil_all[j,1]/tempc,vil_all[j,1])
  cyc_density[j,1]<-ifelse(tp>0,log(tp),tp)
  
  tempc<-mean(vil_indiv_w1_bl_diff_hh_kcy3[j,])
  tp<-ifelse(tempc>1,vil_all[j,2]/tempc,vil_all[j,2])
  cyc_density[j,2]<-ifelse(tp>0,log(tp),tp)
  
  tempc<-mean(vil_indiv_w1_bl_diff_hh_kcy4[j,])
  tp<-ifelse(tempc>1,vil_all[j,3]/tempc,vil_all[j,3])
  cyc_density[j,3]<-ifelse(tp>0,log(tp),tp)
  
  tempc<-mean(vil_indiv_w1_bl_diff_hh_kcy5[j,])
  tp<-ifelse(tempc>1,vil_all[j,4]/tempc,vil_all[j,4])
  cyc_density[j,4]<-ifelse(tp>0,log(tp),tp)
}

borrow_lend_diff_hh_undir_cycdens_w1_sds <- as.data.frame(cyc_density)

# #------------------ hh borrow lend network
# ##Do 100 random villages per village
# vc <- c(1:176)
# 
# vil_hh_w1_sds_kcy2<-array(NA,dim=c(length(vc),100))
# vil_hh_w1_sds_kcy3<-array(NA,dim=c(length(vc),100))
# vil_hh_w1_sds_kcy4<-array(NA,dim=c(length(vc),100))
# vil_hh_w1_sds_kcy5<-array(NA,dim=c(length(vc),100))
# 
# library(sna, quietly = TRUE)
# 
# for(j in 1:176){
#   edgelist_subset <- subset(hh_borrow_lend_edgelist_w1, village_code_w1==j)
#   g <- graph_from_data_frame(edgelist_subset, directed=TRUE)
#   g <- simplify(g)
#   outd <- igraph::degree(g, mode = "out")
#   ind <- igraph::degree(g, mode = "in")
#   
#   for(kk in 1:100) {
#     sim <- sample_degseq(outd, ind, method = "simple.no.multiple")
#     sim <- asNetwork(sim)
#     # sim_adj = get.adjacency(sim,sparse=FALSE)
#     
#     
#     apot.g=kcycle.census(sim, maxlen=5, mode = "digraph", tabulate.by.vertex =TRUE,
#                          cycle.comembership = c("none"))
#     # detach("package:sna", unload = TRUE)
#     
#     vil_hh_w1_sds_kcy2[j,kk] <- apot.g$cycle.count[1,1]
#     vil_hh_w1_sds_kcy3[j,kk] <- apot.g$cycle.count[2,1]
#     vil_hh_w1_sds_kcy4[j,kk] <- apot.g$cycle.count[3,1]
#     vil_hh_w1_sds_kcy5[j,kk] <- apot.g$cycle.count[4,1]
#   }
# }


#------------------------ saving data ---------------------------
# Saving cleaned data to data folder; saved files are ready for analysis.
save(vil_indiv_w1_bl_diff_hh_kcy2, vil_indiv_w1_bl_diff_hh_kcy3, vil_indiv_w1_bl_diff_hh_kcy4, vil_indiv_w1_bl_diff_hh_kcy5,
     file = "exp_cyc_bl_diff_hh_indiv_w1_20240728.rda")

save(borrow_lend_diff_hh_undir_cycdens_w1_sds,
     file = file.path(data_path, "cycdens_borrow_lend_diff_hh_undir_w1_20240728.rda"))

# save(vil_hh_w1_sds_kcy2, vil_hh_w1_sds_kcy3, vil_hh_w1_sds_kcy4, vil_hh_w1_sds_kcy5,
#      file = file.path(data_path, "exp_cycles_SDS_hh_w1_20240306.rda"))
