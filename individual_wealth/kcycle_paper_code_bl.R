
#Code written By Shivkumar Vishnempet Shridhar (contact @ shivkumar.vishnempet@yale.edu)
#Paper with Selena Lee, Georgios Iosifidis, Yanick Charette, & Nicholas A. Christakis (23 Jan 2025)
#This script is for 1)generating cycles and cyclic participation of every individual 2) Creating cycle composition scores for every individual
#3) Investigating relaitonship between cycle particpation & cycle composition with individual wealth/gain in wealth.
#This script is only for Borrow-lend networks

library(tidyverse)
library(DescTools)
library(skimr)
library(lme4)
library(jtools)
library(labelled)
library(ggcorrplot)
library(dplyr)
library(sna)
library(igraph)
require(foreach)
require(doParallel)
library(pheatmap)
library(grid)

setwd("/K-cycle_MS/")

# Create a data path for saving cleaned / augmented data sets
data_path <- "/data"
figure_path <- "/figures"

# load data
load("/vill_df_*.rda")
load("/hh_df_*.rda")

# Load raw MCA coordinates
mca_w1_df <- read_csv(paste(data_path, "/hh_mca_w1.csv", sep=""))#MCA coordinates Wave 1
mca_w2_df <- read_csv(paste(data_path, "/hh_mca_w2.csv", sep=""))#MCA coordinates Wave 2

# Respondent-level data
resp_w1_original <- read_csv("hw1_resp_v8_*.csv")#Wave 1 respondent data
resp_w2_original <- read_csv("hw3_resp_v3_*.csv")#Wave 2 respondent data

#Filtering and cleaning dataset


#### Optional block
#Renaming necessary columns
resp_w2_original<-resp_w2_original%>% rename(village_code_w2 = village_code_W2)
conn_w2_original<-conn_w2_original%>% rename(village_code_w2 = village_code_W2)
####


# Subsetting network data for all villages across both waves data
conn_w1 <- subset(conn_w1_original, alter_source==1 & village_code_w1!=0 & same_village==1 & same_building==0)
conn_w2 <- subset(conn_w2_original, alter_source==1 & village_code_w2!=0 & same_village==1 & same_building==0)

#--------------- 1: Borrow-lend individual networks ----------------
#Wave 1
nodelist_w1 <- conn_w1 %>%
  subset(select=c(ego, village_code_w1)) %>%
  dplyr::rename(node = ego) %>%
  bind_rows(conn_w1 %>%
              subset(select=c(alter, village_code_w1)) %>%
              dplyr::rename(node = alter)) %>%
  distinct(node, .keep_all=TRUE) %>%
  rowid_to_column("id")

nodelist_w1 <- nodelist_w1 %>%
  mutate(nodelist_w1%>%inner_join(resp_w1_original,by=join_by(node==respondent_master_id)) %>%
  subset(select=c(building_id_w1,gender))%>%
  rename(building_id=building_id_w1))

nodelist_w1 <- nodelist_w1 %>%
  mutate(nodelist_w1%>%left_join(hh_w1,by=join_by(building_id==building_id)) %>%
           subset(select=household_wealth_index_w1)%>%
           rename(HHW=household_wealth_index_w1))


borrow_lend_ind_edgelist_w1 <- conn_w1 %>%
  subset(same_village==1 & relationship %in% c("trust_borrow_money"),
         select = c(ego, alter, relationship, village_code_w1)) %>%
  bind_rows(conn_w1 %>%
              subset(same_village==1 & relationship %in% c("trust_lend_money"),
                     select=c(alter, ego, relationship, village_code_w1)) %>%
              rename(ego = alter,
                     alter = ego)) %>%
  left_join(nodelist_w1 %>%
              select(id, node), by = c("ego" = "node")) %>%
  dplyr::rename(from = id) %>%
  inner_join(nodelist_w1 %>%
               select(id, node), by = c("alter" = "node")) %>%
  dplyr::rename(to = id) %>%
  select(from, to, village_code_w1) %>%
  distinct(from, to, .keep_all = TRUE) %>%
  subset(from != to) %>%
  rename(village_code = village_code_w1)

#Wave 2
nodelist_w2 <- conn_w2 %>%
  subset(select=c(ego, village_code_w2)) %>%
  dplyr::rename(node = ego) %>%
  bind_rows(conn_w2 %>%
              subset(select=c(alter, village_code_w2)) %>%
              dplyr::rename(node = alter)) %>%
  distinct(node, .keep_all=TRUE) %>%
  rowid_to_column("id")

nodelist_w2 <- nodelist_w2 %>%
  mutate(nodelist_w2%>%inner_join(resp_w2_original,by=join_by(node==respondent_master_id)) %>%
           subset(select=c(building_id_w2,gender))%>%
           rename(building_id=building_id_w2))

nodelist_w2 <- nodelist_w2 %>%
  mutate(nodelist_w2%>%left_join(hh_w1,by=join_by(building_id==building_id)) %>%
           subset(select=household_wealth_index_w2)%>%
           rename(HHW=household_wealth_index_w2))

borrow_lend_ind_edgelist_w2 <- conn_w2 %>%
  subset(same_village==1 & relationship %in% c("trust_borrow_money"),
         select = c(ego, alter, relationship, village_code_w2)) %>%
  bind_rows(conn_w2 %>%
              subset(same_village==1 & relationship %in% c("trust_lend_money"),
                     select=c(alter, ego, relationship, village_code_w2)) %>%
              rename(ego = alter,
                     alter = ego)) %>%
  left_join(nodelist_w2 %>%
              select(id, node), by = c("ego" = "node")) %>%
  dplyr::rename(from = id) %>%
  inner_join(nodelist_w2 %>%
               select(id, node), by = c("alter" = "node")) %>%
  dplyr::rename(to = id) %>%
  select(from, to, village_code_w2) %>%
  distinct(from, to, .keep_all = TRUE) %>%
  subset(from != to) %>%
  rename(village_code = village_code_w2)



#Function for computing cycles and paths -- works

compute_network_metrics_w1<-function(bl_node,bl_edge,i){
  #Subsetting edge list by village
  bl_edge_graph <-subset(bl_edge, village_code==i)
  #bl_node_graph<-bl_node
  #bl_node_graph <-subset(bl_node, village_code_w1==i)#Nodes need to be matched!
  #Preparing graph 
  n <- graph_from_data_frame(bl_edge_graph, directed=F)
  n <- simplify(n)
  num_nodes <- length(V(n))
  #Preparing table for network metrics including cycle quantity & quality
  network_factor_names<-c("Node_id","building_id","Household_wealth","Degree","Transitivity","Eigen","Betweeness",
                                                 "kcycle3","kcycle4","kcycle5",
                                                 "kcy3_wealth_avg","kcy4_wealth_avg","kcy5_wealth_avg",
                                                 "path3","path4","path5")
  network_factors<-as.data.frame(array(NA,dim=c(num_nodes,length(network_factor_names))))
  colnames(network_factors)<-network_factor_names
  network_factors$Degree<-degree(n)/num_nodes
  network_factors$Transitivity<-transitivity(n,type="local")
  eigen_n<-eigen_centrality(n,directed = F)
  network_factors$Eigen<-eigen_n$vector
  network_factors$Betweeness<-betweenness(n,directed = F)/num_nodes
  network_factors$Household_wealth<-bl_node$HHW[as.numeric(V(n)$name)]
  network_factors$Node_id<-bl_node$node[as.numeric(V(n)$name)]
  network_factors$building_id<-bl_node$building_id[as.numeric(V(n)$name)]
  
  cycle_n=as_adjacency_matrix(n,sparse=FALSE)
  
  cycle_n_tab=kcycle.census(cycle_n, maxlen=5, mode = "graph", tabulate.by.vertex =TRUE, cycle.comembership = c("bylength"))
  path_n_tab=kpath.census(cycle_n, maxlen=5, mode = "graph", tabulate.by.vertex =TRUE)
  #cycle_n_member<-kcycle.census(cycle_n, maxlen=5, mode = "graph", tabulate.by.vertex =F, cycle.comembership = c("none"))
  
  #No same building ties present
  #Computing network metrics for each node in every village
  
    
  for(j in 1:num_nodes){
    
      #Number of cycles -- since first column is aggregate
      network_factors$kcycle3[j]<-cycle_n_tab$cycle.count[2,j+1]#sum(cycle_n_tab$cycle.comemb[2,j,])#K-cycle 3
      network_factors$kcycle4[j]<-cycle_n_tab$cycle.count[3,j+1]#sum(cycle_n_tab$cycle.comemb[3,j,])#K-cycle 4
      network_factors$kcycle5[j]<-cycle_n_tab$cycle.count[4,j+1]#sum(cycle_n_tab$cycle.comemb[4,j,])#K-cycle 5

      #Co-members in every cycle w.r.t to node
      alter_id_node<-match(bl_edge$to[which(bl_edge$from==j)],colnames(cycle_n))

      #Node indices of cycle-co-members
      cycle_mem3<-as.numeric(which(cycle_n_tab$cycle.comemb[2,j,]>0))
      cycle_mem4<-as.numeric(which(cycle_n_tab$cycle.comemb[3,j,]>0))
      cycle_mem5<-as.numeric(which(cycle_n_tab$cycle.comemb[4,j,]>0))
      
      #Caclulating wealth
      wealth_3<-network_factors$Household_wealth[cycle_mem3]
      network_factors$kcy3_wealth_avg[j]<-ifelse(length(wealth_3)==0,0,mean(wealth_3[complete.cases(wealth_3)]))
      wealth_4<-network_factors$Household_wealth[cycle_mem4]
      network_factors$kcy4_wealth_avg[j]<-ifelse(length(wealth_4)==0,0,mean(wealth_4[complete.cases(wealth_4)]))
      wealth_5<-network_factors$Household_wealth[cycle_mem5]
      network_factors$kcy5_wealth_avg[j]<-ifelse(length(wealth_5)==0,0,mean(wealth_5[complete.cases(wealth_5)]))
      
      #Cacluclating number of paths -- since first column is aggregate
      network_factors$path3[j]<-path_n_tab$path.count[3,j+1]#sum(cycle_n_tab$cycle.comemb[2,j,])#K-cycle 3
      network_factors$path4[j]<-path_n_tab$path.count[4,j+1]#sum(cycle_n_tab$cycle.comemb[3,j,])#K-cycle 4
      network_factors$path5[j]<-path_n_tab$path.count[5,j+1]#sum(cycle_n_tab$cycle.comemb[4,j,])#K-cycle 5
      
  }
  
  
  return(network_factors)
}

stopCluster(cl)

#Compute cycles

#Number of cluster nodes for parallel computing
cl <- makeCluster(20)
registerDoParallel(cl)

village_range<-max(nodelist_w1$village_code_w1,nodelist_w2$village_code_w2)

network_factors_vil_w1 <- foreach(sampleNr=1:village_range, .combine=rbind) %dopar% {
  library(sna)
  library(igraph)
  if(sampleNr%in%nodelist_w1$village_code_w1){
    print(compute_network_metrics_w1(nodelist_w1,borrow_lend_ind_edgelist_w1,sampleNr)) 
  }
  
}

network_factors_vil_w2 <- foreach(sampleNr=1:village_range, .combine=rbind) %dopar% {
  library(sna)
  library(igraph)
  if(sampleNr%in%nodelist_w2$village_code_w2){
    print(compute_network_metrics_w1(nodelist_w2,borrow_lend_ind_edgelist_w2,sampleNr)) 
  }
  #print(compute_network_metrics_w1(nodelist_w2,borrow_lend_ind_edgelist_w2,sampleNr))
  
}



#Wealth measures -- Adding in all wealth measures
#zm4_w1$resp_id<-rownames(zm4_w1)
#Wave 1
network_factors_vil_w1<-network_factors_vil_w1%>%
                        mutate(network_factors_vil_w1%>%left_join(zm4_w1,by=join_by(Node_id==resp_id)) %>%
                                 subset(select=c("mca_coord","l0700","l0800","l0900a","l0900b",
                                          "l0900c","l0900d","l0900e","l0900g","l1000",
                                          "l1100","l1200","l1300","l1400","l1500","l1600",
                                          "l1700","food_insecurity","Flush (toilet type)",
                                          "Furnace/firebox with a chimney (cooking stove)",
                                          "Ceramic (floor)","Cement blocks (roof)",
                                          "No facility (outdoors; toilet type)",
                                          "Furnace/firebox without a chimney (cooking stove)",
                                          "Earth/sand (floor)","There aren't windows")))


network_factors_vil_w1$Transitivity<-ifelse(is.nan(network_factors_vil_w1$Transitivity)==T,0,network_factors_vil_w1$Transitivity)
network_factors_vil_w1<-network_factors_vil_w1[which(network_factors_vil_w1==T),]

#Wave 2
network_factors_vil_w2<-network_factors_vil_w2%>%
  mutate(network_factors_vil_w2%>%left_join(zm4_w2,by=join_by(Node_id==resp_id)) %>%
           subset(select=c("mca_coord","l0700","l0800","l0900a","l0900b",
                           "l0900c","l0900d","l0900e","l0900g","l1000",
                           "l1100","l1200","l1300","l1400","l1500","l1600",
                           "l1700","food_insecurity","Flush (toilet type)",
                           "Furnace/firebox with a chimney (cooking stove)",
                           "Ceramic (floor)","Cement blocks (roof)",
                           "No facility (outdoors; toilet type)",
                           "Furnace/firebox without a chimney (cooking stove)",
                           "Earth/sand (floor)","There aren't windows")))


network_factors_vil_w2$Transitivity<-ifelse(is.nan(network_factors_vil_w2$Transitivity)==T,0,network_factors_vil_w2$Transitivity)
network_factors_vil_w2<-network_factors_vil_w2[which(network_factors_vil_w2==T),]


#Regression

com_ind<-intersect(rownames(nodelist_w1),rownames(nodelist_w2))
nodelist_w13c_use<-nodelist_w1[match(com_ind,rownames(nodelist_w1)),match(norm_list[1:(length(norm_list)-2)],colnames(nodelist_w1))]#
colnames(nodelist_w13c_use)[3:length(colnames(nodelist_w13c_use))]<-paste0(colnames(nodelist_w13c_use)[3:length(colnames(nodelist_w13c_use))],"_W1")
nodelist_w13c_use<-cbind(nodelist_w13c_use,nodelist_w2[match(com_ind,rownames(nodelist_w2)),match(norm_list[3:(length(norm_list)-2)],colnames(nodelist_w2))])#
colnames(nodelist_w13c_use)[11:length(colnames(nodelist_w13c_use))]<-paste0(colnames(nodelist_w13c_use)[11:length(colnames(nodelist_w13c_use))],"_W2")

ns3_items_ez_ind_rev<-array(NA,dim=c(4,20))
colnames(ns3_items_ez_ind_rev)<-c("HHW","MCA coordinate","Food insecurity","Shared toilet","Electricity","Radio","TV","Cell/mobile phone","Non mobile phone","No electronics","Separate rooms","Sleeping rooms","Flush (toilet type)","Furnace/firebox with a chimney (cooking stove)","Ceramic (floor)","Cement blocks (roof)","No facility (outdoors; toilet type)","Furnace/firebox without a chimney (cooking stove)","Earth/sand (floor)","There aren't windows")
rownames(ns3_items_ez_ind_rev)<-c("Wave 1","Wave 3","Wave 3--1","Wave 1+3")
ns4_items_ez_ind_rev<-array(NA,dim=c(4,20))
colnames(ns4_items_ez_ind_rev)<-c("HHW","MCA coordinate","Food insecurity","Shared toilet","Electricity","Radio","TV","Cell/mobile phone","Non mobile phone","No electronics","Separate rooms","Sleeping rooms","Flush (toilet type)","Furnace/firebox with a chimney (cooking stove)","Ceramic (floor)","Cement blocks (roof)","No facility (outdoors; toilet type)","Furnace/firebox without a chimney (cooking stove)","Earth/sand (floor)","There aren't windows")
rownames(ns4_items_ez_ind_rev)<-c("Wave 1","Wave 3","Wave 3--1","Wave 1+3")
ns5_items_ez_ind_rev<-array(NA,dim=c(4,20))
colnames(ns5_items_ez_ind_rev)<-c("HHW","MCA coordinate","Food insecurity","Shared toilet","Electricity","Radio","TV","Cell/mobile phone","Non mobile phone","No electronics","Separate rooms","Sleeping rooms","Flush (toilet type)","Furnace/firebox with a chimney (cooking stove)","Ceramic (floor)","Cement blocks (roof)","No facility (outdoors; toilet type)","Furnace/firebox without a chimney (cooking stove)","Earth/sand (floor)","There aren't windows")
rownames(ns5_items_ez_ind_rev)<-c("Wave 1","Wave 3","Wave 3--1","Wave 1+3")

ns3_items_pv_ind_rev<-array(NA,dim=c(4,20))
colnames(ns3_items_pv_ind_rev)<-c("HHW","MCA coordinate","Food insecurity","Shared toilet","Electricity","Radio","TV","Cell/mobile phone","Non mobile phone","No electronics","Separate rooms","Sleeping rooms","Flush (toilet type)","Furnace/firebox with a chimney (cooking stove)","Ceramic (floor)","Cement blocks (roof)","No facility (outdoors; toilet type)","Furnace/firebox without a chimney (cooking stove)","Earth/sand (floor)","There aren't windows")
rownames(ns3_items_pv_ind_rev)<-c("Wave 1","Wave 3","Wave 3--1","Wave 1+3")
ns4_items_pv_ind_rev<-array(NA,dim=c(4,20))
colnames(ns4_items_pv_ind_rev)<-c("HHW","MCA coordinate","Food insecurity","Shared toilet","Electricity","Radio","TV","Cell/mobile phone","Non mobile phone","No electronics","Separate rooms","Sleeping rooms","Flush (toilet type)","Furnace/firebox with a chimney (cooking stove)","Ceramic (floor)","Cement blocks (roof)","No facility (outdoors; toilet type)","Furnace/firebox without a chimney (cooking stove)","Earth/sand (floor)","There aren't windows")
rownames(ns4_items_pv_ind_rev)<-c("Wave 1","Wave 3","Wave 3--1","Wave 1+3")
ns5_items_pv_ind_rev<-array(NA,dim=c(4,20))
colnames(ns5_items_pv_ind_rev)<-c("HHW","MCA coordinate","Food insecurity","Shared toilet","Electricity","Radio","TV","Cell/mobile phone","Non mobile phone","No electronics","Separate rooms","Sleeping rooms","Flush (toilet type)","Furnace/firebox with a chimney (cooking stove)","Ceramic (floor)","Cement blocks (roof)","No facility (outdoors; toilet type)","Furnace/firebox without a chimney (cooking stove)","Earth/sand (floor)","There aren't windows")
rownames(ns5_items_pv_ind_rev)<-c("Wave 1","Wave 3","Wave 3--1","Wave 1+3")

norm_list<-c("Village","kcycle3","kcycle4","kcycle5","path3","path4","path5","wave","size")
d_list3<-c("HHW","mca_coord","l0800","l0900a","l0900b","l0900c",
           "l0900d","l0900e","l0900g","l1200","l1700","Flush (toilet type)",
           "Furnace/firebox with a chimney (cooking stove)","Ceramic (floor)",
           "Cement blocks (roof)","No facility (outdoors; toilet type)",
           "Furnace/firebox without a chimney (cooking stove)","Earth/sand (floor)","There aren't windows")

i<-0

for(i in 1:length(d_list3)){
  #i<-i+1
  #Wave 1
  tempp1<-nodelist_w1[match(com_ind,rownames(nodelist_w1)),match(c(norm_list2,d_list3[i]),colnames(nodelist_w1))]#
  tp<-complete.cases(tempp1)
  nodelist_w1_use<-tempp1[tp,]
  colnames(nodelist_w1_use)[dim(nodelist_w1_use)[2]]<-"item"
  temp<-cbind(nodelist_w1[match(com_ind,rownames(nodelist_w1))[tp],match(c("HHW","Village"),colnames(nodelist_w1))])#Here
  vil_w<-array(NA,dim=dim(temp)[1])
  for(ij in 1:dim(temp)[1]){
    tp<-temp[which(temp[,2]==temp[ij,2]),1]
    vil_w[ij]<-mean(as.numeric(as.character(tp[which(is.na(tp)==F)])))
  }
  nodelist_w1_use$vil_w<-vil_w
  
  #K-Path removed because it highlt correlates with kcycle , so it brings in multi-collinearity -- instead computing separate model for path
  ns3_items_ez_ind_rev[1,i]<-coef(summary(lmer(as.numeric(as.character(kcycle3))~as.numeric(as.character(item))+as.numeric(as.character(path3))+as.numeric(as.character(Degree))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w1_use)))[2,1]
  ns4_items_ez_ind_rev[1,i]<-coef(summary(lmer(as.numeric(as.character(kcycle4))~as.numeric(as.character(item))+as.numeric(as.character(path4))+as.numeric(as.character(Degree))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w1_use)))[2,1]
  ns5_items_ez_ind_rev[1,i]<-coef(summary(lmer(as.numeric(as.character(kcycle5))~as.numeric(as.character(item))+as.numeric(as.character(path5))+as.numeric(as.character(Degree))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w1_use)))[2,1]
  
  ns3_items_pv_ind_rev[1,i]<-coef(summary(lmer(as.numeric(as.character(kcycle3))~as.numeric(as.character(item))+as.numeric(as.character(path3))+as.numeric(as.character(Degree))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w1_use)))[2,5]
  ns4_items_pv_ind_rev[1,i]<-coef(summary(lmer(as.numeric(as.character(kcycle4))~as.numeric(as.character(item))+as.numeric(as.character(path4))+as.numeric(as.character(Degree))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w1_use)))[2,5]
  ns5_items_pv_ind_rev[1,i]<-coef(summary(lmer(as.numeric(as.character(kcycle5))~as.numeric(as.character(item))+as.numeric(as.character(path5))+as.numeric(as.character(Degree))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w1_use)))[2,5]
  
  #Wave 2
  tempp1<-nodelist_w2[match(com_ind,rownames(nodelist_w2)),match(c(norm_list2,d_list3[i]),colnames(nodelist_w1))]#
  tp<-complete.cases(tempp1)
  nodelist_w2_use<-tempp1[tp,]
  colnames(nodelist_w2_use)[dim(nodelist_w2_use)[2]]<-"item"
  temp<-cbind(nodelist_w1[match(com_ind,rownames(nodelist_w1))[tp],match(c("HHW","Village"),colnames(nodelist_w2))])#Here
  vil_w<-array(NA,dim=dim(temp)[1])
  for(ij in 1:dim(temp)[1]){
    tp<-temp[which(temp[,2]==temp[ij,2]),1]
    vil_w[ij]<-mean(as.numeric(as.character(tp[which(is.na(tp)==F)])))
  }
  nodelist_w2_use$vil_w<-vil_w
  
  ns3_items_ez_ind_rev[2,i]<-coef(summary(lmer(as.numeric(as.character(kcycle3))~as.numeric(as.character(item))+as.numeric(as.character(path3))+as.numeric(as.character(Degree))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w2_use)))[2,1]
  ns4_items_ez_ind_rev[2,i]<-coef(summary(lmer(as.numeric(as.character(kcycle4))~as.numeric(as.character(item))+as.numeric(as.character(path4))+as.numeric(as.character(Degree))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w2_use)))[2,1]
  ns5_items_ez_ind_rev[2,i]<-coef(summary(lmer(as.numeric(as.character(kcycle5))~as.numeric(as.character(item))+as.numeric(as.character(path5))+as.numeric(as.character(Degree))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w2_use)))[2,1]
  
  ns3_items_pv_ind_rev[2,i]<-coef(summary(lmer(as.numeric(as.character(kcycle3))~as.numeric(as.character(item))+as.numeric(as.character(path3))+as.numeric(as.character(Degree))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w2_use)))[2,5]
  ns4_items_pv_ind_rev[2,i]<-coef(summary(lmer(as.numeric(as.character(kcycle4))~as.numeric(as.character(item))+as.numeric(as.character(path4))+as.numeric(as.character(Degree))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w2_use)))[2,5]
  ns5_items_pv_ind_rev[2,i]<-coef(summary(lmer(as.numeric(as.character(kcycle5))~as.numeric(as.character(item))+as.numeric(as.character(path5))+as.numeric(as.character(Degree))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w2_use)))[2,5]
  
  #Wave 3--1
  tempp1<-cbind(nodelist_w13c_use,nodelist_w1[match(com_ind,rownames(nodelist_w1)),match(d_list3[i],colnames(nodelist_w1))])
  colnames(tempp1)[dim(tempp1)[2]]<-"item"
  temp<-cbind(nodelist_w13c_use[,match(c("HHW_W1","Village"),colnames(nodelist_w13s_use))])#Here
  vil_w<-array(NA,dim=dim(temp)[1])
  for(ij in 1:dim(temp)[1]){
    tp<-temp[which(temp[,2]==temp[ij,2]),1]
    vil_w[ij]<-mean(as.numeric(as.character(tp[which(is.na(tp)==F)])))
  }
  tempp1$vil_w<-vil_w
  nodelist_w13c_use2<-tempp1[complete.cases(tempp1),]
  
  ns3_items_ez_ind_rev[3,i]<-coef(summary(lmer(as.numeric(as.character(kcycle3_W2))~as.numeric(as.character(item))+as.numeric(as.character(path3_W1))+as.numeric(as.character(Degree))+as.numeric(as.character(kcycle3_W1))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w13c_use2)))[2,1]
  ns4_items_ez_ind_rev[3,i]<-coef(summary(lmer(as.numeric(as.character(kcycle4_W2))~as.numeric(as.character(item))+as.numeric(as.character(path4_W1))+as.numeric(as.character(Degree))+as.numeric(as.character(kcycle4_W1))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w13c_use2)))[2,1]
  ns5_items_ez_ind_rev[3,i]<-coef(summary(lmer(as.numeric(as.character(kcycle5_W2))~as.numeric(as.character(item))+as.numeric(as.character(path5_W1))+as.numeric(as.character(Degree))+as.numeric(as.character(kcycle5_W1))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w13c_use2)))[2,1]
  
  ns3_items_pv_ind_rev[3,i]<-coef(summary(lmer(as.numeric(as.character(kcycle3_W2))~as.numeric(as.character(item))+as.numeric(as.character(path3_W1))+as.numeric(as.character(Degree))+as.numeric(as.character(kcycle3_W1))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w13c_use2)))[2,5]
  ns4_items_pv_ind_rev[3,i]<-coef(summary(lmer(as.numeric(as.character(kcycle4_W2))~as.numeric(as.character(item))+as.numeric(as.character(path4_W1))+as.numeric(as.character(Degree))+as.numeric(as.character(kcycle4_W1))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w13c_use2)))[2,5]
  ns5_items_pv_ind_rev[3,i]<-coef(summary(lmer(as.numeric(as.character(kcycle5_W2))~as.numeric(as.character(item))+as.numeric(as.character(path5_W1))+as.numeric(as.character(Degree))+as.numeric(as.character(kcycle5_W1))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w13c_use2)))[2,5]
  
  
  
  # #Wave 1+3
  #tempp1<-nodelist_w13ss[,c(3,7:14,18,45,d_list3[i])]
  tempp1<-nodelist_w13ss[,match(c(norm_list2,d_list3[i]),colnames(nodelist_w13ss))]#
  nodelist_w13ss2<-tempp1[complete.cases(tempp1),]
  colnames(nodelist_w13ss2)[dim(nodelist_w13ss2)[2]]<-"item"
  temp<-cbind(nodelist_w13ss[complete.cases(tempp1),match(c("HHW","Village"),colnames(nodelist_w13ss))])#Here
  vil_w<-array(NA,dim=dim(temp)[1])
  for(ij in 1:dim(temp)[1]){
    tp<-temp[which(temp[,2]==temp[ij,2]),1]
    vil_w[ij]<-mean(as.numeric(as.character(tp[which(is.na(tp)==F)])))
  }
  nodelist_w13ss2$vil_w<-vil_w
  
  ns3_items_ez_ind_rev[4,i]<-coef(summary(lmer(as.numeric(as.character(kcycle3))~as.numeric(as.character(item))+as.numeric(as.character(path3))+as.numeric(as.character(Degree))+as.numeric(as.character(wave))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w13ss2)))[2,1]
  ns4_items_ez_ind_rev[4,i]<-coef(summary(lmer(as.numeric(as.character(kcycle4))~as.numeric(as.character(item))+as.numeric(as.character(path4))+as.numeric(as.character(Degree))+as.numeric(as.character(wave))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w13ss2)))[2,1]
  ns5_items_ez_ind_rev[4,i]<-coef(summary(lmer(as.numeric(as.character(kcycle5))~as.numeric(as.character(item))+as.numeric(as.character(path5))+as.numeric(as.character(Degree))+as.numeric(as.character(wave))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w13ss2)))[2,1]
  
  ns3_items_pv_ind_rev[4,i]<-coef(summary(lmer(as.numeric(as.character(kcycle3))~as.numeric(as.character(item))+as.numeric(as.character(path3))+as.numeric(as.character(Degree))+as.numeric(as.character(wave))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w13ss2)))[2,5]
  ns4_items_pv_ind_rev[4,i]<-coef(summary(lmer(as.numeric(as.character(kcycle4))~as.numeric(as.character(item))+as.numeric(as.character(path4))+as.numeric(as.character(Degree))+as.numeric(as.character(wave))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w13ss2)))[2,5]
  ns5_items_pv_ind_rev[4,i]<-coef(summary(lmer(as.numeric(as.character(kcycle5))~as.numeric(as.character(item))+as.numeric(as.character(path5))+as.numeric(as.character(Degree))+as.numeric(as.character(wave))+as.numeric(as.character(vil_w))+(1|Village),
                                               data = nodelist_w13ss2)))[2,5]
  
  
  #}
}







#Plots


#Regression results

# Individual 

items_ind_ez<-rbind(ns3_items_ez_ind,rbind(ns4_items_ez_ind,ns5_items_ez_ind))
items_ind_pv<-rbind(ns3_items_pv_ind,rbind(ns4_items_pv_ind,ns5_items_pv_ind))

items_ind_ez<-items_ind_ez[c(1:3,5,6:8,10,11:13,15),]
items_ind_pv<-items_ind_pv[c(1:3,5,6:8,10,11:13,15),]

imm<-unlist(as.numeric(items_ind_pv))
imm2<-p.adjust(imm,method='BH')

items_ind_fdr<-array(NA,dim=c(dim(items_ind_pv)[1],dim(items_ind_pv)[2]))

for(i in 1:dim(items_ind_fdr)[2]){
  items_ind_fdr[,i]<-imm2[(((i-1)*dim(items_ind_fdr)[1])+1):(i*dim(items_ind_fdr)[1])]
}

disp_fdr<-ifelse(items_ind_fdr<0.05,ifelse(items_ind_ez>0,"+","-"),"")


#Add rownames and column names
items_ind_ez2<-ifelse(disp_fdr!="",items_ind_ez,0)
effect_size_phen_sig_plot<-t(items_ind_ez2)
effect_size_phen_sig_plot<-as.data.frame(effect_size_phen_sig_plot)
colnames(effect_size_phen_sig_plot)[1:4]<-c('Variable (W1) ~ cycles (W1)','Variable (W3) ~ cycles (W3)','Variable (W3) ~ cycles (W1) + Wealth(W1)','Variable (W1+W3) ~ cycles (W1+W3)')
colnames(effect_size_phen_sig_plot)[5:8]<-c('Variable (W1) ~ cycles (W1)','Variable (W3) ~ cycles (W3)','Variable (W3) ~ cycles (W1) + Wealth(W1)','Variable (W1+W3) ~ cycles (W1+W3)')
colnames(effect_size_phen_sig_plot)[9:12]<-c('Variable (W1) ~ cycles (W1)','Variable (W3) ~ cycles (W3)','Variable (W3) ~ cycles (W1) + Wealth(W1)','Variable (W1+W3) ~ cycles (W1+W3)')
#colnames(effect_size_phen_sig_plot)[13:16]<-c('Variable (W1) ~ cycles (W1)','Variable (W3) ~ cycles (W3)','Variable (W3) ~ cycles (W1) + Wealth(W1)','Variable (W1+W3) ~ cycles (W1+W3)')

items_ind_ez2<-as.data.frame(items_ind_ez2)
rownames(effect_size_phen_sig_plot)<-colnames(items_ind_ez2)# was as.character(mb_samp_sp_name[which(pval_chk_all>5)])
rownames(effect_size_phen_sig_plot)[1:2]<-c("Household Wealth Index","MCA coordinate")
rownames(effect_size_phen_sig_plot)[3:dim(items_ind_ez)[2]]<-colnames(items_ind_ez)[3:dim(items_ind_ez)[2]]
rownames(effect_size_phen_sig_plot)[8]<-"No phone"

myb2<-c(0.0000000125,0.0012500000,0.0065105438,0.0112500000,0.0143750000,0.0179000000,0.0342018138,0.0720843541,0.15,1.5)
myBreaks<-c(myb2[c(length(myb2):-1:1)]*-1,myb2)
myColor<-colorRampPalette(c("#cb181d","white", "#33a02c"))(19)


pdf('regression_items_bl.pdf',width=9,height=7)
pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,
         labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,
         border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,
         color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_col = 0,
         treeheight_row = 0,cluster_rows = F,cluster_cols = F)
dev.off()


ind_phen<-c(4,8,12,3,7,11)
ind_reg<-c(1,5:7,13:16)#8:10,17:20,4#No MCA coordinate and FI -- doesnt add anything new 

rowr<-rownames(effect_size_phen_sig_plot)[ind_reg]
rowr[5:8]<-c("Toilet type","Cooking stove type","Floor type","Roof type")

pdf('figure_path/figure_1a.pdf',width=5,height=6)#was ind_phen[4:6]
pheatmap(effect_size_phen_sig_plot[ind_reg,ind_phen],annotation_row = NULL,annotation_names_row = T,
         labels_row = rowr,labels_col = colnames(effect_size_phen_sig_plot)[ind_phen],legend = T,fontsize_number = 15,
         border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,
         color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[ind_phen,ind_reg]),number_color = "black",treeheight_col = 0,
         treeheight_row = 0,cluster_rows = F,cluster_cols = F)
dev.off()




#Reverse regression results

# 

items_ind_ez_rev<-rbind(ns3_items_ez_ind_rev,rbind(ns4_items_ez_ind_rev,ns5_items_ez_ind_rev))
items_ind_pv_rev<-rbind(ns3_items_pv_ind_rev,rbind(ns4_items_pv_ind_rev,ns5_items_pv_ind_rev))

imm<-unlist(as.numeric(items_ind_pv_rev))
imm2<-p.adjust(imm,method='BH')

items_ind_fdr_rev<-array(NA,dim=c(dim(items_ind_pv_rev)[1],dim(items_ind_pv_rev)[2]))

for(i in 1:dim(items_ind_fdr_rev)[2]){
  items_ind_fdr_rev[,i]<-imm2[(((i-1)*dim(items_ind_fdr_rev)[1])+1):(i*dim(items_ind_fdr_rev)[1])]
}

disp_fdr_rev<-ifelse(items_ind_fdr_rev<0.05,ifelse(items_ind_ez_rev>0,"+","-"),"")


#Add rownames and column names
items_ind_ez2_rev<-ifelse(disp_fdr_rev!="",items_ind_ez_rev,0)
effect_size_phen_sig_plot<-t(items_ind_ez2_rev)
effect_size_phen_sig_plot<-as.data.frame(effect_size_phen_sig_plot)
colnames(effect_size_phen_sig_plot)[1:4]<-c('Cycles (W1) ~ Variable (W1)','Cycles (W3) ~ Variable (W3)','Cycles (W3) ~ Variable (W1) + Cycles(W1)','Cycles (W1+W3) ~ Variable (W1+W3)')
colnames(effect_size_phen_sig_plot)[5:8]<-c('Cycles (W1) ~ Variable (W1)','Cycles (W3) ~ Variable (W3)','Cycles (W3) ~ Variable (W1) + Cycles(W1)','Cycles (W1+W3) ~ Variable (W1+W3)')
colnames(effect_size_phen_sig_plot)[9:12]<-c('Cycles (W1) ~ Variable (W1)','Cycles (W3) ~ Variable (W3)','Cycles (W3) ~ Variable (W1) + Cycles(W1)','Cycles (W1+W3) ~ Variable (W1+W3)')
#colnames(effect_size_phen_sig_plot)[13:16]<-c('Cycles (W1) ~ Variable (W1)','Cycles (W3) ~ Variable (W3)','Cycles (W3) ~ Variable (W1) + Cycles(W1)','Cycles (W1+W3) ~ Variable (W1+W3)')

items_ind_ez2_rev<-as.data.frame(items_ind_ez2_rev)
rownames(effect_size_phen_sig_plot)<-colnames(items_ind_ez2_rev)# was as.character(mb_samp_sp_name[which(pval_chk_all>5)])
rownames(effect_size_phen_sig_plot)[1:2]<-c("Household Wealth Index","MCA coordinate")
rownames(effect_size_phen_sig_plot)[3:dim(items_ind_ez_rev)[2]]<-colnames(items_ind_ez_rev)[3:dim(items_ind_ez_rev)[2]]
rownames(effect_size_phen_sig_plot)[8]<-"No phone"

myb2<-c(0.0000000125,0.0012500000,0.0065105438,0.0112500000,0.0143750000,0.0179000000,0.0342018138,0.0720843541,0.15,1.5)#
myBreaks<-c(myb2[c(length(myb2):-1:1)]*-1,myb2)
myColor<-colorRampPalette(c("#cb181d","white", "#33a02c"))(19)


pdf('regression_items_rev_bl.pdf',width=9,height=7)
pheatmap(effect_size_phen_sig_plot/3,annotation_row = NULL,annotation_names_row = T,
         labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,
         border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,
         color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr_rev),number_color = "black",treeheight_col = 0,
         treeheight_row = 0,cluster_rows = F,cluster_cols = F)
dev.off()

ind_phen<-c(4,8,12,3,7,11)
ind_reg<-c(1,5:7,13:16)

rowr<-rownames(effect_size_phen_sig_plot)[ind_reg]
rowr[5:8]<-c("Toilet type","Cooking stove type","Floor type","Roof type")

pdf('figure_path/figure_1b.pdf',width=5,height=6)#was ind_phen[4:6]
pheatmap(effect_size_phen_sig_plot[ind_reg,ind_phen]/3,annotation_row = NULL,annotation_names_row = T,
         labels_row = rowr,labels_col = colnames(effect_size_phen_sig_plot)[ind_phen],legend = T,fontsize_number = 15,
         border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,
         color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr_rev[ind_phen,ind_reg]),number_color = "black",treeheight_col = 0,
         treeheight_row = 0,cluster_rows = F,cluster_cols = F)
dev.off()


##Cycle composition
kcy_cyc_list<-c("kcy3_wealth_avg","kcy4_wealth_avg","kcy5_wealth_avg")

ind_w13<-intersect(network_factors_vil_w1$Node_id,network_factors_vil_w2$Node_id)
cyc_comp_w13<-network_factors_vil_w1[match(ind_w13,rownames(network_factors_vil_w1)),
                                     colnames(network_factors_vil_w2)[match(kcy_cyc_list,colnames(network_factors_vil_w2))]]
cyc_comp_w13<-cbind(cyc_comp_w13,"delta_W"=array(NA,dim=c(dim(cyc_comp_w13)[1],1)))

temp<-cbind(network_factors_vil_w1$HHW[match(rownames(cyc_comp_w13),network_factors_vil_w1$Node_id)],
            zm4_W2$HHW[match(rownames(cyc_comp_w13),rownames(network_factors_vil_w2))])
temp2<-array(NA,dim=dim(temp)[1])

for(i in 1:length(temp2)){
  if((is.na(temp[i,1])==F)&(is.na(temp[i,2])==F)){
    temp2[i]<-temp[i,2]-temp[i,1]
  }
}
cyc_comp_w13$delta_W<-temp2

cyc_comp_w13<-cbind(cyc_comp_w13,network_factors_vil_w1[match(rownames(cyc_comp_w13),rownames(network_factors_vil_w1)),match("HHW",colnames(network_factors_vil_w1))])
colnames(cyc_comp_w13)[dim(cyc_comp_w13)[2]]<-"Ego_HHW"
cyc_comp_list2<-c("kcycle3","kcycle4","kcycle5")
cyc_comp_w13<-cbind(cyc_comp_w13,network_factors_vil_w1[match(rownames(cyc_comp_w13),rownames(network_factors_vil_w1)),match(cyc_comp_list2,colnames(network_factors_vil_w1))])


cyc_comp_w13_l<-cyc_comp_w13%>%
  pivot_longer(!c(delta_W,Ego_HHW,shape),names_to='Phenotype',values_to='value')

cyc_comp_w13_l2<-cyc_comp_w13_l[which(complete.cases(cyc_comp_w13_l)==T),]

cyc_comp_w13_l2$delta_W_cond<-ifelse(cyc_comp_w13_l2$delta_W>0,"delta_HHW>0",ifelse(cyc_comp_w13_l2$delta_W==0,"delta_HHW=0","delta_HHW<0"))

length(which(cyc_comp_w13_l2$delta_W_cond=="delta_HHW=0"))
length(which(cyc_comp_w13_l2$delta_W_cond=="delta_HHW>0"))
length(which(cyc_comp_w13_l2$delta_W_cond=="delta_HHW<0"))

length(which(cyc_comp_w13_l2$Ego_HHW%in%2:4))


#Relative phenotypes


pv_del_com_ind<-array(NA,dim=c(11,length(unique(cyc_comp_w13_l2$phen))))
rownames(pv_del_com_ind)<-c("ego_1_>0_=0","ego_2_>0_=0","ego_2_>0_<0","ego_2_<0_=0","ego_3_>0_=0","ego_3_>0_<0","ego_3_<0_=0","ego_4_>0_=0","ego_4_>0_<0","ego_4_<0_=0","ego_5_<0_=0")
colnames(pv_del_com_ind)<-unique(cyc_comp_w13_l2$phen)
pv_del_com_ind<-as.data.frame(pv_del_com_ind)



for(i in 1:length(unique(cyc_comp_w13_l2$phen))){
  
  name<-unique(cyc_comp_w13_l2$phen)[i]
  cyc_comp_w13_l3<-cyc_comp_w13_l2[which(cyc_comp_w13_l2$phen%in%name),]
  cyc_comp_w13_l3$delta_W_cond<-factor(cyc_comp_w13_l3$delta_W_cond,levels=c("delta_HHW<0","delta_HHW=0","delta_HHW>0"))
  
  library(ggplot2)
  
  #pdf('comp_delta_w_kcy3.pdf',width=11,height=8)
  #pdf(paste0('comp_del_',name,'.pdf'),width=11,height=8)
  ggplot(cyc_comp_w13_l3,aes(x=as.character(delta_W_cond), y=as.numeric(as.character(value)),fill=as.character(Ego_HHW))) +
    geom_jitter(size=1, alpha=0.5) +
    geom_violin(alpha=0.7) +
    theme(text = element_text(size=17),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5),plot.title = element_text(hjust = 0.5))+
    #ggtitle("Average wealth in cycles") +
    xlab("HHW (W3) - HHW(W1)")+ylab(name)+
    scale_fill_manual(values=c("#1a9641", "#a6d96a", "#ffffbf","#fdae61","#d7191c"))+
    guides(fill=guide_legend(title="HHW (Wave 1)"))+#+ylim(0.8,6.5)
    theme(plot.margin=margin(c(3,0.2,0.2,0.2), unit = "cm"))
  ggsave(paste0('figure_path/comp_del_ind_',name,'.pdf'),width=11,height=8)
  #dev.off()
  
  pv_del_com_ind[1,i]<-wilcox.test(as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW=0")&(cyc_comp_w13_l3$Ego_HHW==1))])),as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW>0")&(cyc_comp_w13_l3$Ego_HHW==1))])))$p.value
  pv_del_com_ind[2,i]<-wilcox.test(as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW=0")&(cyc_comp_w13_l3$Ego_HHW==2))])),as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW>0")&(cyc_comp_w13_l3$Ego_HHW==2))])))$p.value
  pv_del_com_ind[3,i]<-wilcox.test(as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW<0")&(cyc_comp_w13_l3$Ego_HHW==2))])),as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW>0")&(cyc_comp_w13_l3$Ego_HHW==2))])))$p.value
  pv_del_com_ind[4,i]<-wilcox.test(as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW<0")&(cyc_comp_w13_l3$Ego_HHW==2))])),as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW=0")&(cyc_comp_w13_l3$Ego_HHW==2))])))$p.value
  pv_del_com_ind[5,i]<-wilcox.test(as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW=0")&(cyc_comp_w13_l3$Ego_HHW==3))])),as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW>0")&(cyc_comp_w13_l3$Ego_HHW==3))])))$p.value
  pv_del_com_ind[6,i]<-wilcox.test(as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW<0")&(cyc_comp_w13_l3$Ego_HHW==3))])),as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW>0")&(cyc_comp_w13_l3$Ego_HHW==3))])))$p.value
  pv_del_com_ind[7,i]<-wilcox.test(as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW<0")&(cyc_comp_w13_l3$Ego_HHW==3))])),as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW=0")&(cyc_comp_w13_l3$Ego_HHW==3))])))$p.value
  pv_del_com_ind[8,i]<-wilcox.test(as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW=0")&(cyc_comp_w13_l3$Ego_HHW==4))])),as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW>0")&(cyc_comp_w13_l3$Ego_HHW==4))])))$p.value
  pv_del_com_ind[9,i]<-wilcox.test(as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW<0")&(cyc_comp_w13_l3$Ego_HHW==4))])),as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW>0")&(cyc_comp_w13_l3$Ego_HHW==4))])))$p.value
  pv_del_com_ind[10,i]<-wilcox.test(as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW<0")&(cyc_comp_w13_l3$Ego_HHW==4))])),as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW=0")&(cyc_comp_w13_l3$Ego_HHW==4))])))$p.value
  pv_del_com_ind[11,i]<-wilcox.test(as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW=0")&(cyc_comp_w13_l3$Ego_HHW==5))])),as.numeric(as.character(cyc_comp_w13_l3$value[which((cyc_comp_w13_l3$delta_W_cond=="delta_HHW<0")&(cyc_comp_w13_l3$Ego_HHW==5))])))$p.value
  
  
}

pv_del_com_ind<-pv_del_com_ind[,c("kcy3_wealth_avg","kcy4_wealth_avg","kcy5_wealth_avg",
                                  "kcycle3","kcycle4","kcycle5")]

pv_del_com_ind$Comparison<-rownames(pv_del_com_ind)
pv_del_com_ind$Ego_wealth<-c(1:5)
pv_del_com_ind$shape<-23

pv_del_com_ind_long<-pv_del_com_ind%>%
  pivot_longer(!c(Comparison,Ego_wealth,shape),names_to='Phenotype',values_to='value')

pv_del_com_ind_long$Phenotype<-factor(pv_del_com_ind_long$Phenotype,levels=c("kcy3_wealth_avg","kcy4_wealth_avg","kcy5_wealth_avg",
                                                                             "kcycle3","kcycle4","kcycle5"))
pv_del_com_ind_long$value<--log10(as.numeric(as.character(pv_del_com_ind_long$value)))
pv_del_com_ind_long$Ego_wealth<-ifelse(as.numeric(as.character(pv_del_com_ind_long$value))<(-log10(0.05)),0,pv_del_com_ind_long$Ego_wealth)

ggplot(pv_del_com_ind_long,aes(x=Phenotype,y=as.numeric(as.character(value)),fill=as.factor(Ego_wealth),shape=as.factor(shape)))+#,shape=as.factor(Comparison)
  geom_point(size=3,cex=1,pch=as.numeric(as.character(pv_del_com_ind_long$shape)),color="black")+ylab("-log10(p-value)")+
  scale_fill_manual(values=c("#BDBDBD","#1a9641", "#a6d96a", "#ffffbf","#fdae61","#d7191c"))+
  theme(text = element_text(size=15),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black",angle=60,hjust=1),
        panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept=(-log10(0.05)),linetype='dashed')+
  scale_x_discrete(labels=rep(c("Kcycle-3","Kcycle-4","Kcycle-5"),4))+
  annotate("text",x=2.5,y=94,label="Number of cycles")+#ylim(-1,60)+
  annotate("text",x=6.5,y=92,label="Alter's wealth average")+#ylim(-1,60)+
  annotate("text",x=10.5,y=94,label="Alter's weighted wealth average")+ylim(-1,20)+
  annotate("text",x=14.5,y=92,label="Gini coefficient")+#ylim(-1,35)+
  geom_rect(xmin=0.8,xmax=4.2,ymin=88,ymax=89,fill="#1f78b4")+
  geom_rect(xmin=4.8,xmax=8.2,ymin=88,ymax=89,fill="#cab2d6")+
  geom_rect(xmin=8.8,xmax=12.2,ymin=88,ymax=89,fill="#fb8072")+
  geom_rect(xmin=12.8,xmax=16.2,ymin=88,ymax=89,fill="#b3de69")+
  ggtitle("Comparisons (Individual cycles)")


ggsave('figure_path/fig_3d_manhattan.pdf',width=7,height=5)











