library(tidyverse)
library(dplyr)
library(igraph)
library(intergraph)
library(tidygraph)
library(ggraph)
library(sna)

library(ggplot2)

#------------------------ NOTES ---------------------------
# The aim of this file is to investigate the observed cycles in undirected friend vs
# borrow-lend networks where ties in the same household have been removed from the cycles analysis.

# This file calculates the observed number of cycles using data source 1 and 2,
# meaning the data includes individuals who completed the census, but not the w1 
# survey. 6% of the total household population doesn't have a household wealth index
# and are from data source 2 households, as opposed to only 1% of the household
# population not having a household wealth index from data source 1 households.

# NAME GENERATOR QUESTIONS:
# 1. Who would you feel comfortable asking to borrow 200 lempiras from if you needed
#    them for the day?
# 2. Who do you think would be comfortable asking you to borrow 200 lempiras for the day?


#--------------------- loading data ------------------------
setwd("/scratch")

# Create a data path for saving cleaned / augmented data sets
data_path <- "/scratch"

# Respondents data (original to have full village networks)
load(paste(c(data_path, "/original_netwrks_resp_data_*.rda"), collapse=""))

# Subset data
conn_w1 <- subset(conn_w1_original, alter_source==1 & village_code_w1!=0 & same_village==1)
conn_w3 <- subset(conn_w3_original, alter_source==1 & village_code_w3!=0 & same_village==1)

#--------------- 1: Friend indiv networks (wave 1) ----------------
# Create a node list of all unique egos and alters. Create unique integer id's for all respondent_master_id's.
# node_attributes_w1 <- resp_w1_original %>%
#   subset(data_source_w1==1 & village_code_w1!=0) %>%
#   select(respondent_master_id, building_id_w1, gender, age_at_survey) %>%
#   mutate(female = ifelse(gender=="female", 1, 0)) %>%
#   left_join(hh_w1 %>%
#               select(building_id, household_wealth_index_w1),
#             by = c("building_id_w1"="building_id"))
# 
nodelist_w1 <- conn_w1 %>%
  subset(select=c(ego, village_code_w1)) %>%
  dplyr::rename(node = ego) %>%
  bind_rows(conn_w1 %>%
              subset(select=c(alter, village_code_w1)) %>%
              dplyr::rename(node = alter)) %>%
  distinct(node, .keep_all=TRUE) %>%
  rowid_to_column("id")

friend_diff_hh_edgelist_w1 <- conn_w1 %>%
  subset(same_village==1 & same_building==0 &
           relationship %in% c("closest_friend", "personal_private", "free_time"),
         select = c(ego, alter, relationship, village_code_w1)) %>%
  # bind_rows(conn_w1 %>%
  #             subset(same_village==1 & relationship %in% c("trust_lend_money"),
  #                    select=c(alter, ego, relationship, village_code_w1)) %>%
  #             rename(ego = alter,
  #                    alter = ego)) %>%
  left_join(nodelist_w1 %>%
              select(id, node), by = c("ego" = "node")) %>%
  dplyr::rename(from = id) %>%
  inner_join(nodelist_w1 %>%
               select(id, node), by = c("alter" = "node")) %>%
  dplyr::rename(to = id) %>%
  select(from, to, village_code_w1) %>%
  distinct(from, to, .keep_all = TRUE) %>%
  subset(from != to)

# Pull network attributes
# Create list to hold village-level network statistics for each village
vill_net_attr_ls <- vector(mode = "list", length = 176)

for (i in 1:176) {
  edgelist_subset <- subset(friend_diff_hh_edgelist_w1, village_code_w1==i)
  n <- graph_from_data_frame(edgelist_subset, directed=FALSE)
  n <- simplify(n)
  isolated <- which(degree(asNetwork(n))==0)
  n <- igraph::delete.vertices(n, isolated)
  
  vill_net_attr_ls[[i]] <- data.frame(village_code = i,
                                      wave = "w1") %>%
    as_tibble %>%
    mutate(num_nodes = vcount(n),
           num_edges = ecount(n),
           diameter = diameter(n, directed = FALSE, unconnected = TRUE, weights = NULL),
           distance = mean_distance(n),
           ave_degree = num_edges / num_nodes,
           density = edge_density(n, loops=FALSE),
           transitivity = transitivity(n, type="global"),
           reciprocity = reciprocity(n, ignore.loops = TRUE))
}

vill_net_attr_friend_diff_hh_w1_df <- do.call(rbind, vill_net_attr_ls)

# Initialize a list to hold kcycle stats (village_level)
kcycle_list.friend <- vector(mode = "list", length = 176)

for (i in 1:176) {
  edgelist_subset <- subset(friend_diff_hh_edgelist_w1, village_code_w1==i)
  n <- graph_from_data_frame(edgelist_subset, directed=FALSE)
  n <- simplify(n)
  n <- asNetwork(n)
  
  kcycle_list.friend[[i]] <- kcycle.census(n, maxlen = 5, mode = "graph", 
                                                tabulate.by.vertex = FALSE, 
                                                cycle.comembership = "none")
  print(i)
}

# Initialize a list to hold kpath stats (village_level)
kpath_list.friend <- vector(mode = "list", length = 176)

for (i in 1:176) {
  edgelist_subset <- subset(friend_diff_hh_edgelist_w1, village_code_w1==i)
  n <- graph_from_data_frame(edgelist_subset, directed=FALSE)
  n <- simplify(n)
  n <- asNetwork(n)
  
  kpath_list.friend[[i]] <- kpath.census(n, maxlen = 5, mode = "graph", 
                                              tabulate.by.vertex = FALSE, 
                                              path.comembership = "none", 
                                              dyadic.tabulation = "none")
  print(i)
}

kstats_friend_diff_hh_w1_df <- do.call(rbind, unlist(kcycle_list.friend, recursive=FALSE)) %>%
  as_tibble() %>%
  rowid_to_column("village_code") %>%
  remove_rownames() %>%
  dplyr::rename(kcycle2 = `2`,
                kcycle3 = `3`,
                kcycle4 = `4`,
                kcycle5 = `5`) %>%
  cbind(relationship = 'borrow') %>%
  left_join(do.call(rbind, unlist(kpath_list.friend, recursive = FALSE)) %>%
              as_tibble() %>%
              rowid_to_column("village_code") %>%
              remove_rownames() %>%
              dplyr::rename(kpath1 = `1`,
                            kpath2 = `2`,
                            kpath3 = `3`,
                            kpath4 = `4`,
                            kpath5 = `5`))

#-------------- 2: Friend indiv networks (wave 3) -------------------

all_vill_w3_vec <- c(1:176)
vill_w3_vec <- all_vill_w3_vec

# Use the same node list as indiv networks using borrow and lend ties (nodelist_w3)
nodelist_w3 <- conn_w3 %>%
  subset(select=c(ego, village_code_w3)) %>%
  dplyr::rename(node = ego) %>%
  bind_rows(conn_w3 %>%
              subset(select=c(alter, village_code_w3)) %>%
              dplyr::rename(node = alter)) %>%
  distinct(node, .keep_all=TRUE) %>%
  rowid_to_column("id") %>%
  subset(village_code_w3 %in% vill_w3_vec)

friend_diff_hh_edgelist_w3 <- conn_w3 %>%
  subset(same_village==1 & same_building==0 &
           relationship %in% c("closest_friend", "personal_private", "free_time"),
         select = c(ego, alter, relationship, village_code_w3)) %>%
  left_join(nodelist_w3 %>%
              select(id, node), by = c("ego" = "node")) %>%
  dplyr::rename(from = id) %>%
  inner_join(nodelist_w3 %>%
               select(id, node), by = c("alter" = "node")) %>%
  dplyr::rename(to = id) %>%
  select(from, to, village_code_w3) %>%
  distinct(from, to, .keep_all = TRUE) %>%
  subset(from != to & village_code_w3 %in% vill_w3_vec)

# Pull network attributes
# Create list to hold village-level network statistics for each village
vill_net_attr_ls <- vector(mode = "list", length = length(vill_w3_vec))

for (i in 1:length(vill_w3_vec)) {
  vill_code <- vill_w3_vec[i]
  edgelist_subset <- subset(friend_diff_hh_edgelist_w3, village_code_w3==vill_code)
  n <- graph_from_data_frame(edgelist_subset, directed=FALSE)
  n <- simplify(n)
  
  vill_net_attr_ls[[i]] <- data.frame(village_code = vill_code,
                                      wave = "w3") %>%
    as_tibble %>%
    mutate(num_nodes = vcount(n),
           num_edges = ecount(n),
           diameter = diameter(n, directed = FALSE, unconnected = TRUE, weights = NULL),
           distance = mean_distance(n),
           ave_degree = num_edges / num_nodes,
           density = edge_density(n, loops=FALSE),
           transitivity = transitivity(n, type="global"),
           reciprocity = reciprocity(n, ignore.loops = TRUE))
}

vill_net_attr_friend_diff_hh_w3_df <- do.call(rbind, vill_net_attr_ls)

# Initialize a list to hold kcycle stats (village_level)
kcycle_list.friend <- vector(mode = "list", length = length(vill_w3_vec))

for (i in 1:length(vill_w3_vec)) {
  vill_code <- vill_w3_vec[i]
  edgelist_subset <- subset(friend_diff_hh_edgelist_w3, village_code_w3==vill_code)
  n <- graph_from_data_frame(edgelist_subset, directed=FALSE)
  n <- simplify(n)
  n <- asNetwork(n)

  kcycle_list.friend[[i]] <- kcycle.census(n, maxlen = 5, mode = "graph", 
                                                tabulate.by.vertex = FALSE, 
                                                cycle.comembership = "none")
  print(i)
}

# Initialize a list to hold kpath stats (village_level)
kpath_list.friend <- vector(mode = "list", length = length(vill_w3_vec))

for (i in 1:length(vill_w3_vec)) {
  vill_code <- vill_w3_vec[i]
  edgelist_subset <- subset(friend_diff_hh_edgelist_w3, village_code_w3==vill_code)
  n <- graph_from_data_frame(edgelist_subset, directed=FALSE)
  n <- simplify(n)
  n <- asNetwork(n)
  
  # CHECK!!!! this should actually be kpath_list.friend[[i]]
  kpath_list.friend[[i]] <- kpath.census(n, maxlen = 5, mode = "graph", 
                                              tabulate.by.vertex = FALSE, 
                                              path.comembership = "none", 
                                              dyadic.tabulation = "none")
  print(i)
}

kstats_friend_diff_hh_w3_df <- do.call(rbind, unlist(kcycle_list.friend, recursive=FALSE)) %>%
  as_tibble() %>%
  mutate(village_code = vill_w3_vec) %>%
  remove_rownames() %>%
  dplyr::rename(kcycle2 = `2`,
                kcycle3 = `3`,
                kcycle4 = `4`,
                kcycle5 = `5`) %>%
  cbind(relationship = 'borrow') %>%
  left_join(do.call(rbind, unlist(kpath_list.friend, recursive = FALSE)) %>%
              as_tibble() %>%
              mutate(village_code = vill_w3_vec) %>%
              remove_rownames() %>%
              dplyr::rename(kpath1 = `1`,
                            kpath2 = `2`,
                            kpath3 = `3`,
                            kpath4 = `4`,
                            kpath5 = `5`))




#------------------------ saving data ---------------------------
# Save data frame of network attributes
save(vill_net_attr_friend_diff_hh_w1_df, vill_net_attr_friend_diff_hh_w3_df,
     file = file.path(data_path, "vill_net_attr_friend_diff_hh_w1w3_df_*.rda"))

# Save data frames of network features.
save(kstats_friend_diff_hh_w1_df, kstats_friend_diff_hh_w3_df, 
     file = file.path(data_path, "net_attr_friend_ds2_w1w3_*.rda"))

# Save hh borrow edge list and node list.
save(nodelist_w1, friend_diff_hh_edgelist_w1, nodelist_w3, friend_diff_hh_edgelist_w3,
     file = file.path(data_path, "friend_diff_hh_edge_node_list_w1w3_*.rda"))

