library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)

library(labelled)


#--------------------- loading data ------------------------
#setwd("/WORKAREA/work/HONDURAS_GATES/BORROWING")
setwd("C:/Users/Shiv/Documents/Honduras/Kcycle and networks 04-19-23/kcycle_new_s/00_selena_final_code")
# Create a data path for saving cleaned / augmented data sets
data_path <- "00_data"
borrow_data_path <- paste0(data_path, "/mccleary_output/borrow_lend")
friend_data_path <- paste0(data_path, "/mccleary_output/friend")
figure_path <- paste0(data_path, "/figures")

# Load raw expected cycles count (same degree distribution)
load(paste0(borrow_data_path, "/exp_cyc_bl_diff_hh_indiv_w1_20240728.rda"))
load(paste0(borrow_data_path, "/exp_cyc_bl_diff_hh_indiv_w3_20240728.rda"))
load(paste0(friend_data_path, "/exp_cyc_indiv_friend_diff_hh_w1_20240728.rda"))
load(paste0(friend_data_path, "/exp_cyc_indiv_friend_diff_hh_w3_20240728.rda"))

# Load village-level data
load(paste(c(data_path, "/vill_df_20230911.rda"), collapse=""))
load(paste(c(data_path, "/gini_df_20240327.rda"), collapse=""))
load(paste(c(data_path, "/wealth_measures_df_20240327.rda"), collapse=""))

# Load network attributes data
load(paste0(borrow_data_path, "/net_attr_borrow_lend_undir_ds2_w1w3_20240728.rda"))
load(paste0(friend_data_path, "/net_attr_friend_ds2_w1w3_20240728.rda"))

# # Load cyclic and path density measure for correlations (individual)
# load(paste(c(data_path, "/03_merged_data/analysis_20240613.rda"), collapse=""))


#---------- distribution of raw cycles ----------------
### Compare raw distributions of cycles (randomly generated vs observed)
# Convert matrix of # of cycles in randomly generated networks to single distribution
# Do not convert for kcycle 2 (no cycles for undirected networks)
kcycle3_fr <- c(c(as.matrix(vil_indiv_friend_w3_diff_hh_kcy3)), c(as.matrix(vil_indiv_friend_w3_diff_hh_kcy3)))
kcycle4_fr <- c(c(as.matrix(vil_indiv_friend_w3_diff_hh_kcy4)), c(as.matrix(vil_indiv_friend_w3_diff_hh_kcy4)))
kcycle5_fr <- c(c(as.matrix(vil_indiv_friend_w3_diff_hh_kcy5)), c(as.matrix(vil_indiv_friend_w3_diff_hh_kcy5)))

kcycle3_bl <- c(c(as.matrix(vil_indiv_w1_bl_diff_hh_kcy3)), c(as.matrix(vil_indiv_w3_bl_diff_hh_kcy3)))
kcycle4_bl <- c(c(as.matrix(vil_indiv_w1_bl_diff_hh_kcy4)), c(as.matrix(vil_indiv_w3_bl_diff_hh_kcy4)))
kcycle5_bl <- c(c(as.matrix(vil_indiv_w1_bl_diff_hh_kcy5)), c(as.matrix(vil_indiv_w3_bl_diff_hh_kcy5)))

vil_all_sds_fr <- data.frame(kcycle3_fr, kcycle4_fr, kcycle5_fr)
vil_all_sds_bl <- data.frame(kcycle3_bl, kcycle4_bl, kcycle5_bl)

vil_all_fr_w1w3 <- kstats_friend_diff_hh_w1_df %>%
  select(kcycle3, kcycle4, kcycle5) %>%
  bind_rows(kstats_friend_diff_hh_w3_df %>%
              select(kcycle3, kcycle4, kcycle5))

vil_all_bl_w1w3 <- kstats_borrow_lend_diff_hh_undir_w1_df %>%
  select(kcycle3, kcycle4, kcycle5) %>%
  bind_rows(kstats_borrow_lend_diff_hh_undir_w3_df %>%
              select(kcycle3, kcycle4, kcycle5))

raw_cyc_distr_fr_df <- bind_rows(list(observed = vil_all_fr_w1w3, sds = vil_all_sds_fr), .id="network_type")
raw_cyc_distr_bl_df <- bind_rows(list(observed = vil_all_bl_w1w3, sds = vil_all_sds_bl), .id="network_type")

# Calculate the z-score for the number of raw cycles observed in each village as
# compared with the expected number of cycles.
# z-score = (obs - mean) / sd
# For friend networks (w1):
vc <- c(1:176)
vil_all <- kstats_friend_diff_hh_w1_df %>%
  select(kcycle3, kcycle4, kcycle5)

cyc_z_score<-array(NA,dim=c(length(vc),3))
colnames(cyc_z_score)<-c("kcycle3","kcycle4","kcycle5")

for(j in 1:length(vc)){
  mean_cyc3 <- mean(vil_indiv_friend_w1_diff_hh_kcy3[j,])
  sd_cyc3 <- sd(vil_indiv_friend_w1_diff_hh_kcy3[j,])
  cyc_z_score[j,1] <- (vil_all[j,1] - mean_cyc3) / sd_cyc3
  
  mean_cyc4 <- mean(vil_indiv_friend_w1_diff_hh_kcy4[j,])
  sd_cyc4 <- sd(vil_indiv_friend_w1_diff_hh_kcy4[j,])
  cyc_z_score[j,2] <- (vil_all[j,2] - mean_cyc4) / sd_cyc4
  
  mean_cyc5 <- mean(vil_indiv_friend_w1_diff_hh_kcy5[j,])
  sd_cyc5 <- sd(vil_indiv_friend_w1_diff_hh_kcy5[j,])
  cyc_z_score[j,3] <- (vil_all[j,3] - mean_cyc5) / sd_cyc5
}

sds_cyc_distr_fr_w1 <- as.data.frame(cyc_z_score) %>%
  bind_cols(vill_df %>%
              select(vill_size))

# For friend networks (w3):
all_vill_w3_vec <- c(1:176)
vill_w3_vec <- all_vill_w3_vec[!all_vill_w3_vec %in% c(155, 156)]

vil_all <- kstats_friend_diff_hh_w3_df %>%
  select(kcycle3, kcycle4, kcycle5) %>%
  drop_na()

cyc_z_score<-array(NA,dim=c(length(vill_w3_vec),3))
colnames(cyc_z_score)<-c("kcycle3","kcycle4","kcycle5")

for(j in 1:length(vill_w3_vec)){
  mean_cyc3 <- mean(vil_indiv_friend_w3_diff_hh_kcy3[j,])
  sd_cyc3 <- sd(vil_indiv_friend_w3_diff_hh_kcy3[j,])
  cyc_z_score[j,1] <- (vil_all[j,1] - mean_cyc3) / sd_cyc3
  
  mean_cyc4 <- mean(vil_indiv_friend_w3_diff_hh_kcy4[j,])
  sd_cyc4 <- sd(vil_indiv_friend_w3_diff_hh_kcy4[j,])
  cyc_z_score[j,2] <- (vil_all[j,2] - mean_cyc4) / sd_cyc4
  
  mean_cyc5 <- mean(vil_indiv_friend_w3_diff_hh_kcy5[j,])
  sd_cyc5 <- sd(vil_indiv_friend_w3_diff_hh_kcy5[j,])
  cyc_z_score[j,3] <- (vil_all[j,3] - mean_cyc5) / sd_cyc5
}

sds_cyc_distr_fr_w3 <- as.data.frame(cyc_z_score) %>%
  bind_cols(vill_df %>%
              subset(!(village_code %in% c(155, 156))) %>%
              select(vill_size))

# For borrow-lend networks (w1):
vc <- c(1:176)
vil_all <- kstats_borrow_lend_diff_hh_undir_w1_df %>%
  select(kcycle3, kcycle4, kcycle5) %>%
  drop_na()

cyc_z_score<-array(NA,dim=c(length(vc),3))
colnames(cyc_z_score)<-c("kcycle3","kcycle4","kcycle5")

for(j in 1:length(vc)){
  mean_cyc3 <- mean(vil_indiv_w1_bl_diff_hh_kcy3[j,])
  sd_cyc3 <- sd(vil_indiv_w1_bl_diff_hh_kcy3[j,])
  cyc_z_score[j,1] <- (vil_all[j,1] - mean_cyc3) / sd_cyc3
  
  mean_cyc4 <- mean(vil_indiv_w1_bl_diff_hh_kcy4[j,])
  sd_cyc4 <- sd(vil_indiv_w1_bl_diff_hh_kcy4[j,])
  cyc_z_score[j,2] <- (vil_all[j,2] - mean_cyc4) / sd_cyc4
  
  mean_cyc5 <- mean(vil_indiv_w1_bl_diff_hh_kcy5[j,])
  sd_cyc5 <- sd(vil_indiv_w1_bl_diff_hh_kcy5[j,])
  cyc_z_score[j,3] <- (vil_all[j,3] - mean_cyc5) / sd_cyc5
}

sds_cyc_distr_bl_w1 <- as.data.frame(cyc_z_score) %>%
  bind_cols(vill_df %>%
              select(vill_size))

# For borrow-lend networks (w3):
all_vill_w3_vec <- c(1:176)
vill_w3_vec <- all_vill_w3_vec[!all_vill_w3_vec %in% c(155, 156)]

vil_all <- kstats_borrow_lend_diff_hh_undir_w3_df %>%
  select(kcycle3, kcycle4, kcycle5) %>%
  drop_na()

vil_all_w1 <- kstats_borrow_lend_diff_hh_undir_w1_df %>%
  select(kcycle3, kcycle4, kcycle5) %>%
  drop_na()

cyc_z_score<-array(NA,dim=c(length(vill_w3_vec),3))
colnames(cyc_z_score)<-c("kcycle3","kcycle4","kcycle5")

for(j in 1:length(vill_w3_vec)){
  mean_cyc3 <- mean(vil_indiv_w3_bl_diff_hh_kcy3[j,])
  sd_cyc3 <- sd(vil_indiv_w3_bl_diff_hh_kcy3[j,])
  cyc_z_score[j,1] <- (vil_all[j,1] - mean_cyc3) / sd_cyc3
  
  mean_cyc4 <- mean(vil_indiv_w3_bl_diff_hh_kcy4[j,])
  sd_cyc4 <- sd(vil_indiv_w3_bl_diff_hh_kcy4[j,])
  cyc_z_score[j,2] <- (vil_all[j,2] - mean_cyc4) / sd_cyc4#Here
  
  mean_cyc5 <- mean(vil_indiv_w3_bl_diff_hh_kcy5[j,])
  sd_cyc5 <- sd(vil_indiv_w3_bl_diff_hh_kcy5[j,])
  cyc_z_score[j,3] <- (vil_all[j,3] - mean_cyc5) / sd_cyc5
}

sds_cyc_distr_bl_w3 <- as.data.frame(cyc_z_score) %>%
  bind_cols(vill_df %>%
              subset(!(village_code %in% c(155, 156))) %>%#Missing villages
              select(vill_size))


# Merge z-scores across waves and name-generator types for visualization
cyc_distr_sds_df <- sds_cyc_distr_bl_w1 %>%
  mutate(wave = 1,
         name_generator = "borrow_lend") %>%
  bind_rows(sds_cyc_distr_bl_w3 %>%
              mutate(wave = 3,
                     name_generator = "borrow_lend")) %>%
  bind_rows(sds_cyc_distr_fr_w1 %>%
              mutate(wave = 1,
                     name_generator = "friend")) %>%
  bind_rows(sds_cyc_distr_fr_w3 %>%
              mutate(wave = 3,
                     name_generator = "friend"))

cyc_distr_sds_long_df <- cyc_distr_sds_df %>%
  pivot_longer(cols = -c(wave, vill_size, name_generator),
               names_to = "cyc_length",
               values_to = "cyc_z_score") %>%
  mutate(cyc_length = case_when(cyc_length=="kcycle3" ~ "3",
                                cyc_length=="kcycle4" ~ "4",
                                cyc_length=="kcycle5" ~ "5"))

cyc_z_score_bl_long_df <- subset(cyc_distr_sds_long_df, name_generator=="borrow_lend")
cyc_z_score_fr_long_df <- subset(cyc_distr_sds_long_df, name_generator=="friend")

### Density distributions of cycle count z-scores:
network_names <- as_labeller(
  c(borrow_lend = "Borrow-lend networks", friend = "Friend networks"))

ggplot(data = cyc_distr_sds_long_df,
       aes(x = cyc_z_score, 
           group = cyc_length, 
           fill = cyc_length)) + 
  geom_vline(xintercept = 1.96, color = "red") +
  geom_density(adjust=1.5, alpha=.4) +
  labs(x = "Normalized observed vs expected cycles (Z-scores)",#title = "Observed to expected K-cycle",
       fill = "Cycle length") +
  facet_wrap(~name_generator, labeller=network_names) +
  theme_classic()
ggsave("cyc_zscore_distr.png", width = 20, height = 10, units = "cm")

#Village statistics -- for MS -- 
#BL
bl_w1_5<-cyc_distr_sds_long_df[which((cyc_distr_sds_long_df$name_generator=="borrow_lend")&
                                                       (cyc_distr_sds_long_df$wave==1)&
                                                       (cyc_distr_sds_long_df$cyc_length==5)),]
bl_w3_5<-cyc_distr_sds_long_df[which((cyc_distr_sds_long_df$name_generator=="borrow_lend")&
                                       (cyc_distr_sds_long_df$wave==3)&
                                       (cyc_distr_sds_long_df$cyc_length==5)),]

(length(which(bl_w1_5$cyc_z_score>1.96))+length(which(bl_w3_5$cyc_z_score>1.96)))/(dim(bl_w1_5)[1]+dim(bl_w3_5)[1])*100


#Frnd
fr_w1_5<-cyc_distr_sds_long_df[which((cyc_distr_sds_long_df$name_generator=="friend")&
                                       (cyc_distr_sds_long_df$wave==1)&
                                       (cyc_distr_sds_long_df$cyc_length==5)),]
fr_w3_5<-cyc_distr_sds_long_df[which((cyc_distr_sds_long_df$name_generator=="friend")&
                                       (cyc_distr_sds_long_df$wave==3)&
                                       (cyc_distr_sds_long_df$cyc_length==5)),]

(length(which(fr_w1_5$cyc_z_score>1.96))+length(which(fr_w3_5$cyc_z_score>1.96)))/(dim(fr_w1_5)[1]+dim(fr_w3_5)[1])*100



#ggsave("cyc_zscore_distr.pdf", width = 20, height = 10, units = "cm")

ggplot(data = cyc_distr_sds_long_df,
       aes(x = cyc_z_score, 
           group = cyc_length, 
           fill = cyc_length)) + 
  geom_vline(xintercept = 1.96, color = "red") +
  geom_density(adjust=1.5, alpha=.4) +
  labs(x = "Normalized observed to expected cycles (Z-scores)",#title = "Z-scores for observed compared vs expected cycles",
       fill = "Cycle length",y="Probability density") +
  facet_wrap(~name_generator, labeller=network_names) +
  theme(text = element_text(size=15),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),legend.position="none",plot.title = element_text(hjust = 0.3))
ggsave("cyc_dist_zscore_sel_fig1.pdf",width=10,height=5)

#kcy3 > expected

temp<-cyc_distr_sds_long_df[which(cyc_distr_sds_long_df$name_generator=="borrow_lend"),]
length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==3)))/length(which((temp$cyc_length==3)))
length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==4)))/length(which((temp$cyc_length==4)))
length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==5)))/length(which((temp$cyc_length==5)))

temp<-cyc_distr_sds_long_df[which(cyc_distr_sds_long_df$name_generator=="friend"),]
length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==3)))/length(which((temp$cyc_length==3)))
length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==4)))/length(which((temp$cyc_length==4)))
length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==5)))/length(which((temp$cyc_length==5)))

## Split by wave

cyc_distr_sds_long_df2<-cyc_distr_sds_long_df
cyc_distr_sds_long_df2$name_wave<-paste0(cyc_distr_sds_long_df2$name_generator,"_",cyc_distr_sds_long_df2$wave)

network_waves <- as_labeller(
  c(borrow_lend_1 = "Borrow/lend networks (Wave 1)",borrow_lend_3 = "Borrow/lend networks (Wave 3)", 
    friend_1 = "Friendship networks (Wave 1)",friend_3 = "Friendship networks (Wave 3)"))

ggplot(data = cyc_distr_sds_long_df2,
       aes(x = cyc_z_score, 
           group = cyc_length, 
           fill = cyc_length)) + 
  geom_vline(xintercept = 1.96, color = "red") +
  geom_density(adjust=1.5, alpha=.4) +
  labs(title = "Z-scores for observed compared vs expected cycles",
       x = "Z-score for observed cycles count per village",
       fill = "Cycle length",y="Probability density") +
  facet_wrap(~name_wave, labeller=network_waves,nrow=1) +
  theme(text = element_text(size=15),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),legend.position="none",plot.title = element_text(hjust = 0.3))
ggsave("cyc_dist_zscore_sel_wv.pdf",width=15,height=5)

temp<-cyc_distr_sds_long_df2
length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==3)&(temp$name_wave=="borrow_lend_1")))/length(which((temp$cyc_length==3)&(temp$name_wave=="borrow_lend_1")))
length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==4)&(temp$name_wave=="borrow_lend_1")))/length(which((temp$cyc_length==4)&(temp$name_wave=="borrow_lend_1")))
length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==5)&(temp$name_wave=="borrow_lend_1")))/length(which((temp$cyc_length==5)&(temp$name_wave=="borrow_lend_1")))

length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==3)&(temp$name_wave=="borrow_lend_3")))/length(which((temp$cyc_length==3)&(temp$name_wave=="borrow_lend_3")))
length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==4)&(temp$name_wave=="borrow_lend_3")))/length(which((temp$cyc_length==4)&(temp$name_wave=="borrow_lend_3")))
length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==5)&(temp$name_wave=="borrow_lend_3")))/length(which((temp$cyc_length==5)&(temp$name_wave=="borrow_lend_3")))

length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==3)&(temp$name_wave=="friend_1")))/length(which((temp$cyc_length==3)&(temp$name_wave=="friend_1")))
length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==4)&(temp$name_wave=="friend_1")))/length(which((temp$cyc_length==4)&(temp$name_wave=="friend_1")))
length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==5)&(temp$name_wave=="friend_1")))/length(which((temp$cyc_length==5)&(temp$name_wave=="friend_1")))

length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==3)&(temp$name_wave=="friend_3")))/length(which((temp$cyc_length==3)&(temp$name_wave=="friend_3")))
length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==4)&(temp$name_wave=="friend_3")))/length(which((temp$cyc_length==4)&(temp$name_wave=="friend_3")))
length(which((temp$cyc_z_score>1.96)&(temp$cyc_length==5)&(temp$name_wave=="friend_3")))/length(which((temp$cyc_length==5)&(temp$name_wave=="friend_3")))


#############################################################
# Borrow network only
ggplot(data = cyc_z_score_bl_long_df,
       aes(x = cyc_z_score, 
           group = cyc_length, 
           fill = cyc_length)) + 
  geom_vline(xintercept = 1.96, color = "red") +
  geom_density(adjust=1.5, alpha=.4) +
  labs(title = "Distribution of z-scores for observed compared to expected cycles count",
       x = "Z-score for observed cycles count per village",
       fill = "Cycle length") +
  theme_classic()
ggsave("cyc_zscore_bl_distr.png", path = figure_path, width = 16, height = 10, units = "cm")

### Scatter plots of z-scores and size (indiv)
cyc_distr_zsc_size_cy2 <- ggplot(data = cyc_distr_sds_indiv_df, 
                                 aes(x = vill_size, 
                                     y = kcycle2)) +
  geom_point() +
  labs(title = "Kcycle 2") +
  xlab("Z-score of observed network (kcycle 2)") +
  ylab("Village size") +
  theme_classic()

cyc_distr_zsc_size_cy3 <- ggplot(data = cyc_distr_sds_indiv_df, 
                                 aes(x = vill_size, 
                                     y = kcycle3)) +
  geom_point() +
  labs(title = "Kcycle 3") +
  xlab("Z-score of observed network (kcycle 3)") +
  ylab("Village size") +
  theme_classic()

cyc_distr_zsc_size_cy4 <- ggplot(data = cyc_distr_sds_indiv_df, 
                                 aes(x = vill_size, 
                                     y = kcycle4)) +
  geom_point() +
  labs(title = "Kcycle 4") +
  xlab("Z-score of observed network (kcycle 4)") +
  ylab("Village size") +
  theme_classic()

cyc_distr_zsc_size_cy5 <- ggplot(data = cyc_distr_sds_indiv_df, 
                                 aes(x = vill_size, 
                                     y = kcycle5)) +
  geom_point() +
  labs(title = "Kcycle 5") +
  xlab("Z-score of observed network (kcycle 5)") +
  ylab("Village size") +
  theme_classic()

ggarrange(cyc_distr_zsc_size_cy2, cyc_distr_zsc_size_cy3, cyc_distr_zsc_size_cy4, cyc_distr_zsc_size_cy5,
          ncol = 2, nrow = 2)
ggsave("cyc_distr_sds_villsize.png", path = figure_path, width = 40, height = 25, units = "cm")



# Identify % of villages with higher than expected cycles
sum(sds_cyc_distr_indiv_w1$kcycle2 > 1.96, na.rm=TRUE)
sum(sds_cyc_distr_indiv_w1$kcycle3 > 1.96, na.rm=TRUE)
sum(sds_cyc_distr_indiv_w1$kcycle4 > 1.96, na.rm=TRUE)
sum(sds_cyc_distr_indiv_w1$kcycle5 > 1.96, na.rm=TRUE)

sum(sds_cyc_distr_indiv_w3$kcycle2 > 1.96, na.rm=TRUE)
sum(sds_cyc_distr_indiv_w3$kcycle3 > 1.96, na.rm=TRUE)
sum(sds_cyc_distr_indiv_w3$kcycle4 > 1.96, na.rm=TRUE)
sum(sds_cyc_distr_indiv_w3$kcycle5 > 1.96, na.rm=TRUE)

# Identify % of villages with fewer than expected cycles
sum(sds_cyc_distr_indiv_w1$kcycle2 < -1.96, na.rm=TRUE)
sum(sds_cyc_distr_indiv_w1$kcycle3 < -1.96, na.rm=TRUE)
sum(sds_cyc_distr_indiv_w1$kcycle4 < -1.96, na.rm=TRUE)
sum(sds_cyc_distr_indiv_w1$kcycle5 < -1.96, na.rm=TRUE)

sum(sds_cyc_distr_indiv_w3$kcycle2 < -1.96, na.rm=TRUE)
sum(sds_cyc_distr_indiv_w3$kcycle3 < -1.96, na.rm=TRUE)
sum(sds_cyc_distr_indiv_w3$kcycle4 < -1.96, na.rm=TRUE)
sum(sds_cyc_distr_indiv_w3$kcycle5 < -1.96, na.rm=TRUE)


### Density distributions of raw cycles (indiv)
raw_cyc_distr_cy2 <- ggplot(data = raw_cyc_distr_df, 
                             aes(x = kcycle2, 
                                 group = network_type, 
                                 fill = network_type)) +
  xlim(0, 250) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(title = "Distribution of raw cycles for cycles of length 2") +
  theme_classic()

raw_cyc_distr_cy3 <- ggplot(data = raw_cyc_distr_df, 
                            aes(x = kcycle3, 
                                group = network_type, 
                                fill = network_type)) +
  xlim(0, 250) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(title = "Distribution of raw cycles for cycles of length 3") +
  theme_classic()

raw_cyc_distr_cy4 <- ggplot(data = raw_cyc_distr_df, 
                            aes(x = kcycle4, 
                                group = network_type, 
                                fill = network_type)) +
  xlim(0, 250) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(title = "Distribution of raw cycles for cycles of length 4") +
  theme_classic()

raw_cyc_distr_cy5 <- ggplot(data = raw_cyc_distr_df, 
                            aes(x = kcycle5, 
                                group = network_type, 
                                fill = network_type)) +
  xlim(0, 250) +
  geom_density(adjust=1.5, alpha=.4) +
  labs(title = "Distribution of raw cycles for cycles of length 5") +
  theme_classic()

ggarrange(raw_cyc_distr_cy2, raw_cyc_distr_cy3, raw_cyc_distr_cy4, raw_cyc_distr_cy5,
          ncol = 2, nrow = 2)
ggsave("raw_cyc_distr_w1w3.png", path = figure_path, width = 40, height = 25, units = "cm")


# #---------- distributions of wealth, gini --------------
# wealth_distr_w1 <- ggplot(vill_df, aes(x=village_wealth_index_w1)) +
#   geom_histogram(color="black", fill="gray") +
#   labs(title = "Village wealth distribution (wave 1), n=176",
#        x = "Village wealth index (avg household wealth quintile, wave 1)") +
#   theme_classic()
# 
# wealth_distr_w3 <- ggplot(vill_df, aes(x=village_wealth_index_w3)) +
#   geom_histogram(color="black", fill="gray") +
#   labs(title = "Village wealth distribution (wave 3), n=174",
#        x = "Village wealth index (avg household wealth quintile, wave 3)") +
#   theme_classic()
# 
# gini_distr_w1 <- ggplot(gini_df, aes(x=gini_w1)) +
#   geom_histogram(color="black", fill="gray") +
#   labs(title = "Village gini distribution (wave 1), n=176",
#        x = "Village gini from household wealth (wave 1)") +
#   theme_classic()
# 
# gini_distr_w3 <- ggplot(gini_df, aes(x=gini_w3)) +
#   geom_histogram(color="black", fill="gray") +
#   labs(title = "Village gini distribution (wave 3), n=174",
#        x = "Village gini from houeshold wealth (wave 3)") +
#   theme_classic()
# 
# ggarrange(wealth_distr_w1, wealth_distr_w3, gini_distr_w1, gini_distr_w3,
#           ncol = 2, nrow = 2)
# ggsave("wealth_gini_distr.png", path = figure_path, width = 40, height = 25, units = "cm")
# 
# 
# #---- visualize different cyclic density metrics ------
# cycdens_er_df <- bind_rows("Indiv. borrow ties" = cycdens, 
#                             "Indiv. borrow and lend ties" = borrow_lend_cycdens,
#                             "Household borrow ties" = hh_cycdens, 
#                             "Household borrow and lend ties" = hh_borrow_lend_cycdens, 
#                             .id = "network_type") %>%
#   mutate(rand_graph_method = "Erdos-Renyi")
# cycdens_sds_df <- bind_rows("Indiv. borrow ties" = cycdens_sds,
#                             "Indiv. borrow and lend ties" = borrow_lend_cycdens_sds,
#                             "Household borrow ties" = hh_cycdens_sds,
#                             "Household borrow and lend ties" = hh_borrow_lend_cycdens_sds,
#                             .id = "network_type") %>%
#   mutate(rand_graph_method = "Same degree distribution")
# cycdens_pl_df <- bind_rows("Indiv. borrow ties" = cycdens_pl, 
#                             "Indiv. borrow and lend ties" = borrow_lend_cycdens_pl,
#                             "Household borrow ties" = hh_cycdens_pl, 
#                             "Household borrow and lend ties" = hh_borrow_lend_cycdens_pl, 
#                             .id = "network_type") %>%
#   mutate(rand_graph_method = "Power law degree distribution")
# cycdens_ba_df <- bind_rows("Indiv. borrow ties" = cycdens_ba, 
#                             "Indiv. borrow and lend ties" = borrow_lend_cycdens_ba,
#                             "Household borrow ties" = hh_cycdens_ba, 
#                             "Household borrow and lend ties" = hh_borrow_lend_cycdens_ba, 
#                             .id = "network_type") %>%
#   mutate(rand_graph_method = "Barabasi-Albert")
# 
# cycdens_df <- cycdens_er_df %>%
#   bind_rows(cycdens_sds_df) %>%
#   bind_rows(cycdens_pl_df) %>%
#   bind_rows(cycdens_ba_df)
# 
# var_label(cycdens_df) <- list(
#   network_type = "Type of network",
#   rand_graph_method = "Random graph method"
# )
# 
# # cycles of length 2
# cycdens2_sds <- ggplot(data=cycdens_sds_df, aes(x=kcycle2, color=network_type, fill=network_type)) +
#   #geom_histogram(aes(y=..density..), position="identity", alpha=0.4) +
#   geom_density(alpha=0.3) +
#   theme_classic() +
#   xlim(-3,7) +
#   ylim(0,1) +
#   labs(title="Cycles of length 2",
#        x = "Cyclic densities (same degree distribution)")
# 
# # cycles of length 3
# cycdens3_sds <- ggplot(data=cycdens_sds_df, aes(x=kcycle3, color=network_type, fill=network_type)) +
#   #geom_histogram(aes(y=..density..), position="identity", alpha=0.4) +
#   geom_density(alpha=0.3) +
#   theme_classic() +
#   xlim(-3,6) +
#   ylim(0,1.2) +
#   labs(title="Cycles of length 3",
#        x = "Cyclic densities (same degree distribution)")
# 
# # cycles of length 4
# cycdens4_sds <- ggplot(data=cycdens_sds_df, aes(x=kcycle4, color=network_type, fill=network_type)) +
#   #geom_histogram(aes(y=..density..), position="identity", alpha=0.4) +
#   geom_density(alpha=0.3) +
#   theme_classic() +
#   xlim(-3,8) +
#   ylim(0,3) +
#   labs(title="Cycles of length 4",
#        x = "Cyclic densities (same degree distribution)")
# 
# # cycles of length 5
# cycdens5_sds <- ggplot(data=cycdens_sds_df, aes(x=kcycle5, color=network_type, fill=network_type)) +
#   #geom_histogram(aes(y=..density..), position="identity", alpha=0.4) +
#   geom_density(alpha=0.3) +
#   theme_classic() +
#   xlim(-4,9) +
#   ylim(0,3.5) +
#   labs(title="Cycles of length 5",
#        x = "Cyclic densities (same degree distribution)")
# 
# ggarrange(cycdens2, cycdens3, cycdens4, cycdens5,
#           ncol = 2, nrow = 2)
# 
# 
# # #------- observed raw cycle counts by network type --------
# # kstats_combined_df <- bind_rows("Indiv. borrow ties" = kstats_df %>%
# #                             select(-c(village_code, relationship)), 
# #                             "Indiv. borrow and lend ties" = kstats_borrow_lend_df %>%
# #                             select(-c(village_code, relationship)),
# #                             "Household borrow ties" = kstats_hh_df %>%
# #                             select(-c(village_code, relationship)), 
# #                             "Household borrow and lend ties" = kstats_hh_borrow_lend_df %>%
# #                             select(-c(village_code, relationship)), 
# #                             .id = "network_type")
# # 
# # var_label(kstats_combined_df) <- list(
# #   network_type = "Type of network"
# # )
# # 
# # cyc_count_df <- kstats_combined_df %>%
# #   select(network_type, kcycle2, kcycle3, kcycle4, kcycle5)
# # 
# # cyccount2 <- ggplot(data=cyc_count_df, aes(x=kcycle2, color=network_type, fill=network_type)) +
# #   #geom_histogram(aes(y=..density..), position="identity", alpha=0.4) +
# #   geom_density(alpha=0.3) +
# #   theme_classic() +
# #   labs(title="Cycles of length 2",
# #        x = "Cyclic densities (erdos-renyi)")
# # 
# # cyccount3 <- ggplot(data=cyc_count_df, aes(x=kcycle3, color=network_type, fill=network_type)) +
# #   #geom_histogram(aes(y=..density..), position="identity", alpha=0.4) +
# #   geom_density(alpha=0.3) +
# #   theme_classic() +
# #   labs(title="Cycles of length 3",
# #        x = "Cyclic densities (erdos-renyi)")
# # 
# # cyccount4 <- ggplot(data=cyc_count_df, aes(x=kcycle4, color=network_type, fill=network_type)) +
# #   #geom_histogram(aes(y=..density..), position="identity", alpha=0.4) +
# #   geom_density(alpha=0.3) +
# #   theme_classic() +
# #   labs(title="Cycles of length 4",
# #        x = "Cyclic densities (erdos-renyi)")
# # 
# # cyccount5 <- ggplot(data=cyc_count_df, aes(x=kcycle5, color=network_type, fill=network_type)) +
# #   #geom_histogram(aes(y=..density..), position="identity", alpha=0.4) +
# #   geom_density(alpha=0.3) +
# #   theme_classic() +
# #   labs(title="Cycles of length 5",
# #        x = "Cyclic densities (erdos-renyi)")
# # 
# # ggarrange(cyccount2, cyccount3, cyccount4, cyccount5,
# #           ncol = 2, nrow = 2)
