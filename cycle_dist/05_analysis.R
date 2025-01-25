# code by Selena Lee; work for kcycles in undirected networks (borrow-lend) in Honduras (11 Sept 2023)
# How might the cyclical structure of social relations, at a village and individual level, affect
# the accumulation of community wealth over time in rural Honduran villages?

# Shared folder with Shiv: setwd("/SCRATCH/tmp/kcyle networks")

library(tidyverse)
library(dplyr)
library(lme4)
library(nnet)
library(dfidx)
library(mlogit)
library(mclogit)
library(scales)

library(jtools)
library(skimr)

library(ggpubr)
library(gtsummary)

library(ggplot2)
library(pheatmap)

library(performance)
library(see)

#-------------------------- Load data --------------------------------
#setwd("/WORKAREA/work/HONDURAS_GATES/BORROWING")
setwd("C:/Users/Shiv/Documents/Honduras/Kcycle and networks 04-19-23/kcycle_new_s/00_selena_final_code")

# Create a data path for saving cleaned / augmented data sets
#data_path <- "/WORKAREA/work/HONDURAS_GATES/BORROWING/00_selena_final_code/00_data"
figure_path <- paste0(data_path, "/figures")

# Load data for analysis
load("00_data/analysis_20240730.rda")

load("00_data/original_netwrks_resp_data_20230911.rda")
#------------------------------------------------------


# save(kstats_indiv_borrow_lend_no_kin_df,
#      file = file.path(data_path, "kstats_undir_borrow_lend_no_kin_20240516.rda"))
#----------------------------- To-do -----------------------------
# Re-run kcycle density calculations and run with model (also include density as control var in lm's)
# Run model looking at food security questions (vill-level avg) and presence of kcycle4
# Viz to do: distr of vill size, distr of vill wealth, distr of proportion of nodes in borrowing network
#   vs total nodes in vill, distr of proportion of ties in borrowing network vs total vill
# Log kcycle measures for kcycle density calculations (can also run using quadratic if relationship
#   is more pronounced)
# Create random network graphs to control for kcycles above what would be expected to be generated
#   by a random network

##### Models with both path and cyclic density included
# village wealth - quintile avg
# lm_cycpath2_vill_wealth_sds <- lmer(village_wealth_index ~ cycdens2_sds + pathdens2_sds + time +
#                                     (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_vill_wealth_sds <- lmer(village_wealth_index ~ cycdens3_sds + pathdens3_sds + time + 
#                                     (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_vill_wealth_sds <- lmer(village_wealth_index ~ cycdens4_sds + pathdens4_sds + time + 
#                                     (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_vill_wealth_sds <- lmer(village_wealth_index ~ cycdens5_sds + pathdens5_sds + time +
#                                     (1|village_code), data = kstats_borrow_lend_dir_kin_df)


#-------------------- household wealth -----------------------


#--------------------- village controls ----------------------
borrow_lend_df <- kstats_ud_bl_dh_df %>%
  mutate(time_to_main_road = ifelse(wave==1, time_to_main_road_w1, time_to_main_road_w3))

friend_df <- kstats_ud_fr_dh_df %>%
  mutate(time_to_main_road = ifelse(wave==1, time_to_main_road_w1, time_to_main_road_w3))

# village wealth - MCA
lm_cycpath3_ud_bl_dh_df_mca <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
                                     vill_size + time_to_main_road +
                                     (1|village_code), data = borrow_lend_df)
lm_cycpath4_ud_bl_dh_df_mca <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
                                     vill_size + time_to_main_road +
                                     (1|village_code), data = borrow_lend_df)
lm_cycpath5_ud_bl_dh_df_mca <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
                                     vill_size + time_to_main_road +
                                     (1|village_code), data = borrow_lend_df)

lm_cycpath3_ud_fr_dh_df_mca <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
                                     vill_size + time_to_main_road +
                                     (1|village_code), data = friend_df)
lm_cycpath4_ud_fr_dh_df_mca <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
                                     vill_size + time_to_main_road +
                                     (1|village_code), data = friend_df)
lm_cycpath5_ud_fr_dh_df_mca <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
                                     vill_size + time_to_main_road +
                                     (1|village_code), data = friend_df)


summ(lm_cycpath3_ud_bl_dh_df_mca, digits=3)
summ(lm_cycpath4_ud_bl_dh_df_mca, digits=3)
summ(lm_cycpath5_ud_bl_dh_df_mca, digits=3)

summ(lm_cycpath3_ud_fr_dh_df_mca, digits=3)
summ(lm_cycpath4_ud_fr_dh_df_mca, digits=3)
summ(lm_cycpath5_ud_fr_dh_df_mca, digits=3)

lm_cycpath3_bl_mca_tbl <- tbl_regression(lm_cycpath3_ud_bl_dh_df_mca)
lm_cycpath4_bl_mca_tbl <- tbl_regression(lm_cycpath4_ud_bl_dh_df_mca)
lm_cycpath5_bl_mca_tbl <- tbl_regression(lm_cycpath5_ud_bl_dh_df_mca)

lm_cycpath3_fr_mca_tbl <- tbl_regression(lm_cycpath3_ud_fr_dh_df_mca)
lm_cycpath4_fr_mca_tbl <- tbl_regression(lm_cycpath4_ud_fr_dh_df_mca)
lm_cycpath5_fr_mca_tbl <- tbl_regression(lm_cycpath5_ud_fr_dh_df_mca)

library(cardx)
tbl_cycpath3_ud_bl_dh_df_mca <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
                                     vill_size + time_to_main_road +
                                     (1|village_code), data = borrow_lend_df) %>%
  tbl_regression(tidy_fun = broom.mixed::tidy,
                 conf.int=FALSE) %>%
  add_global_p()
tbl_cycpath4_ud_bl_dh_df_mca <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
                                     vill_size + time_to_main_road +
                                     (1|village_code), data = borrow_lend_df) %>%
  tbl_regression(tidy_fun = broom.mixed::tidy,
                 conf.int=FALSE) %>%
  add_global_p()
tbl_cycpath5_ud_bl_dh_df_mca <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
                                     vill_size + time_to_main_road +
                                     (1|village_code), data = borrow_lend_df) %>%
  tbl_regression(tidy_fun = broom.mixed::tidy,
                 conf.int=FALSE) %>%
  add_global_p()

tbl_cycpath3_ud_fr_dh_df_mca <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
                                     vill_size + time_to_main_road +
                                     (1|village_code), data = friend_df) %>%
  tbl_regression(tidy_fun = broom.mixed::tidy,
                 conf.int=FALSE) %>%
  add_global_p()
tbl_cycpath4_ud_fr_dh_df_mca <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
                                     vill_size + time_to_main_road +
                                     (1|village_code), data = friend_df) %>%
  tbl_regression(tidy_fun = broom.mixed::tidy,
                 conf.int=FALSE) %>%
  add_global_p()
tbl_cycpath5_ud_fr_dh_df_mca <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
                                     vill_size + time_to_main_road +
                                     (1|village_code), data = friend_df) %>%
  tbl_regression(tidy_fun = broom.mixed::tidy,
                 conf.int=FALSE) %>%
  add_global_p()

# Table-formatted output (borrow-lend network)
tbl_borrow_lend <-
  tbl_merge(
    tbls = list(tbl_cycpath3_ud_bl_dh_df_mca, tbl_cycpath4_ud_bl_dh_df_mca, tbl_cycpath5_ud_bl_dh_df_mca),
    tab_spanner = c("**Cycles of length 3**", "**Cycles of length 4**", "**Cycles of length 5**")
  )

# Table-formatted output (friend network)
tbl_friend <-
  tbl_merge(
    tbls = list(tbl_cycpath3_ud_fr_dh_df_mca, tbl_cycpath4_ud_fr_dh_df_mca, tbl_cycpath5_ud_fr_dh_df_mca),
    tab_spanner = c("**Cycles of length 3**", "**Cycles of length 4**", "**Cycles of length 5**")
  )

#---------------------- outcome: village wealth (MCA) -------------------------
#### Cannot be run after this 08-17-24 -- SHIV

# village wealth - MCA
lm_cycpath2_d_bl_k_df_mca <- lmer(vill_wealth_mca ~ cycdens2_sds + pathdens2_sds + time +
                                      (1|village_code), data = kstats_ud_bl_k_df)
lm_cycpath3_d_bl_k_df_mca <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
                                      (1|village_code), data = kstats_d_bl_k_df)
lm_cycpath4_d_bl_k_df_mca <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
                                      (1|village_code), data = kstats_d_bl_k_df)
lm_cycpath5_d_bl_k_df_mca <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
                                      (1|village_code), data = kstats_d_bl_k_df)
# village wealth - MCA
lm_cycpath2_d_bl_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens2_sds + pathdens2_sds + time +
                                     (1|village_code), data = kstats_d_bl_nk_df)
lm_cycpath3_d_bl_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
                                     (1|village_code), data = kstats_d_bl_nk_df)
lm_cycpath4_d_bl_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
                                     (1|village_code), data = kstats_d_bl_nk_df)
lm_cycpath5_d_bl_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
                                     (1|village_code), data = kstats_d_bl_nk_df)
# village wealth - MCA
lm_cycpath2_ud_bl_k_df_mca <- lmer(vill_wealth_mca ~ cycdens2_sds + pathdens2_sds + time +
                                     (1|village_code), data = kstats_ud_bl_k_df)
lm_cycpath3_ud_bl_k_df_mca <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
                                     (1|village_code), data = kstats_ud_bl_k_df)
lm_cycpath4_ud_bl_k_df_mca <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
                                     (1|village_code), data = kstats_ud_bl_k_df)
lm_cycpath5_ud_bl_k_df_mca <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
                                     (1|village_code), data = kstats_ud_bl_k_df)
# village wealth - MCA
lm_cycpath2_ud_bl_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens2_sds + pathdens2_sds + time +
                                      (1|village_code), data = kstats_ud_bl_nk_df)
lm_cycpath3_ud_bl_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
                                      (1|village_code), data = kstats_ud_bl_nk_df)
lm_cycpath4_ud_bl_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
                                      (1|village_code), data = kstats_ud_bl_nk_df)
lm_cycpath5_ud_bl_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
                                      (1|village_code), data = kstats_ud_bl_nk_df)
# village wealth - MCA
lm_cycpath2_d_fr_k_df_mca <- lmer(vill_wealth_mca ~ cycdens2_sds + pathdens2_sds + time +
                                    (1|village_code), data = kstats_d_fr_k_df)
lm_cycpath3_d_fr_k_df_mca <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
                                    (1|village_code), data = kstats_d_fr_k_df)
lm_cycpath4_d_fr_k_df_mca <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
                                    (1|village_code), data = kstats_d_fr_k_df)
lm_cycpath5_d_fr_k_df_mca <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
                                    (1|village_code), data = kstats_d_fr_k_df)
# village wealth - MCA
lm_cycpath2_d_fr_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens2_sds + pathdens2_sds + time +
                                     (1|village_code), data = kstats_d_fr_nk_df)
lm_cycpath3_d_fr_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
                                     (1|village_code), data = kstats_d_fr_nk_df)
lm_cycpath4_d_fr_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
                                     (1|village_code), data = kstats_d_fr_nk_df)
lm_cycpath5_d_fr_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
                                     (1|village_code), data = kstats_d_fr_nk_df)
# village wealth - MCA
lm_cycpath2_ud_fr_k_df_mca <- lmer(vill_wealth_mca ~ cycdens2_sds + pathdens2_sds + time +
                                      (1|village_code), data = kstats_ud_fr_k_df)
lm_cycpath3_ud_fr_k_df_mca <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
                                      (1|village_code), data = kstats_ud_fr_k_df)
lm_cycpath4_ud_fr_k_df_mca <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
                                      (1|village_code), data = kstats_ud_fr_k_df)
lm_cycpath5_ud_fr_k_df_mca <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
                                      (1|village_code), data = kstats_ud_fr_k_df)
# village wealth - MCA
lm_cycpath2_ud_fr_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens2_sds + pathdens2_sds + time +
                                      (1|village_code), data = kstats_ud_fr_nk_df)
lm_cycpath3_ud_fr_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
                                      (1|village_code), data = kstats_ud_fr_nk_df)
lm_cycpath4_ud_fr_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
                                      (1|village_code), data = kstats_ud_fr_nk_df)
lm_cycpath5_ud_fr_nk_df_mca <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
                                      (1|village_code), data = kstats_ud_fr_nk_df)

# #--------------- outcome: village wealth (MCA) with random slope (time) --------------------
# # I tried to run the model with  random slope and random intercept, but ran into an error.
# # Because I would need to include both a random slope and random intercept in the same model,
# # there would be one random effect per observation. Because of variance in the initial timepoint
# # in villages, I shouldn't use only a random slope without a random intercept.
#
# # Per post on stack overflow: https://stackoverflow.com/questions/26465215/random-slope-for-time-in-subject-not-working-in-lme4
# # The sensible thing to do is to recognize that the variation among slopes is unidentifiable; 
# # there may be among-individual variation among slopes, but you just can't estimate it with 
# # this model. Don't try; fit a random-intercept model and let the implicit/default random error 
# # term take care of the variation among slopes.
# 
# kstats_d_bl_k_df_subset <- subset(kstats_d_bl_k_df, !(village_code %in% c(155, 156)))
# kstats_d_bl_nk_df_subset <- subset(kstats_d_bl_nk_df, !(village_code %in% c(155, 156)))
# kstats_ud_bl_k_df_subset <- subset(kstats_ud_bl_k_df, !(village_code %in% c(155, 156)))
# kstats_ud_bl_nk_df_subset <- subset(kstats_ud_bl_nk_df, !(village_code %in% c(155, 156)))
# kstats_d_fr_k_df_subset <- subset(kstats_d_fr_k_df, !(village_code %in% c(155, 156)))
# kstats_d_fr_nk_df_subset <- subset(kstats_d_fr_nk_df, !(village_code %in% c(155, 156)))
# # kstats_ud_fr_k_df_subset <- subset(kstats_ud_fr_k_df, !(village_code %in% c(155, 156)))
# kstats_ud_fr_nk_df_subset <- subset(kstats_ud_fr_nk_df, !(village_code %in% c(155, 156)))
# 
# # village wealth - MCA
# lm_cycpath2_d_bl_k_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens2_sds + pathdens2_sds + time +
#                                     (time|village_code), data = kstats_d_bl_k_df_subset)
# lm_cycpath3_d_bl_k_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
#                                     (1+time|village_code), data = kstats_d_bl_k_df_subset)
# lm_cycpath4_d_bl_k_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
#                                     (1+time|village_code), data = kstats_d_bl_k_df_subset)
# lm_cycpath5_d_bl_k_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
#                                     (1+time|village_code), data = kstats_d_bl_k_df_subset)
# # village wealth - MCA
# lm_cycpath2_d_bl_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens2_sds + pathdens2_sds + time +
#                                      (1+time|village_code), data = kstats_d_bl_nk_df_subset)
# lm_cycpath3_d_bl_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
#                                      (1+time|village_code), data = kstats_d_bl_nk_df_subset)
# lm_cycpath4_d_bl_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
#                                      (1+time|village_code), data = kstats_d_bl_nk_df_subset)
# lm_cycpath5_d_bl_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
#                                      (1+time|village_code), data = kstats_d_bl_nk_df_subset)
# # village wealth - MCA
# lm_cycpath2_ud_bl_k_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens2_sds + pathdens2_sds + time +
#                                      (1+time|village_code), data = kstats_ud_bl_k_df_subset)
# lm_cycpath3_ud_bl_k_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
#                                      (1+time|village_code), data = kstats_ud_bl_k_df_subset)
# lm_cycpath4_ud_bl_k_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
#                                      (1+time|village_code), data = kstats_ud_bl_k_df_subset)
# lm_cycpath5_ud_bl_k_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
#                                      (1+time|village_code), data = kstats_ud_bl_k_df_subset)
# # village wealth - MCA
# lm_cycpath2_ud_bl_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens2_sds + pathdens2_sds + time +
#                                       (1+time|village_code), data = kstats_ud_bl_nk_df_subset)
# lm_cycpath3_ud_bl_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
#                                       (1+time|village_code), data = kstats_ud_bl_nk_df_subset)
# lm_cycpath4_ud_bl_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
#                                       (1+time|village_code), data = kstats_ud_bl_nk_df_subset)
# lm_cycpath5_ud_bl_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
#                                       (1+time|village_code), data = kstats_ud_bl_nk_df_subset)
# # village wealth - MCA
# lm_cycpath2_d_fr_k_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens2_sds + pathdens2_sds + time +
#                                     (1+time|village_code), data = kstats_d_fr_k_df_subset)
# lm_cycpath3_d_fr_k_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
#                                     (1+time|village_code), data = kstats_d_fr_k_df_subset)
# lm_cycpath4_d_fr_k_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
#                                     (1+time|village_code), data = kstats_d_fr_k_df_subset)
# lm_cycpath5_d_fr_k_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
#                                     (1+time|village_code), data = kstats_d_fr_k_df_subset)
# # village wealth - MCA
# lm_cycpath2_d_fr_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens2_sds + pathdens2_sds + time +
#                                      (1+time|village_code), data = kstats_d_fr_nk_df_subset)
# lm_cycpath3_d_fr_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
#                                      (1+time|village_code), data = kstats_d_fr_nk_df_subset)
# lm_cycpath4_d_fr_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
#                                      (1+time|village_code), data = kstats_d_fr_nk_df_subset)
# lm_cycpath5_d_fr_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
#                                      (1+time|village_code), data = kstats_d_fr_nk_df_subset)
# # village wealth - MCA
# lm_cycpath2_ud_fr_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens2_sds + pathdens2_sds + time +
#                                       (1+time|village_code), data = kstats_ud_fr_nk_df_subset)
# lm_cycpath3_ud_fr_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens3_sds + pathdens3_sds + time + 
#                                       (1+time|village_code), data = kstats_ud_fr_nk_df_subset)
# lm_cycpath4_ud_fr_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens4_sds + pathdens4_sds + time + 
#                                       (1+time|village_code), data = kstats_ud_fr_nk_df_subset)
# lm_cycpath5_ud_fr_nk_df_mca_randslope <- lmer(vill_wealth_mca ~ cycdens5_sds + pathdens5_sds + time +
#                                       (1+time|village_code), data = kstats_ud_fr_nk_df_subset)

#-------------------- outcome: village wealth (avg quintiles) -----------------------
# village wealth - avg quintiles
lm_cycpath2_d_bl_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens2_sds + pathdens2_sds + time +
                                            (1|village_code), data = kstats_d_bl_k_df)
lm_cycpath3_d_bl_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens3_sds + pathdens3_sds + time + 
                                            (1|village_code), data = kstats_d_bl_k_df)
lm_cycpath4_d_bl_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens4_sds + pathdens4_sds + time + 
                                            (1|village_code), data = kstats_d_bl_k_df)
lm_cycpath5_d_bl_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens5_sds + pathdens5_sds + time +
                                            (1|village_code), data = kstats_d_bl_k_df)
# village wealth - avg quintiles
lm_cycpath2_d_bl_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens2_sds + pathdens2_sds + time +
                                            (1|village_code), data = kstats_d_bl_nk_df)
lm_cycpath3_d_bl_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens3_sds + pathdens3_sds + time + 
                                            (1|village_code), data = kstats_d_bl_nk_df)
lm_cycpath4_d_bl_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens4_sds + pathdens4_sds + time + 
                                            (1|village_code), data = kstats_d_bl_nk_df)
lm_cycpath5_d_bl_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens5_sds + pathdens5_sds + time +
                                            (1|village_code), data = kstats_d_bl_nk_df)
# village wealth - avg quintiles
lm_cycpath2_ud_bl_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens2_sds + pathdens2_sds + time +
                                             (1|village_code), data = kstats_ud_bl_k_df)
lm_cycpath3_ud_bl_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens3_sds + pathdens3_sds + time + 
                                             (1|village_code), data = kstats_ud_bl_k_df)
lm_cycpath4_ud_bl_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens4_sds + pathdens4_sds + time + 
                                             (1|village_code), data = kstats_ud_bl_k_df)
lm_cycpath5_ud_bl_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens5_sds + pathdens5_sds + time +
                                             (1|village_code), data = kstats_ud_bl_k_df)
# village wealth - avg quintiles
lm_cycpath2_ud_bl_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens2_sds + pathdens2_sds + time +
                                              (1|village_code), data = kstats_ud_bl_nk_df)
lm_cycpath3_ud_bl_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens3_sds + pathdens3_sds + time + 
                                              (1|village_code), data = kstats_ud_bl_nk_df)
lm_cycpath4_ud_bl_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens4_sds + pathdens4_sds + time + 
                                              (1|village_code), data = kstats_ud_bl_nk_df)
lm_cycpath5_ud_bl_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens5_sds + pathdens5_sds + time +
                                              (1|village_code), data = kstats_ud_bl_nk_df)
# village wealth - avg quintiles
lm_cycpath2_d_fr_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens2_sds + pathdens2_sds + time +
                                            (1|village_code), data = kstats_d_fr_k_df)
lm_cycpath3_d_fr_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens3_sds + pathdens3_sds + time + 
                                            (1|village_code), data = kstats_d_fr_k_df)
lm_cycpath4_d_fr_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens4_sds + pathdens4_sds + time + 
                                            (1|village_code), data = kstats_d_fr_k_df)
lm_cycpath5_d_fr_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens5_sds + pathdens5_sds + time +
                                            (1|village_code), data = kstats_d_fr_k_df)
# village wealth - avg quintiles
lm_cycpath2_d_fr_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens2_sds + pathdens2_sds + time +
                                             (1|village_code), data = kstats_d_fr_nk_df)
lm_cycpath3_d_fr_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens3_sds + pathdens3_sds + time + 
                                             (1|village_code), data = kstats_d_fr_nk_df)
lm_cycpath4_d_fr_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens4_sds + pathdens4_sds + time + 
                                             (1|village_code), data = kstats_d_fr_nk_df)
lm_cycpath5_d_fr_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens5_sds + pathdens5_sds + time +
                                             (1|village_code), data = kstats_d_fr_nk_df)
# village wealth - avg quintiles
lm_cycpath2_ud_fr_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens2_sds + pathdens2_sds + time +
                                              (1|village_code), data = kstats_ud_fr_k_df)
lm_cycpath3_ud_fr_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens3_sds + pathdens3_sds + time + 
                                              (1|village_code), data = kstats_ud_fr_k_df)
lm_cycpath4_ud_fr_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens4_sds + pathdens4_sds + time + 
                                              (1|village_code), data = kstats_ud_fr_k_df)
lm_cycpath5_ud_fr_k_df_vill_wealth <- lmer(village_wealth_index ~ cycdens5_sds + pathdens5_sds + time +
                                              (1|village_code), data = kstats_ud_fr_k_df)
# village wealth - avg quintiles
lm_cycpath2_ud_fr_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens2_sds + pathdens2_sds + time +
                                              (1|village_code), data = kstats_ud_fr_nk_df)
lm_cycpath3_ud_fr_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens3_sds + pathdens3_sds + time + 
                                              (1|village_code), data = kstats_ud_fr_nk_df)
lm_cycpath4_ud_fr_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens4_sds + pathdens4_sds + time + 
                                              (1|village_code), data = kstats_ud_fr_nk_df)
lm_cycpath5_ud_fr_nk_df_vill_wealth <- lmer(village_wealth_index ~ cycdens5_sds + pathdens5_sds + time +
                                              (1|village_code), data = kstats_ud_fr_nk_df)

#----------------------- outcome: village gini (MCA) ---------------------------
# gini - from MCA coord
lm_cycpath2_d_bl_k_df_gini <- lmer(gini_frm_cont ~ cycdens2_sds + pathdens2_sds + time +
                                     (1|village_code), data = kstats_d_bl_k_df)
lm_cycpath3_d_bl_k_df_gini <- lmer(gini_frm_cont ~ cycdens3_sds + pathdens3_sds + time + 
                                     (1|village_code), data = kstats_d_bl_k_df)
lm_cycpath4_d_bl_k_df_gini <- lmer(gini_frm_cont ~ cycdens4_sds + pathdens4_sds + time + 
                                     (1|village_code), data = kstats_d_bl_k_df)
lm_cycpath5_d_bl_k_df_gini <- lmer(gini_frm_cont ~ cycdens5_sds + pathdens5_sds + time +
                                     (1|village_code), data = kstats_d_bl_k_df)
# gini - from MCA coord
lm_cycpath2_d_bl_nk_df_gini <- lmer(gini_frm_cont ~ cycdens2_sds + pathdens2_sds + time +
                                     (1|village_code), data = kstats_d_bl_nk_df)
lm_cycpath3_d_bl_nk_df_gini <- lmer(gini_frm_cont ~ cycdens3_sds + pathdens3_sds + time + 
                                     (1|village_code), data = kstats_d_bl_nk_df)
lm_cycpath4_d_bl_nk_df_gini <- lmer(gini_frm_cont ~ cycdens4_sds + pathdens4_sds + time + 
                                     (1|village_code), data = kstats_d_bl_nk_df)
lm_cycpath5_d_bl_nk_df_gini <- lmer(gini_frm_cont ~ cycdens5_sds + pathdens5_sds + time +
                                     (1|village_code), data = kstats_d_bl_nk_df)
# gini - from MCA coord
lm_cycpath2_ud_bl_k_df_gini <- lmer(gini_frm_cont ~ cycdens2_sds + pathdens2_sds + time +
                                     (1|village_code), data = kstats_ud_bl_k_df)
lm_cycpath3_ud_bl_k_df_gini <- lmer(gini_frm_cont ~ cycdens3_sds + pathdens3_sds + time + 
                                     (1|village_code), data = kstats_ud_bl_k_df)
lm_cycpath4_ud_bl_k_df_gini <- lmer(gini_frm_cont ~ cycdens4_sds + pathdens4_sds + time + 
                                     (1|village_code), data = kstats_ud_bl_k_df)
lm_cycpath5_ud_bl_k_df_gini <- lmer(gini_frm_cont ~ cycdens5_sds + pathdens5_sds + time +
                                     (1|village_code), data = kstats_ud_bl_k_df)
# gini - from MCA coord
lm_cycpath2_ud_bl_nk_df_gini <- lmer(gini_frm_cont ~ cycdens2_sds + pathdens2_sds + time +
                                      (1|village_code), data = kstats_ud_bl_nk_df)
lm_cycpath3_ud_bl_nk_df_gini <- lmer(gini_frm_cont ~ cycdens3_sds + pathdens3_sds + time + 
                                      (1|village_code), data = kstats_ud_bl_nk_df)
lm_cycpath4_ud_bl_nk_df_gini <- lmer(gini_frm_cont ~ cycdens4_sds + pathdens4_sds + time + 
                                      (1|village_code), data = kstats_ud_bl_nk_df)
lm_cycpath5_ud_bl_nk_df_gini <- lmer(gini_frm_cont ~ cycdens5_sds + pathdens5_sds + time +
                                      (1|village_code), data = kstats_ud_bl_nk_df)
# gini - from MCA coord
lm_cycpath2_d_fr_k_df_gini <- lmer(gini_frm_cont ~ cycdens2_sds + pathdens2_sds + time +
                                     (1|village_code), data = kstats_d_fr_k_df)
lm_cycpath3_d_fr_k_df_gini <- lmer(gini_frm_cont ~ cycdens3_sds + pathdens3_sds + time + 
                                     (1|village_code), data = kstats_d_fr_k_df)
lm_cycpath4_d_fr_k_df_gini <- lmer(gini_frm_cont ~ cycdens4_sds + pathdens4_sds + time + 
                                     (1|village_code), data = kstats_d_fr_k_df)
lm_cycpath5_d_fr_k_df_gini <- lmer(gini_frm_cont ~ cycdens5_sds + pathdens5_sds + time +
                                     (1|village_code), data = kstats_d_fr_k_df)
# gini - from MCA coord
lm_cycpath2_d_fr_nk_df_gini <- lmer(gini_frm_cont ~ cycdens2_sds + pathdens2_sds + time +
                                     (1|village_code), data = kstats_d_fr_nk_df)
lm_cycpath3_d_fr_nk_df_gini <- lmer(gini_frm_cont ~ cycdens3_sds + pathdens3_sds + time + 
                                     (1|village_code), data = kstats_d_fr_nk_df)
lm_cycpath4_d_fr_nk_df_gini <- lmer(gini_frm_cont ~ cycdens4_sds + pathdens4_sds + time + 
                                     (1|village_code), data = kstats_d_fr_nk_df)
lm_cycpath5_d_fr_nk_df_gini <- lmer(gini_frm_cont ~ cycdens5_sds + pathdens5_sds + time +
                                     (1|village_code), data = kstats_d_fr_nk_df)
# gini - from MCA coord
lm_cycpath2_ud_fr_k_df_gini <- lmer(gini_frm_cont ~ cycdens2_sds + pathdens2_sds + time +
                                       (1|village_code), data = kstats_ud_fr_k_df)
lm_cycpath3_ud_fr_k_df_gini <- lmer(gini_frm_cont ~ cycdens3_sds + pathdens3_sds + time + 
                                       (1|village_code), data = kstats_ud_fr_k_df)
lm_cycpath4_ud_fr_k_df_gini <- lmer(gini_frm_cont ~ cycdens4_sds + pathdens4_sds + time + 
                                       (1|village_code), data = kstats_ud_fr_k_df)
lm_cycpath5_ud_fr_k_df_gini <- lmer(gini_frm_cont ~ cycdens5_sds + pathdens5_sds + time +
                                       (1|village_code), data = kstats_ud_fr_k_df)
# gini - from MCA coord
lm_cycpath2_ud_fr_nk_df_gini <- lmer(gini_frm_cont ~ cycdens2_sds + pathdens2_sds + time +
                                      (1|village_code), data = kstats_ud_fr_nk_df)
lm_cycpath3_ud_fr_nk_df_gini <- lmer(gini_frm_cont ~ cycdens3_sds + pathdens3_sds + time + 
                                      (1|village_code), data = kstats_ud_fr_nk_df)
lm_cycpath4_ud_fr_nk_df_gini <- lmer(gini_frm_cont ~ cycdens4_sds + pathdens4_sds + time + 
                                      (1|village_code), data = kstats_ud_fr_nk_df)
lm_cycpath5_ud_fr_nk_df_gini <- lmer(gini_frm_cont ~ cycdens5_sds + pathdens5_sds + time +
                                      (1|village_code), data = kstats_ud_fr_nk_df)

#------------------ model summaries ---------------------
##### MCA
# Sig increase in cycdens corresponds to inc vill wealth using wealth quintiles
# (cycdens2, 3, 4 significant)
summ(lm_cycpath2_d_bl_k_df_mca, digits=3)
summ(lm_cycpath2_d_bl_nk_df_mca, digits=3)
summ(lm_cycpath2_ud_bl_k_df_mca, digits=3)
summ(lm_cycpath2_ud_bl_nk_df_mca, digits=3)
summ(lm_cycpath2_d_fr_k_df_mca, digits=3)
summ(lm_cycpath2_d_fr_nk_df_mca, digits=3)
summ(lm_cycpath2_ud_fr_k_df_mca, digits=3)
summ(lm_cycpath2_ud_fr_nk_df_mca, digits=3)

summ(lm_cycpath3_d_bl_k_df_mca, digits=3)
summ(lm_cycpath3_d_bl_nk_df_mca, digits=3)
summ(lm_cycpath3_ud_bl_k_df_mca, digits=3)
summ(lm_cycpath3_ud_bl_nk_df_mca, digits=3)
summ(lm_cycpath3_d_fr_k_df_mca, digits=3)
summ(lm_cycpath3_d_fr_nk_df_mca, digits=3)
summ(lm_cycpath3_ud_fr_k_df_mca, digits=3)
summ(lm_cycpath3_ud_fr_nk_df_mca, digits=3)

summ(lm_cycpath4_d_bl_k_df_mca, digits=3)
summ(lm_cycpath4_d_bl_nk_df_mca, digits=3)
summ(lm_cycpath4_ud_bl_k_df_mca, digits=3)
summ(lm_cycpath4_ud_bl_nk_df_mca, digits=3)
summ(lm_cycpath4_d_fr_k_df_mca, digits=3)
summ(lm_cycpath4_d_fr_nk_df_mca, digits=3)
summ(lm_cycpath4_ud_fr_k_df_mca, digits=3)
summ(lm_cycpath4_ud_fr_nk_df_mca, digits=3)

summ(lm_cycpath5_d_bl_k_df_mca, digits=3)
summ(lm_cycpath5_d_bl_nk_df_mca, digits=3)
summ(lm_cycpath5_ud_bl_k_df_mca, digits=3)
summ(lm_cycpath5_ud_bl_nk_df_mca, digits=3)
summ(lm_cycpath5_d_fr_k_df_mca, digits=3)
summ(lm_cycpath5_d_fr_nk_df_mca, digits=3)
summ(lm_cycpath5_ud_fr_k_df_mca, digits=3)
summ(lm_cycpath5_ud_fr_nk_df_mca, digits=3)

##### Village wealth (avg quintiles)
# Sig increase in cycdens corresponds to inc vill wealth using wealth quintiles
# (cycdens2, 3, 4 significant)
summ(lm_cycpath2_d_bl_k_df_vill_wealth, digits=3)
summ(lm_cycpath2_d_bl_nk_df_vill_wealth, digits=3)
summ(lm_cycpath2_ud_bl_k_df_vill_wealth, digits=3)
summ(lm_cycpath2_ud_bl_nk_df_vill_wealth, digits=3)
summ(lm_cycpath2_d_fr_k_df_vill_wealth, digits=3)
summ(lm_cycpath2_d_fr_nk_df_vill_wealth, digits=3)
summ(lm_cycpath2_ud_fr_k_df_vill_wealth, digits=3)
summ(lm_cycpath2_ud_fr_nk_df_vill_wealth, digits=3)

summ(lm_cycpath3_d_bl_k_df_vill_wealth, digits=3)
summ(lm_cycpath3_d_bl_nk_df_vill_wealth, digits=3)
summ(lm_cycpath3_ud_bl_k_df_vill_wealth, digits=3)
summ(lm_cycpath3_ud_bl_nk_df_vill_wealth, digits=3)
summ(lm_cycpath3_d_fr_k_df_vill_wealth, digits=3)
summ(lm_cycpath3_d_fr_nk_df_vill_wealth, digits=3)
summ(lm_cycpath3_ud_fr_k_df_vill_wealth, digits=3)
summ(lm_cycpath3_ud_fr_nk_df_vill_wealth, digits=3)

summ(lm_cycpath4_d_bl_k_df_vill_wealth, digits=3)
summ(lm_cycpath4_d_bl_nk_df_vill_wealth, digits=3)
summ(lm_cycpath4_ud_bl_k_df_vill_wealth, digits=3)
summ(lm_cycpath4_ud_bl_nk_df_vill_wealth, digits=3)
summ(lm_cycpath4_d_fr_k_df_vill_wealth, digits=3)
summ(lm_cycpath4_d_fr_nk_df_vill_wealth, digits=3)
summ(lm_cycpath4_ud_fr_k_df_vill_wealth, digits=3)
summ(lm_cycpath4_ud_fr_nk_df_vill_wealth, digits=3)

summ(lm_cycpath5_d_bl_k_df_vill_wealth, digits=3)
summ(lm_cycpath5_d_bl_nk_df_vill_wealth, digits=3)
summ(lm_cycpath5_ud_bl_k_df_vill_wealth, digits=3)
summ(lm_cycpath5_ud_bl_nk_df_vill_wealth, digits=3)
summ(lm_cycpath5_d_fr_k_df_vill_wealth, digits=3)
summ(lm_cycpath5_d_fr_nk_df_vill_wealth, digits=3)
summ(lm_cycpath5_ud_fr_k_df_vill_wealth, digits=3)
summ(lm_cycpath5_ud_fr_nk_df_vill_wealth, digits=3)

##### vill gini
# Sig increase in cycdens corresponds to inc vill wealth using wealth quintiles
# (cycdens2, 3, 4 significant)
summ(lm_cycpath2_d_bl_k_df_gini, digits=3)
summ(lm_cycpath2_d_bl_nk_df_gini, digits=3)
summ(lm_cycpath2_ud_bl_k_df_gini, digits=3)
summ(lm_cycpath2_ud_bl_nk_df_gini, digits=3)
summ(lm_cycpath2_d_fr_k_df_gini, digits=3)
summ(lm_cycpath2_d_fr_nk_df_gini, digits=3)
summ(lm_cycpath2_ud_fr_k_df_gini, digits=3)
summ(lm_cycpath2_ud_fr_nk_df_gini, digits=3)

summ(lm_cycpath3_d_bl_k_df_gini, digits=3)
summ(lm_cycpath3_d_bl_nk_df_gini, digits=3)
summ(lm_cycpath3_ud_bl_k_df_gini, digits=3)
summ(lm_cycpath3_ud_bl_nk_df_gini, digits=3)
summ(lm_cycpath3_d_fr_k_df_gini, digits=3)
summ(lm_cycpath3_d_fr_nk_df_gini, digits=3)
summ(lm_cycpath3_ud_fr_k_df_gini, digits=3)
summ(lm_cycpath3_ud_fr_nk_df_gini, digits=3)

summ(lm_cycpath4_d_bl_k_df_gini, digits=3)
summ(lm_cycpath4_d_bl_nk_df_gini, digits=3)
summ(lm_cycpath4_ud_bl_k_df_gini, digits=3)
summ(lm_cycpath4_ud_bl_nk_df_gini, digits=3)
summ(lm_cycpath4_d_fr_k_df_gini, digits=3)
summ(lm_cycpath4_d_fr_nk_df_gini, digits=3)
summ(lm_cycpath4_ud_fr_k_df_gini, digits=3)
summ(lm_cycpath4_ud_fr_nk_df_gini, digits=3)

summ(lm_cycpath5_d_bl_k_df_gini, digits=3)
summ(lm_cycpath5_d_bl_nk_df_gini, digits=3)
summ(lm_cycpath5_ud_bl_k_df_gini, digits=3)
summ(lm_cycpath5_ud_bl_nk_df_gini, digits=3)
summ(lm_cycpath5_d_fr_k_df_gini, digits=3)
summ(lm_cycpath5_d_fr_nk_df_gini, digits=3)
summ(lm_cycpath5_ud_fr_k_df_gini, digits=3)
summ(lm_cycpath5_ud_fr_nk_df_gini, digits=3)

fig_vill_d_bl_k_cyc4_mca <- ggplot(kstats_d_bl_k_df, aes(cycdens4_sds, vill_wealth_mca)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_classic() +
  labs(title="lm_cycpath4_d_bl_k_df_mca")

fig_vill_d_bl_nk_cyc4_mca <- ggplot(kstats_d_bl_nk_df, aes(cycdens4_sds, vill_wealth_mca)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_classic() +
  labs(title="lm_cycpath4_d_bl_nk_df_mca")

fig_vill_ud_bl_k_cyc4_mca <- ggplot(kstats_ud_bl_k_df, aes(cycdens4_sds, vill_wealth_mca)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_classic() +
  labs(title="lm_cycpath4_ud_bl_k_df_mca")

fig_vill_ud_bl_nk_cyc4_mca <- ggplot(kstats_ud_bl_nk_df, aes(cycdens4_sds, vill_wealth_mca)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_classic() +
  labs(title="lm_cycpath4_ud_bl_nk_df_mca")

fig_vill_d_fr_k_cyc4_mca <- ggplot(kstats_d_fr_k_df, aes(cycdens4_sds, vill_wealth_mca)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_classic() +
  labs(title="lm_cycpath4_d_fr_k_df_mca")

fig_vill_d_fr_nk_cyc4_mca <- ggplot(kstats_d_fr_nk_df, aes(cycdens4_sds, vill_wealth_mca)) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_classic() +
  labs(title="lm_cycpath4_d_fr_nk_df_mca")

ggsave("vill_d_bl_k_cyc4_mca.png", plot=fig_vill_d_bl_k_cyc4_mca, path = figure_path,
       width = 20, height = 15, units = "cm")
ggsave("vill_d_bl_nk_cyc4_mca.png", plot=fig_vill_d_bl_nk_cyc4_mca, path = figure_path,
       width = 20, height = 15, units = "cm")
ggsave("vill_ud_bl_k_cyc4_mca.png", plot=fig_vill_ud_bl_k_cyc4_mca, path = figure_path,
       width = 20, height = 15, units = "cm")
ggsave("vill_ud_bl_nk_cyc4_mca.png", plot=fig_vill_ud_bl_nk_cyc4_mca, path = figure_path,
       width = 20, height = 15, units = "cm")
ggsave("vill_d_fr_k_cyc4_mca.png", plot=fig_vill_d_fr_k_cyc4_mca, path = figure_path,
       width = 20, height = 15, units = "cm")
ggsave("vill_d_fr_nk_cyc4_mca.png", plot=fig_vill_d_fr_nk_cyc4_mca, path = figure_path,
       width = 20, height = 15, units = "cm")


# Check model assumptions
check_model(lm_cycpath4_d_bl_k_df_mca)
check_model(lm_cycpath4_d_bl_nk_df_mca)
check_model(lm_cycpath4_ud_bl_k_df_mca)
check_model(lm_cycpath4_ud_bl_nk_df_mca)
check_model(lm_cycpath4_d_fr_k_df_mca)
check_model(lm_cycpath4_d_fr_nk_df_mca)

summary(lm_cycpath4_d_bl_k_df_mca)
summary(lm_cycpath4_d_bl_nk_df_mca)
summary(lm_cycpath4_ud_bl_k_df_mca)
summary(lm_cycpath4_ud_bl_nk_df_mca)
summary(lm_cycpath4_d_fr_k_df_mca)
summary(lm_cycpath4_d_fr_nk_df_mca)




# # food insecurity
# lm_cycpath2_food_insc_sds <- lmer(food_insc ~ cycdens2_sds + pathdens2_sds + time +
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_food_insc_sds <- lmer(food_insc ~ cycdens3_sds + pathdens3_sds + time + 
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_food_insc_sds <- lmer(food_insc ~ cycdens4_sds + pathdens4_sds + time + 
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_food_insc_sds <- lmer(food_insc ~ cycdens5_sds + pathdens5_sds + time +
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# lm_cycpath2_sev_food_insc_sds <- lmer(severe_food_insc ~ cycdens2_sds + pathdens2_sds + time +
#                                     (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_sev_food_insc_sds <- lmer(severe_food_insc ~ cycdens3_sds + pathdens3_sds + time + 
#                                          (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_sev_food_insc_sds <- lmer(severe_food_insc ~ cycdens4_sds + pathdens4_sds + time + 
#                                          (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_sev_food_insc_sds <- lmer(severe_food_insc ~ cycdens5_sds + pathdens5_sds + time +
#                                          (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# ## HH attributes
# # electricity
# lm_cycpath2_elctr_pct_sds <- lmer(elctr_pct ~ cycdens2_sds + pathdens2_sds + time +
#                                     (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_elctr_pct_sds <- lmer(elctr_pct ~ cycdens3_sds + pathdens3_sds + time + 
#                                     (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_elctr_pct_sds <- lmer(elctr_pct ~ cycdens4_sds + pathdens4_sds + time + 
#                                     (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_elctr_pct_sds <- lmer(elctr_pct ~ cycdens5_sds + pathdens5_sds + time +
#                                     (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# # tv
# lm_cycpath2_tv_pct_sds <- lmer(tv_pct ~ cycdens2_sds + pathdens2_sds + time +
#                                     (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_tv_pct_sds <- lmer(tv_pct ~ cycdens3_sds + pathdens3_sds + time + 
#                                     (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_tv_pct_sds <- lmer(tv_pct ~ cycdens4_sds + pathdens4_sds + time + 
#                                     (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_tv_pct_sds <- lmer(tv_pct ~ cycdens5_sds + pathdens5_sds + time +
#                                     (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# # cell / mobile phone
# lm_cycpath2_cell_pct_sds <- lmer(cell_pct ~ cycdens2_sds + pathdens2_sds + time +
#                                  (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_cell_pct_sds <- lmer(cell_pct ~ cycdens3_sds + pathdens3_sds + time + 
#                                  (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_cell_pct_sds <- lmer(cell_pct ~ cycdens4_sds + pathdens4_sds + time + 
#                                  (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_cell_pct_sds <- lmer(cell_pct ~ cycdens5_sds + pathdens5_sds + time +
#                                  (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# # no phone
# lm_cycpath2_no_phone_sds <- lmer(no_phone ~ cycdens2_sds + pathdens2_sds + time +
#                                    (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_no_phone_sds <- lmer(no_phone ~ cycdens3_sds + pathdens3_sds + time + 
#                                    (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_no_phone_sds <- lmer(no_phone ~ cycdens4_sds + pathdens4_sds + time + 
#                                    (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_no_phone_sds <- lmer(no_phone ~ cycdens5_sds + pathdens5_sds + time +
#                                    (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# # no electronics
# lm_cycpath2_no_electronics_sds <- lmer(no_electronics ~ cycdens2_sds + pathdens2_sds + time +
#                                    (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_no_electronics_sds <- lmer(no_electronics ~ cycdens3_sds + pathdens3_sds + time + 
#                                    (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_no_electronics_sds <- lmer(no_electronics ~ cycdens4_sds + pathdens4_sds + time + 
#                                    (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_no_electronics_sds <- lmer(no_electronics ~ cycdens5_sds + pathdens5_sds + time +
#                                    (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# # flush toilet
# lm_cycpath2_flush_tlt_pct_sds <- lmer(flush_tlt_pct ~ cycdens2_sds + pathdens2_sds + time +
#                                  (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_flush_tlt_pct_sds <- lmer(flush_tlt_pct ~ cycdens3_sds + pathdens3_sds + time + 
#                                  (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_flush_tlt_pct_sds <- lmer(flush_tlt_pct ~ cycdens4_sds + pathdens4_sds + time + 
#                                  (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_flush_tlt_pct_sds <- lmer(flush_tlt_pct ~ cycdens5_sds + pathdens5_sds + time +
#                                  (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# # outdoor toilet
# lm_cycpath2_outdoor_tlt_pct_sds <- lmer(outdoor_tlt_pct ~ cycdens2_sds + pathdens2_sds + time +
#                                         (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_outdoor_tlt_pct_sds <- lmer(outdoor_tlt_pct ~ cycdens3_sds + pathdens3_sds + time + 
#                                         (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_outdoor_tlt_pct_sds <- lmer(outdoor_tlt_pct ~ cycdens4_sds + pathdens4_sds + time + 
#                                         (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_outdoor_tlt_pct_sds <- lmer(outdoor_tlt_pct ~ cycdens5_sds + pathdens5_sds + time +
#                                         (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# # separate kitchen
# lm_cycpath2_sep_kitch_pct_sds <- lmer(sep_kitch_pct ~ cycdens2_sds + pathdens2_sds + time +
#                                           (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_sep_kitch_pct_sds <- lmer(sep_kitch_pct ~ cycdens3_sds + pathdens3_sds + time + 
#                                           (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_sep_kitch_pct_sds <- lmer(sep_kitch_pct ~ cycdens4_sds + pathdens4_sds + time + 
#                                           (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_sep_kitch_pct_sds <- lmer(sep_kitch_pct ~ cycdens5_sds + pathdens5_sds + time +
#                                           (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# # stove with chimney
# lm_cycpath2_stove_w_chmny_pct_sds <- lmer(stove_w_chmny_pct ~ cycdens2_sds + pathdens2_sds + time +
#                                         (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_stove_w_chmny_pct_sds <- lmer(stove_w_chmny_pct ~ cycdens3_sds + pathdens3_sds + time + 
#                                         (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_stove_w_chmny_pct_sds <- lmer(stove_w_chmny_pct ~ cycdens4_sds + pathdens4_sds + time + 
#                                         (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_stove_w_chmny_pct_sds <- lmer(stove_w_chmny_pct ~ cycdens5_sds + pathdens5_sds + time +
#                                         (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# # stove without chimney
# lm_cycpath2_stove_no_chmny_pct_sds <- lmer(stove_no_chmny_pct ~ cycdens2_sds + pathdens2_sds + time +
#                                         (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_stove_no_chmny_pct_sds <- lmer(stove_no_chmny_pct ~ cycdens3_sds + pathdens3_sds + time + 
#                                         (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_stove_no_chmny_pct_sds <- lmer(stove_no_chmny_pct ~ cycdens4_sds + pathdens4_sds + time + 
#                                         (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_stove_no_chmny_pct_sds <- lmer(stove_no_chmny_pct ~ cycdens5_sds + pathdens5_sds + time +
#                                         (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# # earth / sand floor
# lm_cycpath2_earth_flr_pct_sds <- lmer(earth_flr_pct ~ cycdens2_sds + pathdens2_sds + time +
#                                              (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_earth_flr_pct_sds <- lmer(earth_flr_pct ~ cycdens3_sds + pathdens3_sds + time + 
#                                              (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_earth_flr_pct_sds <- lmer(earth_flr_pct ~ cycdens4_sds + pathdens4_sds + time + 
#                                              (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_earth_flr_pct_sds <- lmer(earth_flr_pct ~ cycdens5_sds + pathdens5_sds + time +
#                                              (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# # ceramic floor
# lm_cycpath2_ceramic_flr_pct_sds <- lmer(ceramic_flr_pct ~ cycdens2_sds + pathdens2_sds + time +
#                                              (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_ceramic_flr_pct_sds <- lmer(ceramic_flr_pct ~ cycdens3_sds + pathdens3_sds + time + 
#                                              (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_ceramic_flr_pct_sds <- lmer(ceramic_flr_pct ~ cycdens4_sds + pathdens4_sds + time + 
#                                              (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_ceramic_flr_pct_sds <- lmer(ceramic_flr_pct ~ cycdens5_sds + pathdens5_sds + time +
#                                              (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# # no windows
# lm_cycpath2_no_window_pct_sds <- lmer(no_window_pct ~ cycdens2_sds + pathdens2_sds + time +
#                                           (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_no_window_pct_sds <- lmer(no_window_pct ~ cycdens3_sds + pathdens3_sds + time + 
#                                           (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_no_window_pct_sds <- lmer(no_window_pct ~ cycdens4_sds + pathdens4_sds + time + 
#                                           (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_no_window_pct_sds <- lmer(no_window_pct ~ cycdens5_sds + pathdens5_sds + time +
#                                           (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# # cement walls
# lm_cycpath2_cement_walls_pct_sds <- lmer(cement_walls_pct ~ cycdens2_sds + pathdens2_sds + time +
#                                          (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_cement_walls_pct_sds <- lmer(cement_walls_pct ~ cycdens3_sds + pathdens3_sds + time + 
#                                          (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_cement_walls_pct_sds <- lmer(cement_walls_pct ~ cycdens4_sds + pathdens4_sds + time + 
#                                          (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_cement_walls_pct_sds <- lmer(cement_walls_pct ~ cycdens5_sds + pathdens5_sds + time +
#                                          (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# # sleeping rooms
# lm_cycpath2_over_1_sleeprm_pct_sds <- lmer(over_1_sleeprm_pct ~ cycdens2_sds + pathdens2_sds + time +
#                                           (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath3_over_1_sleeprm_pct_sds <- lmer(over_1_sleeprm_pct ~ cycdens3_sds + pathdens3_sds + time + 
#                                           (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath4_over_1_sleeprm_pct_sds <- lmer(over_1_sleeprm_pct ~ cycdens4_sds + pathdens4_sds + time + 
#                                           (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# lm_cycpath5_over_1_sleeprm_pct_sds <- lmer(over_1_sleeprm_pct ~ cycdens5_sds + pathdens5_sds + time +
#                                           (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# # Label wealth variables
# var_label(wealth_corr) <- list(
#   village_wealth_index = "Village wealth (avg hh quintile)",
#   vill_wealth_mca = "Village wealth (avg MCA coord)",
#   food_insc = "% food insecure",
#   severe_food_insc = "% severely food insecure",
#   elctr_pct = "% with electricity",
#   tv_pct = "% with TV",
#   cell_pct = "% with cell/mobile phone", 
#   no_phone = "% with no phone (cell or landline)",
#   no_electronics = "% with no electronics",
#   flush_tlt_pct = "% with flush toilet",
#   outdoor_tlt_pct = "% using outdoors for toilet facility",
#   sep_kitch_pct = "% with separate kitchen room",
#   stove_w_chmny_pct = "% have a stove with chimney", 
#   stove_no_chmny_pct = "% have a stove but no chimney",
#   ceramic_flr_pct = "% with ceramic floors",
#   earth_flr_pct = "% with earth/sand floors",
#   no_window_pct = "% with no window at home", 
#   cement_walls_pct = "% with cement walls",
#   over_1_sleeprm_pct = "% with >1 room used for sleeping"
# )
# # Create heatmap of results
# # Pull coefficients from the mixed models
# cycdens2_coef <- c(coef(summary(lm_cycpath2_vill_wealth_sds))[2,1], coef(summary(lm_cycpath2_vill_wealth_mca_sds))[2,1],
#                    coef(summary(lm_cycpath2_food_insc_sds))[2,1], coef(summary(lm_cycpath2_sev_food_insc_sds))[2,1],
#                    coef(summary(lm_cycpath2_elctr_pct_sds))[2,1], coef(summary(lm_cycpath2_tv_pct_sds))[2,1],
#                    coef(summary(lm_cycpath2_cell_pct_sds))[2,1], coef(summary(lm_cycpath2_no_phone_sds))[2,1],
#                    coef(summary(lm_cycpath2_no_electronics_sds))[2,1], coef(summary(lm_cycpath2_flush_tlt_pct_sds))[2,1],
#                    coef(summary(lm_cycpath2_outdoor_tlt_pct_sds))[2,1], coef(summary(lm_cycpath2_sep_kitch_pct_sds))[2,1],
#                    coef(summary(lm_cycpath2_stove_w_chmny_pct_sds))[2,1], coef(summary(lm_cycpath2_stove_no_chmny_pct_sds))[2,1],
#                    coef(summary(lm_cycpath2_ceramic_flr_pct_sds))[2,1], coef(summary(lm_cycpath2_earth_flr_pct_sds))[2,1],
#                    coef(summary(lm_cycpath2_no_window_pct_sds))[2,1], coef(summary(lm_cycpath2_cement_walls_pct_sds))[2,1],
#                    coef(summary(lm_cycpath2_over_1_sleeprm_pct_sds))[2,1])
# cycdens3_coef <- c(coef(summary(lm_cycpath3_vill_wealth_sds))[2,1], coef(summary(lm_cycpath3_vill_wealth_mca_sds))[2,1],
#                    coef(summary(lm_cycpath3_food_insc_sds))[2,1], coef(summary(lm_cycpath3_sev_food_insc_sds))[2,1],
#                    coef(summary(lm_cycpath3_elctr_pct_sds))[2,1], coef(summary(lm_cycpath3_tv_pct_sds))[2,1],
#                    coef(summary(lm_cycpath3_cell_pct_sds))[2,1], coef(summary(lm_cycpath3_no_phone_sds))[2,1],
#                    coef(summary(lm_cycpath3_no_electronics_sds))[2,1], coef(summary(lm_cycpath3_flush_tlt_pct_sds))[2,1],
#                    coef(summary(lm_cycpath3_outdoor_tlt_pct_sds))[2,1], coef(summary(lm_cycpath3_sep_kitch_pct_sds))[2,1],
#                    coef(summary(lm_cycpath3_stove_w_chmny_pct_sds))[2,1], coef(summary(lm_cycpath3_stove_no_chmny_pct_sds))[2,1],
#                    coef(summary(lm_cycpath3_ceramic_flr_pct_sds))[2,1], coef(summary(lm_cycpath3_earth_flr_pct_sds))[2,1],
#                    coef(summary(lm_cycpath3_no_window_pct_sds))[2,1], coef(summary(lm_cycpath3_cement_walls_pct_sds))[2,1],
#                    coef(summary(lm_cycpath3_over_1_sleeprm_pct_sds))[2,1])
# cycdens4_coef <- c(coef(summary(lm_cycpath4_vill_wealth_sds))[2,1], coef(summary(lm_cycpath4_vill_wealth_mca_sds))[2,1],
#                    coef(summary(lm_cycpath4_food_insc_sds))[2,1], coef(summary(lm_cycpath4_sev_food_insc_sds))[2,1],
#                    coef(summary(lm_cycpath4_elctr_pct_sds))[2,1], coef(summary(lm_cycpath4_tv_pct_sds))[2,1],
#                    coef(summary(lm_cycpath4_cell_pct_sds))[2,1], coef(summary(lm_cycpath4_no_phone_sds))[2,1],
#                    coef(summary(lm_cycpath4_no_electronics_sds))[2,1], coef(summary(lm_cycpath4_flush_tlt_pct_sds))[2,1],
#                    coef(summary(lm_cycpath4_outdoor_tlt_pct_sds))[2,1], coef(summary(lm_cycpath4_sep_kitch_pct_sds))[2,1],
#                    coef(summary(lm_cycpath4_stove_w_chmny_pct_sds))[2,1], coef(summary(lm_cycpath4_stove_no_chmny_pct_sds))[2,1],
#                    coef(summary(lm_cycpath4_ceramic_flr_pct_sds))[2,1], coef(summary(lm_cycpath4_earth_flr_pct_sds))[2,1],
#                    coef(summary(lm_cycpath4_no_window_pct_sds))[2,1], coef(summary(lm_cycpath4_cement_walls_pct_sds))[2,1],
#                    coef(summary(lm_cycpath4_over_1_sleeprm_pct_sds))[2,1])
# cycdens5_coef <- c(coef(summary(lm_cycpath5_vill_wealth_sds))[2,1], coef(summary(lm_cycpath5_vill_wealth_mca_sds))[2,1],
#                    coef(summary(lm_cycpath5_food_insc_sds))[2,1], coef(summary(lm_cycpath5_sev_food_insc_sds))[2,1],
#                    coef(summary(lm_cycpath5_elctr_pct_sds))[2,1], coef(summary(lm_cycpath5_tv_pct_sds))[2,1],
#                    coef(summary(lm_cycpath5_cell_pct_sds))[2,1], coef(summary(lm_cycpath5_no_phone_sds))[2,1],
#                    coef(summary(lm_cycpath5_no_electronics_sds))[2,1], coef(summary(lm_cycpath5_flush_tlt_pct_sds))[2,1],
#                    coef(summary(lm_cycpath5_outdoor_tlt_pct_sds))[2,1], coef(summary(lm_cycpath5_sep_kitch_pct_sds))[2,1],
#                    coef(summary(lm_cycpath5_stove_w_chmny_pct_sds))[2,1], coef(summary(lm_cycpath5_stove_no_chmny_pct_sds))[2,1],
#                    coef(summary(lm_cycpath5_ceramic_flr_pct_sds))[2,1], coef(summary(lm_cycpath5_earth_flr_pct_sds))[2,1],
#                    coef(summary(lm_cycpath5_no_window_pct_sds))[2,1], coef(summary(lm_cycpath5_cement_walls_pct_sds))[2,1],
#                    coef(summary(lm_cycpath5_over_1_sleeprm_pct_sds))[2,1])
# 
# # Pull t-values from models and convert to p-values (approximate using z distr)
# cycdens2_sig <- c(coef(summary(lm_cycpath2_vill_wealth_sds))[2,3], coef(summary(lm_cycpath2_vill_wealth_mca_sds))[2,3],
#                    coef(summary(lm_cycpath2_food_insc_sds))[2,3], coef(summary(lm_cycpath2_sev_food_insc_sds))[2,3],
#                    coef(summary(lm_cycpath2_elctr_pct_sds))[2,3], coef(summary(lm_cycpath2_tv_pct_sds))[2,3],
#                    coef(summary(lm_cycpath2_cell_pct_sds))[2,3], coef(summary(lm_cycpath2_no_phone_sds))[2,3],
#                    coef(summary(lm_cycpath2_no_electronics_sds))[2,3], coef(summary(lm_cycpath2_flush_tlt_pct_sds))[2,3],
#                    coef(summary(lm_cycpath2_outdoor_tlt_pct_sds))[2,3], coef(summary(lm_cycpath2_sep_kitch_pct_sds))[2,3],
#                    coef(summary(lm_cycpath2_stove_w_chmny_pct_sds))[2,3], coef(summary(lm_cycpath2_stove_no_chmny_pct_sds))[2,3],
#                    coef(summary(lm_cycpath2_ceramic_flr_pct_sds))[2,3], coef(summary(lm_cycpath2_earth_flr_pct_sds))[2,3],
#                    coef(summary(lm_cycpath2_no_window_pct_sds))[2,3], coef(summary(lm_cycpath2_cement_walls_pct_sds))[2,3],
#                    coef(summary(lm_cycpath2_over_1_sleeprm_pct_sds))[2,3])
# cycdens3_sig <- c(coef(summary(lm_cycpath3_vill_wealth_sds))[2,3], coef(summary(lm_cycpath3_vill_wealth_mca_sds))[2,3],
#                    coef(summary(lm_cycpath3_food_insc_sds))[2,3], coef(summary(lm_cycpath3_sev_food_insc_sds))[2,3],
#                    coef(summary(lm_cycpath3_elctr_pct_sds))[2,3], coef(summary(lm_cycpath3_tv_pct_sds))[2,3],
#                    coef(summary(lm_cycpath3_cell_pct_sds))[2,3], coef(summary(lm_cycpath3_no_phone_sds))[2,3],
#                    coef(summary(lm_cycpath3_no_electronics_sds))[2,3], coef(summary(lm_cycpath3_flush_tlt_pct_sds))[2,3],
#                    coef(summary(lm_cycpath3_outdoor_tlt_pct_sds))[2,3], coef(summary(lm_cycpath3_sep_kitch_pct_sds))[2,3],
#                    coef(summary(lm_cycpath3_stove_w_chmny_pct_sds))[2,3], coef(summary(lm_cycpath3_stove_no_chmny_pct_sds))[2,3],
#                    coef(summary(lm_cycpath3_ceramic_flr_pct_sds))[2,3], coef(summary(lm_cycpath3_earth_flr_pct_sds))[2,3],
#                    coef(summary(lm_cycpath3_no_window_pct_sds))[2,3], coef(summary(lm_cycpath3_cement_walls_pct_sds))[2,3],
#                    coef(summary(lm_cycpath3_over_1_sleeprm_pct_sds))[2,3])
# cycdens4_sig <- c(coef(summary(lm_cycpath4_vill_wealth_sds))[2,3], coef(summary(lm_cycpath4_vill_wealth_mca_sds))[2,3],
#                    coef(summary(lm_cycpath4_food_insc_sds))[2,3], coef(summary(lm_cycpath4_sev_food_insc_sds))[2,3],
#                    coef(summary(lm_cycpath4_elctr_pct_sds))[2,3], coef(summary(lm_cycpath4_tv_pct_sds))[2,3],
#                    coef(summary(lm_cycpath4_cell_pct_sds))[2,3], coef(summary(lm_cycpath4_no_phone_sds))[2,3],
#                    coef(summary(lm_cycpath4_no_electronics_sds))[2,3], coef(summary(lm_cycpath4_flush_tlt_pct_sds))[2,3],
#                    coef(summary(lm_cycpath4_outdoor_tlt_pct_sds))[2,3], coef(summary(lm_cycpath4_sep_kitch_pct_sds))[2,3],
#                    coef(summary(lm_cycpath4_stove_w_chmny_pct_sds))[2,3], coef(summary(lm_cycpath4_stove_no_chmny_pct_sds))[2,3],
#                    coef(summary(lm_cycpath4_ceramic_flr_pct_sds))[2,3], coef(summary(lm_cycpath4_earth_flr_pct_sds))[2,3],
#                    coef(summary(lm_cycpath4_no_window_pct_sds))[2,3], coef(summary(lm_cycpath4_cement_walls_pct_sds))[2,3],
#                    coef(summary(lm_cycpath4_over_1_sleeprm_pct_sds))[2,3])
# cycdens5_sig <- c(coef(summary(lm_cycpath5_vill_wealth_sds))[2,3], coef(summary(lm_cycpath5_vill_wealth_mca_sds))[2,3],
#                    coef(summary(lm_cycpath5_food_insc_sds))[2,3], coef(summary(lm_cycpath5_sev_food_insc_sds))[2,3],
#                    coef(summary(lm_cycpath5_elctr_pct_sds))[2,3], coef(summary(lm_cycpath5_tv_pct_sds))[2,3],
#                    coef(summary(lm_cycpath5_cell_pct_sds))[2,3], coef(summary(lm_cycpath5_no_phone_sds))[2,3],
#                    coef(summary(lm_cycpath5_no_electronics_sds))[2,3], coef(summary(lm_cycpath5_flush_tlt_pct_sds))[2,3],
#                    coef(summary(lm_cycpath5_outdoor_tlt_pct_sds))[2,3], coef(summary(lm_cycpath5_sep_kitch_pct_sds))[2,3],
#                    coef(summary(lm_cycpath5_stove_w_chmny_pct_sds))[2,3], coef(summary(lm_cycpath5_stove_no_chmny_pct_sds))[2,3],
#                    coef(summary(lm_cycpath5_ceramic_flr_pct_sds))[2,3], coef(summary(lm_cycpath5_earth_flr_pct_sds))[2,3],
#                    coef(summary(lm_cycpath5_no_window_pct_sds))[2,3], coef(summary(lm_cycpath5_cement_walls_pct_sds))[2,3],
#                    coef(summary(lm_cycpath5_over_1_sleeprm_pct_sds))[2,3])
# 
# vill_wealth_heatmap_sig_df <- data.frame(cycdens2_tv = cycdens2_sig, cycdens3_tv = cycdens3_sig, 
#                                           cycdens4_tv = cycdens4_sig, cycdens5_tv = cycdens5_sig) %>%
#   # Convert t-value to p-value using z distribution as an approximation
#   # Correct for multiple testing
#   mutate(cycdens2_pv = p.adjust(2 * (1 - pnorm(abs(cycdens2_tv))), method="BH"),
#          cycdens3_pv = p.adjust(2 * (1 - pnorm(abs(cycdens3_tv))), method="BH"),
#          cycdens4_pv = p.adjust(2 * (1 - pnorm(abs(cycdens4_tv))), method="BH"),
#          cycdens5_pv = p.adjust(2 * (1 - pnorm(abs(cycdens5_tv))), method="BH")) %>%
#   select(cycdens2_pv, cycdens3_pv, cycdens4_pv, cycdens5_pv)
# 
# vill_wealth_heatmap_coef_df <- data.frame(cycdens2 = cycdens2_coef, cycdens3 = cycdens3_coef, 
#                                           cycdens4 = cycdens4_coef, cycdens5 = cycdens5_coef)
# colnames(vill_wealth_heatmap_coef_df) <- c("Cyclic density 2", "Cyclic density 3", "Cyclic density 4",
#                                            "Cyclic density 5")
# rownames(vill_wealth_heatmap_coef_df) <- c("Village wealth (avg hh quintile)", "Village wealth (avg MCA coord)", "Food insecure",
#                                            "Severely food insecure","Electricity","TV","Cell/mobile phone", 
#                                            "No phone (cell or landline)","No electronics","Flush toilet",
#                                            "Outdoors as toilet facility","Separate kitchen room","Stove with chimney", 
#                                            "Stove but no chimney","Ceramic floors","Earth/sand floors",
#                                            "No window at home", "Cement walls",">1 room used for sleeping")
# 
# # Color scale specifications (matching Shiv's)
# myb2<-c(0.0000000125,0.0012500000,0.0065105438,0.0112500000,0.0143750000,0.0179000000,0.0242018138,0.0520843541,0.1)
# myBreaks<-c(myb2[c(length(myb2):-1:1)]*-1,myb2)
# myColor<-colorRampPalette(c("#cb181d","white", "#33a02c"))(17)
# 
# vill_wealth_heatmap <- pheatmap(vill_wealth_heatmap_coef_df, annotation_row = NULL, annotation_names_row = T,
#          labels_row = rownames(vill_wealth_heatmap_coef_df), legend = T,
#          fontsize_number = 15, border_color = "#EEEEEE", na_col = "white",
#          fontsize_col = 10, angle_col = 45, fontsize_row = 10, #fontface="bold",
#          color=myColor, breaks=myBreaks, #display_numbers = t(disp_fdr),
#          number_color = "black", treeheight_col = 0, treeheight_row = 0, cluster_rows = F,
#          cluster_cols = F)
# ggsave("vill_wealth_heatmap_undir_borrow_lend_no_kin.png", plot=vill_wealth_heatmap, path = figure_path,
#        width = 20, height = 15, units = "cm")
# 
# # Sig increase in cycdens corresponds to inc vill wealth using wealth quintiles
# # (cycdens2, 3, 4 significant)
# summ(lm_cycpath2_vill_wealth_sds, digits=3)
# summ(lm_cycpath3_vill_wealth_sds, digits=3)
# summ(lm_cycpath4_vill_wealth_sds, digits=3)
# summ(lm_cycpath5_vill_wealth_sds, digits=3)
# 
# # Sig increase in cycdens corresponds to inc vill wealth using MCA avg
# # (cycdens2, 3, 4 significant)
# summ(lm_cycpath2_vill_wealth_mca_sds, digits=3)
# summ(lm_cycpath3_vill_wealth_mca_sds, digits=3)
# summ(lm_cycpath4_vill_wealth_mca_sds, digits=3)
# summ(lm_cycpath5_vill_wealth_mca_sds, digits=3)
# 
# # Increase in cyc density, sig decrease in food insecurity (cycdens 3, 4)
# summ(lm_cycpath2_food_insc_sds, digits=3)
# summ(lm_cycpath3_food_insc_sds, digits=3)
# summ(lm_cycpath4_food_insc_sds, digits=3)
# summ(lm_cycpath5_food_insc_sds, digits=3)
# 
# # sig inc in cycdens 3 correspons with dec severe food insecurity
# summ(lm_cycpath2_sev_food_insc_sds, digits=3)
# summ(lm_cycpath3_sev_food_insc_sds, digits=3)
# summ(lm_cycpath4_sev_food_insc_sds, digits=3)
# summ(lm_cycpath5_sev_food_insc_sds, digits=3)
# 
# # electricity: only cycdens2 sig (pos corr)
# summ(lm_cycpath2_elctr_pct_sds, digits=3)
# summ(lm_cycpath3_elctr_pct_sds, digits=3)
# summ(lm_cycpath4_elctr_pct_sds, digits=3)
# summ(lm_cycpath5_elctr_pct_sds, digits=3)
# 
# # tv: cycdens2, 3, 4 sig (pos corr)
# summ(lm_cycpath2_tv_pct_sds, digits=3)
# summ(lm_cycpath3_tv_pct_sds, digits=3)
# summ(lm_cycpath4_tv_pct_sds, digits=3)
# summ(lm_cycpath5_tv_pct_sds, digits=3)
# 
# # cell/mobile phone: inc significantly (cycdens 2)
# summ(lm_cycpath2_cell_pct_sds, digits=3)
# summ(lm_cycpath3_cell_pct_sds, digits=3)
# summ(lm_cycpath4_cell_pct_sds, digits=3)
# summ(lm_cycpath5_cell_pct_sds, digits=3)
# 
# # no phone:
# summ(lm_cycpath2_no_phone_sds, digits=3)
# summ(lm_cycpath3_no_phone_sds, digits=3)
# summ(lm_cycpath4_no_phone_sds, digits=3)
# summ(lm_cycpath5_no_phone_sds, digits=3)
# 
# # no electronics:
# summ(lm_cycpath2_no_electronics_sds, digits=3)
# summ(lm_cycpath3_no_electronics_sds, digits=3)
# summ(lm_cycpath4_no_electronics_sds, digits=3)
# summ(lm_cycpath5_no_electronics_sds, digits=3)
# 
# # flush toilet: only cycdens2 sig (pos cor)
# summ(lm_cycpath2_flush_tlt_pct_sds, digits=3)
# summ(lm_cycpath3_flush_tlt_pct_sds, digits=3)
# summ(lm_cycpath4_flush_tlt_pct_sds, digits=3)
# summ(lm_cycpath5_flush_tlt_pct_sds, digits=3)
# 
# # outdoor toilet: no sig corr
# summ(lm_cycpath2_outdoor_tlt_pct_sds, digits=3)
# summ(lm_cycpath3_outdoor_tlt_pct_sds, digits=3)
# summ(lm_cycpath4_outdoor_tlt_pct_sds, digits=3)
# summ(lm_cycpath5_outdoor_tlt_pct_sds, digits=3)
# 
# # separate kitchen: cycdens and time sig neg corr; pathdens sig pos corr; 
# summ(lm_cycpath2_sep_kitch_pct_sds, digits=3)
# summ(lm_cycpath3_sep_kitch_pct_sds, digits=3)
# summ(lm_cycpath4_sep_kitch_pct_sds, digits=3)
# summ(lm_cycpath5_sep_kitch_pct_sds, digits=3)
# 
# # stove with chimney: cycdens2 sig neg corr; time sig pos corr; pathdens5 sig neg corr
# summ(lm_cycpath2_stove_w_chmny_pct_sds, digits=3)
# summ(lm_cycpath3_stove_w_chmny_pct_sds, digits=3)
# summ(lm_cycpath4_stove_w_chmny_pct_sds, digits=3)
# summ(lm_cycpath5_stove_w_chmny_pct_sds, digits=3)
# 
# # stove without chimney: time sig neg corr
# summ(lm_cycpath2_stove_no_chmny_pct_sds, digits=3)
# summ(lm_cycpath3_stove_no_chmny_pct_sds, digits=3)
# summ(lm_cycpath4_stove_no_chmny_pct_sds, digits=3)
# summ(lm_cycpath5_stove_no_chmny_pct_sds, digits=3)
# 
# # ceramic floor: pathdens sig neg corr
# summ(lm_cycpath2_ceramic_flr_pct_sds, digits=3)
# summ(lm_cycpath3_ceramic_flr_pct_sds, digits=3)
# summ(lm_cycpath4_ceramic_flr_pct_sds, digits=3)
# summ(lm_cycpath5_ceramic_flr_pct_sds, digits=3)
# 
# # earth / sand floor: cycdens and time sig neg corr
# summ(lm_cycpath2_earth_flr_pct_sds, digits=3)
# summ(lm_cycpath3_earth_flr_pct_sds, digits=3)
# summ(lm_cycpath4_earth_flr_pct_sds, digits=3)
# summ(lm_cycpath5_earth_flr_pct_sds, digits=3)
# 
# # no windows: no sig corr (only time sig neg corr for cycles2)
# summ(lm_cycpath2_no_window_pct_sds, digits=3)
# summ(lm_cycpath3_no_window_pct_sds, digits=3)
# summ(lm_cycpath4_no_window_pct_sds, digits=3)
# summ(lm_cycpath5_no_window_pct_sds, digits=3)
# 
# # cement walls: no sig corr
# summ(lm_cycpath2_cement_walls_pct_sds, digits=3)
# summ(lm_cycpath3_cement_walls_pct_sds, digits=3)
# summ(lm_cycpath4_cement_walls_pct_sds, digits=3)
# summ(lm_cycpath5_cement_walls_pct_sds, digits=3)
# 
# # sleeping rooms: cycdens2,3 4, 5 sig (pos cor)
# summ(lm_cycpath2_over_1_sleeprm_pct_sds, digits=3)
# summ(lm_cycpath3_over_1_sleeprm_pct_sds, digits=3)
# summ(lm_cycpath4_over_1_sleeprm_pct_sds, digits=3)
# summ(lm_cycpath5_over_1_sleeprm_pct_sds, digits=3)
# 
# # Models separating cycles and paths for descriptive comparison of coefficients
# ### Raw cycles and paths (separate models)
# lm_cyc2_w1w3_sds <- lmer(food_insc ~ kcycle2 + num_nodes + num_edges + time + 
#                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # <0.001 sig
# lm_cyc3_w1w3_sds <- lmer(food_insc ~ kcycle3 + num_nodes + num_edges + time + 
#                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.003 sig
# lm_cyc4_w1w3_sds <- lmer(food_insc ~ kcycle4 + num_nodes + num_edges + time + 
#                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.004 significance
# lm_cyc5_w1w3_sds <- lmer(food_insc ~ kcycle5 + num_nodes + num_edges + time + 
#                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.08 significance
# 
# lm_path2_w1w3_sds <- lmer(food_insc ~ kpath2 + num_nodes + num_edges + time + 
#                             (1|village_code), data = kstats_borrow_lend_dir_kin_df) # <0.001 sig
# lm_path3_w1w3_sds <- lmer(food_insc ~ kpath3 + num_nodes + num_edges + time + 
#                             (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.003 sig
# lm_path4_w1w3_sds <- lmer(food_insc ~ kpath4 + num_nodes + num_edges + time + 
#                             (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.004 significance
# lm_path5_w1w3_sds <- lmer(food_insc ~ kpath5 + num_nodes + num_edges + time + 
#                             (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# summ(lm_cyc2_w1w3_sds, digits=3)
# summ(lm_cyc3_w1w3_sds, digits=3)
# summ(lm_cyc4_w1w3_sds, digits=3)
# summ(lm_cyc5_w1w3_sds, digits=3)
# 
# summ(lm_path2_w1w3_sds, digits=3)
# summ(lm_path3_w1w3_sds, digits=3)
# summ(lm_path4_w1w3_sds, digits=3)
# summ(lm_path5_w1w3_sds, digits=3)
# 
# 
# ### Cyclic and path density (separate models)
# lm_cycdens2_w1w3_sds <- lmer(food_insc ~ cycdens2_sds + time + num_nodes + num_edges +
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # <0.001 sig
# lm_cycdens3_w1w3_sds <- lmer(food_insc ~ cycdens3_sds + time + num_nodes + num_edges +
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.003 sig
# lm_cycdens4_w1w3_sds <- lmer(food_insc ~ cycdens4_sds + time + num_nodes + num_edges +
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.004 significance
# lm_cycdens5_w1w3_sds <- lmer(food_insc ~ cycdens5_sds + time + num_nodes + num_edges +
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.08 significance
# 
# lm_pathdens2_w1w3_sds <- lmer(food_insc ~ pathdens2_sds + time + num_nodes + num_edges +
#                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # <0.001 sig
# lm_pathdens3_w1w3_sds <- lmer(food_insc ~ pathdens3_sds + time + num_nodes + num_edges +
#                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.003 sig
# lm_pathdens4_w1w3_sds <- lmer(food_insc ~ pathdens4_sds + time + num_nodes + num_edges +
#                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.004 significance
# lm_pathdens5_w1w3_sds <- lmer(food_insc ~ pathdens5_sds + time + num_nodes + num_edges +
#                            (1|village_code), data = kstats_borrow_lend_dir_kin_df)
# 
# summ(lm_cycdens2_w1w3_sds, digits=3)
# summ(lm_cycdens3_w1w3_sds, digits=3)
# summ(lm_cycdens4_w1w3_sds, digits=3)
# summ(lm_cycdens5_w1w3_sds, digits=3)
# 
# summ(lm_pathdens2_w1w3_sds, digits=3)
# summ(lm_pathdens3_w1w3_sds, digits=3)
# summ(lm_pathdens4_w1w3_sds, digits=3)
# summ(lm_pathdens5_w1w3_sds, digits=3)
# 
# # Including cycles and paths in the same model
# ### Raw cycle and path counts
# lm_cycpath2_w1w3_sds <- lmer(food_insc ~ kcycle2 + kpath2 + num_nodes + num_edges + time + 
#                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # <0.001 sig
# lm_cycpath3_w1w3_sds <- lmer(food_insc ~ kcycle3 + kpath3 + num_nodes + num_edges + time + 
#                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.003 sig
# lm_cycpath4_w1w3_sds <- lmer(food_insc ~ kcycle4 + kpath4 + num_nodes + num_edges + time + 
#                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.004 significance
# lm_cycpath5_w1w3_sds <- lmer(food_insc ~ kcycle5 + kpath5 + num_nodes + num_edges + time + 
#                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.08 significance
# 
# summ(lm_cycpath2_w1w3_sds, digits=3)
# summ(lm_cycpath3_w1w3_sds, digits=3)
# summ(lm_cycpath4_w1w3_sds, digits=3)
# summ(lm_cycpath5_w1w3_sds, digits=3)
# 
# # # Including both cyclic and path densities in the model (same degree distr for random networks)
# lm_cycpath2_borrow_lend_w1w3_sds <- lmer(food_insc ~ cycdens2_sds + pathdens2_sds + time +
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # <0.001 sig
# lm_cycpath3_borrow_lend_w1w3_sds <- lmer(food_insc ~ cycdens3_sds + pathdens3_sds + time + 
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.003 sig
# lm_cycpath4_borrow_lend_w1w3_sds <- lmer(food_insc ~ cycdens4_sds + pathdens4_sds + time + 
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.004 significance
# lm_cycpath5_borrow_lend_w1w3_sds <- lmer(food_insc ~ cycdens5_sds + pathdens5_sds + time +
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.08 significance
# 
# # lm_cycpath2_hh_borrow_lend_w1w3_sds <- lmer(village_wealth_index ~ cycdens2_sds + pathdens2_sds + time + (1|village_code), data = kstats_hh_borrow_lend_df) # <0.001 sig
# # lm_cycpath3_hh_borrow_lend_w1w3_sds <- lmer(village_wealth_index ~ cycdens3_sds + pathdens3_sds + time + (1|village_code), data = kstats_hh_borrow_lend_df) # <0.001 sig
# # lm_cycpath4_hh_borrow_lend_w1w3_sds <- lmer(village_wealth_index ~ cycdens4_sds + pathdens4_sds + time + (1|village_code), data = kstats_hh_borrow_lend_df)
# # lm_cycpath5_hh_borrow_lend_w1w3_sds <- lmer(village_wealth_index ~ cycdens5_sds + pathdens5_sds + time + (1|village_code), data = kstats_hh_borrow_lend_df)
# 
# # Same degree distribution
# summ(lm_cycpath2_borrow_lend_w1w3_sds, digits=3)
# summ(lm_cycpath3_borrow_lend_w1w3_sds, digits=3)
# summ(lm_cycpath4_borrow_lend_w1w3_sds, digits=3)
# summ(lm_cycpath5_borrow_lend_w1w3_sds, digits=3)
# 
# # summ(lm_cycpath2_hh_borrow_lend_w1w3_sds, digits=3)
# # summ(lm_cycpath3_hh_borrow_lend_w1w3_sds, digits=3)
# # summ(lm_cycpath4_hh_borrow_lend_w1w3_sds, digits=3)
# # summ(lm_cycpath5_hh_borrow_lend_w1w3_sds, digits=3)
# 
# 
# #-------- Check hh and combined borrow/lend network cycles w vill food insecurity (waves 1 and 3) ----------
# # Using erdos-renyi random graphs
# # Using the same degree distribution as observed for random graphs
# lm_cycpath2_borrow_lend_w1w3_sds <- lmer(food_insc ~ cycdens2_sds + pathdens2_sds + time + 
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # <0.001 sig
# lm_cycpath3_borrow_lend_w1w3_sds <- lmer(food_insc ~ cycdens3_sds + pathdens3_sds + time + 
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.003 sig
# lm_cycpath4_borrow_lend_w1w3_sds <- lmer(food_insc ~ cycdens4_sds + pathdens4_sds + time + 
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.004 significance
# lm_cycpath5_borrow_lend_w1w3_sds <- lmer(food_insc ~ cycdens5_sds + pathdens5_sds + time + 
#                                            (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.08 significance
# 
# # Same degree distribution
# summ(lm_cycpath2_borrow_lend_w1w3_sds, digits=3)
# summ(lm_cycpath3_borrow_lend_w1w3_sds, digits=3)
# summ(lm_cycpath4_borrow_lend_w1w3_sds, digits=3)
# summ(lm_cycpath5_borrow_lend_w1w3_sds, digits=3)
# 
# 
# 
# #-------- Check hh and combined borrow/lend network cycles w vill gini (waves 1 and 3) ----------
# # Using the same degree distribution as observed for random graphs
# lm_cycpath2_borrow_lend_gini_w1w3_sds <- lmer(gini ~ cycdens2_sds + pathdens2_sds + time + (1|village_code), data = kstats_borrow_lend_dir_kin_df) # <0.001 sig
# lm_cycpath2_hh_borrow_lend_gini_w1w3_sds <- lmer(gini ~ cycdens2_sds + pathdens2_sds + time + (1|village_code), data = kstats_hh_borrow_lend_df) # <0.001 sig
# 
# lm_cycpath3_borrow_lend_gini_w1w3_sds <- lmer(gini ~ cycdens3_sds + pathdens3_sds + time + (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.003 sig
# lm_cycpath3_hh_borrow_lend_gini_w1w3_sds <- lmer(gini ~ cycdens3_sds + pathdens3_sds + time + (1|village_code), data = kstats_hh_borrow_lend_df) # <0.001 sig
# 
# lm_cycpath4_borrow_lend_gini_w1w3_sds <- lmer(gini ~ cycdens4_sds + pathdens4_sds + time + (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.004 significance
# lm_cycpath4_hh_borrow_lend_gini_w1w3_sds <- lmer(gini ~ cycdens4_sds + pathdens4_sds + time + (1|village_code), data = kstats_hh_borrow_lend_df)
# 
# lm_cycpath5_borrow_lend_gini_w1w3_sds <- lmer(gini ~ cycdens5_sds + pathdens5_sds + time + (1|village_code), data = kstats_borrow_lend_dir_kin_df) # 0.08 significance
# lm_cycpath5_hh_borrow_lend_gini_w1w3_sds <- lmer(gini ~ cycdens5_sds + pathdens5_sds + time + (1|village_code), data = kstats_hh_borrow_lend_df)
# 
# # Same degree distribution
# summ(lm_cycpath2_borrow_lend_gini_w1w3_sds, digits=3)
# summ(lm_cycpath3_borrow_lend_gini_w1w3_sds, digits=3)
# summ(lm_cycpath4_borrow_lend_gini_w1w3_sds, digits=3)
# summ(lm_cycpath5_borrow_lend_gini_w1w3_sds, digits=3)
# 
# summ(lm_cycpath2_hh_borrow_lend_gini_w1w3_sds, digits=3)
# summ(lm_cycpath3_hh_borrow_lend_gini_w1w3_sds, digits=3)
# summ(lm_cycpath4_hh_borrow_lend_gini_w1w3_sds, digits=3)
# summ(lm_cycpath5_hh_borrow_lend_gini_w1w3_sds, digits=3)
# 
# 
# 
# 
# 
# 
# 
# #-------- Check hh and combined borrow/lend network cycles w vill wealth (wave 3) ----------
# indiv_df_w3 <- subset(kstats_borrow_lend_dir_kin_df, time==2)
# hh_df_w3 <- subset(kstats_hh_borrow_lend_df, time==2)
# 
# # Using erdos-renyi rndom graphs
# lm_cycpath2_borrow_lend_w3_er <- lm(village_wealth_index ~ cycdens2_er + pathdens2_er, data = indiv_df_w3) # <0.001 sig
# lm_cycpath2_hh_borrow_lend_w3_er <- lm(village_wealth_index ~ cycdens2_er + pathdens2_er, data = hh_df_w3) # <0.001 sig
# 
# lm_cycpath3_borrow_lend_w3_er <- lm(village_wealth_index ~ cycdens3_er + pathdens3_er, data = indiv_df_w3) # 0.003 sig
# lm_cycpath3_hh_borrow_lend_w3_er <- lm(village_wealth_index ~ cycdens3_er + pathdens3_er, data = hh_df_w3) # <0.001 sig
# 
# lm_cycpath4_borrow_lend_w3_er <- lm(village_wealth_index ~ cycdens4_er + pathdens4_er, data = indiv_df_w3) # 0.004 significance
# lm_cycpath4_hh_borrow_lend_w3_er <- lm(village_wealth_index ~ cycdens4_er + pathdens4_er, data = hh_df_w3)
# 
# lm_cycpath5_borrow_lend_w3_er <- lm(village_wealth_index ~ cycdens5_er + pathdens5_er, data = indiv_df_w3) # 0.08 significance
# lm_cycpath5_hh_borrow_lend_w3_er <- lm(village_wealth_index ~ cycdens5_er + pathdens5_er, data = hh_df_w3)
# 
# # Using barabasi-albert random graphs
# lm_cycpath2_borrow_lend_w3_ba <- lm(village_wealth_index ~ cycdens2_ba + pathdens2_ba, data = indiv_df_w3) # <0.001 sig
# lm_cycpath2_hh_borrow_lend_w3_ba <- lm(village_wealth_index ~ cycdens2_ba + pathdens2_ba, data = hh_df_w3) # <0.001 sig
# 
# lm_cycpath3_borrow_lend_w3_ba <- lm(village_wealth_index ~ cycdens3_ba + pathdens3_ba, data = indiv_df_w3) # 0.003 sig
# lm_cycpath3_hh_borrow_lend_w3_ba <- lm(village_wealth_index ~ cycdens3_ba + pathdens3_ba, data = hh_df_w3) # <0.001 sig
# 
# lm_cycpath4_borrow_lend_w3_ba <- lm(village_wealth_index ~ cycdens4_ba + pathdens4_ba, data = indiv_df_w3) # 0.004 significance
# lm_cycpath4_hh_borrow_lend_w3_ba <- lm(village_wealth_index ~ cycdens4_ba + pathdens4_ba, data = hh_df_w3)
# 
# lm_cycpath5_borrow_lend_w3_ba <- lm(village_wealth_index ~ cycdens5_ba + pathdens5_ba, data = indiv_df_w3) # 0.08 significance
# lm_cycpath5_hh_borrow_lend_w3_ba <- lm(village_wealth_index ~ cycdens5_ba + pathdens5_ba, data = hh_df_w3)
# 
# # Using the power law distribution for random graphs
# lm_cycpath2_borrow_lend_w3_pl <- lm(village_wealth_index ~ cycdens2_pl + pathdens2_pl, data = indiv_df_w3) # <0.001 sig
# lm_cycpath2_hh_borrow_lend_w3_pl <- lm(village_wealth_index ~ cycdens2_pl + pathdens2_pl, data = hh_df_w3) # <0.001 sig
# 
# lm_cycpath3_borrow_lend_w3_pl <- lm(village_wealth_index ~ cycdens3_pl + pathdens3_pl, data = indiv_df_w3) # 0.003 sig
# lm_cycpath3_hh_borrow_lend_w3_pl <- lm(village_wealth_index ~ cycdens3_pl + pathdens3_pl, data = hh_df_w3) # <0.001 sig
# 
# lm_cycpath4_borrow_lend_w3_pl <- lm(village_wealth_index ~ cycdens4_pl + pathdens4_pl, data = indiv_df_w3) # 0.004 significance
# lm_cycpath4_hh_borrow_lend_w3_pl <- lm(village_wealth_index ~ cycdens4_pl + pathdens4_pl, data = hh_df_w3)
# 
# lm_cycpath5_borrow_lend_w3_pl <- lm(village_wealth_index ~ cycdens5_pl + pathdens5_pl, data = indiv_df_w3) # 0.08 significance
# lm_cycpath5_hh_borrow_lend_w3_pl <- lm(village_wealth_index ~ cycdens5_pl + pathdens5_pl, data = hh_df_w3)
# 
# # Using the same degree distribution as observed for random graphs
# lm_cycpath2_borrow_lend_w3_sds <- lm(village_wealth_index ~ cycdens2_sds + pathdens2_sds, data = indiv_df_w3) # <0.001 sig
# lm_cycpath2_hh_borrow_lend_w3_sds <- lm(village_wealth_index ~ cycdens2_sds + pathdens2_sds, data = hh_df_w3) # <0.001 sig
# 
# lm_cycpath3_borrow_lend_w3_sds <- lm(village_wealth_index ~ cycdens3_sds + pathdens3_sds, data = indiv_df_w3) # 0.003 sig
# lm_cycpath3_hh_borrow_lend_w3_sds <- lm(village_wealth_index ~ cycdens3_sds + pathdens3_sds, data = hh_df_w3) # <0.001 sig
# 
# lm_cycpath4_borrow_lend_w3_sds <- lm(village_wealth_index ~ cycdens4_sds + pathdens4_sds, data = indiv_df_w3) # 0.004 significance
# lm_cycpath4_hh_borrow_lend_w3_sds <- lm(village_wealth_index ~ cycdens4_sds + pathdens4_sds, data = hh_df_w3)
# 
# lm_cycpath5_borrow_lend_w3_sds <- lm(village_wealth_index ~ cycdens5_sds + pathdens5_sds, data = indiv_df_w3) # 0.08 significance
# lm_cycpath5_hh_borrow_lend_w3_sds <- lm(village_wealth_index ~ cycdens5_sds + pathdens5_sds, data = hh_df_w3)
# 
# # ER graphs
# summary(lm_cycpath2_borrow_lend_w3_er)
# summary(lm_cycpath2_hh_borrow_lend_w3_er)
# 
# summary(lm_cycpath3_borrow_lend_w3_er)
# summary(lm_cycpath3_hh_borrow_lend_w3_er)
# 
# summary(lm_cycpath4_borrow_lend_w3_er)
# summary(lm_cycpath4_hh_borrow_lend_w3_er)
# 
# summary(lm_cycpath5_borrow_lend_w3_er)
# summary(lm_cycpath5_hh_borrow_lend_w3_er)
# 
# # Barabasi-albert graphs
# summary(lm_cycpath2_borrow_lend_w3_ba)
# summary(lm_cycpath2_hh_borrow_lend_w3_ba)
# 
# summary(lm_cycpath3_borrow_lend_w3_ba)
# summary(lm_cycpath3_hh_borrow_lend_w3_ba)
# 
# summary(lm_cycpath4_borrow_lend_w3_ba)
# summary(lm_cycpath4_hh_borrow_lend_w3_ba)
# 
# summary(lm_cycpath5_borrow_lend_w3_ba)
# summary(lm_cycpath5_hh_borrow_lend_w3_ba)
# 
# # Power law distribution
# summary(lm_cycpath2_borrow_lend_w3_pl)
# summary(lm_cycpath2_hh_borrow_lend_w3_pl)
# 
# summary(lm_cycpath3_borrow_lend_w3_pl)
# summary(lm_cycpath3_hh_borrow_lend_w3_pl)
# 
# summary(lm_cycpath4_borrow_lend_w3_pl)
# summary(lm_cycpath4_hh_borrow_lend_w3_pl)
# 
# summary(lm_cycpath5_borrow_lend_w3_pl)
# summary(lm_cycpath5_hh_borrow_lend_w3_pl)
# 
# # Same degree distribution
# summary(lm_cycpath2_borrow_lend_w3_sds)
# summary(lm_cycpath2_hh_borrow_lend_w3_sds)
# 
# summary(lm_cycpath3_borrow_lend_w3_sds)
# summary(lm_cycpath3_hh_borrow_lend_w3_sds)
# 
# summary(lm_cycpath4_borrow_lend_w3_sds)
# summary(lm_cycpath4_hh_borrow_lend_w3_sds)
# 
# summary(lm_cycpath5_borrow_lend_w3_sds)
# summary(lm_cycpath5_hh_borrow_lend_w3_sds)
# 
# 
# #-------- Check hh and combined borrow/lend network cycles w vill wealth ----------
# indiv_df_w1 <- subset(kstats_borrow_lend_dir_kin_df, time==1)
# hh_df_w1 <- subset(kstats_hh_borrow_lend_df, time==1)
# 
# # Using erdos-renyi rndom graphs
# lm_cycpath2_borrow_lend_er <- lm(village_wealth_index ~ cycdens2_er + pathdens2_er, data = indiv_df_w1) # <0.001 sig
# lm_cycpath2_hh_borrow_lend_er <- lm(village_wealth_index ~ cycdens2_er + pathdens2_er, data = hh_df_w1) # <0.001 sig
# 
# lm_cycpath3_borrow_lend_er <- lm(village_wealth_index ~ cycdens3_er + pathdens3_er, data = indiv_df_w1) # 0.003 sig
# lm_cycpath3_hh_borrow_lend_er <- lm(village_wealth_index ~ cycdens3_er + pathdens3_er, data = hh_df_w1) # <0.001 sig
# 
# lm_cycpath4_borrow_lend_er <- lm(village_wealth_index ~ cycdens4_er + pathdens4_er, data = indiv_df_w1) # 0.004 significance
# lm_cycpath4_hh_borrow_lend_er <- lm(village_wealth_index ~ cycdens4_er + pathdens4_er, data = hh_df_w1)
# 
# lm_cycpath5_borrow_lend_er <- lm(village_wealth_index ~ cycdens5_er + pathdens5_er, data = indiv_df_w1) # 0.08 significance
# lm_cycpath5_hh_borrow_lend_er <- lm(village_wealth_index ~ cycdens5_er + pathdens5_er, data = hh_df_w1)
# 
# # Using barabasi-albert random graphs
# lm_cycpath2_borrow_lend_ba <- lm(village_wealth_index ~ cycdens2_ba + pathdens2_ba, data = indiv_df_w1) # <0.001 sig
# lm_cycpath2_hh_borrow_lend_ba <- lm(village_wealth_index ~ cycdens2_ba + pathdens2_ba, data = hh_df_w1) # <0.001 sig
# 
# lm_cycpath3_borrow_lend_ba <- lm(village_wealth_index ~ cycdens3_ba + pathdens3_ba, data = indiv_df_w1) # 0.003 sig
# lm_cycpath3_hh_borrow_lend_ba <- lm(village_wealth_index ~ cycdens3_ba + pathdens3_ba, data = hh_df_w1) # <0.001 sig
# 
# lm_cycpath4_borrow_lend_ba <- lm(village_wealth_index ~ cycdens4_ba + pathdens4_ba, data = indiv_df_w1) # 0.004 significance
# lm_cycpath4_hh_borrow_lend_ba <- lm(village_wealth_index ~ cycdens4_ba + pathdens4_ba, data = hh_df_w1)
# 
# lm_cycpath5_borrow_lend_ba <- lm(village_wealth_index ~ cycdens5_ba + pathdens5_ba, data = indiv_df_w1) # 0.08 significance
# lm_cycpath5_hh_borrow_lend_ba <- lm(village_wealth_index ~ cycdens5_ba + pathdens5_ba, data = hh_df_w1)
# 
# # Using the power law distribution for random graphs
# lm_cycpath2_borrow_lend_pl <- lm(village_wealth_index ~ cycdens2_pl + pathdens2_pl, data = indiv_df_w1) # <0.001 sig
# lm_cycpath2_hh_borrow_lend_pl <- lm(village_wealth_index ~ cycdens2_pl + pathdens2_pl, data = hh_df_w1) # <0.001 sig
# 
# lm_cycpath3_borrow_lend_pl <- lm(village_wealth_index ~ cycdens3_pl + pathdens3_pl, data = indiv_df_w1) # 0.003 sig
# lm_cycpath3_hh_borrow_lend_pl <- lm(village_wealth_index ~ cycdens3_pl + pathdens3_pl, data = hh_df_w1) # <0.001 sig
# 
# lm_cycpath4_borrow_lend_pl <- lm(village_wealth_index ~ cycdens4_pl + pathdens4_pl, data = indiv_df_w1) # 0.004 significance
# lm_cycpath4_hh_borrow_lend_pl <- lm(village_wealth_index ~ cycdens4_pl + pathdens4_pl, data = hh_df_w1)
# 
# lm_cycpath5_borrow_lend_pl <- lm(village_wealth_index ~ cycdens5_pl + pathdens5_pl, data = indiv_df_w1) # 0.08 significance
# lm_cycpath5_hh_borrow_lend_pl <- lm(village_wealth_index ~ cycdens5_pl + pathdens5_pl, data = hh_df_w1)
# 
# # Using the same degree distribution as observed for random graphs
# lm_cycpath2_borrow_lend_sds <- lm(village_wealth_index ~ cycdens2_sds + pathdens2_sds, data = indiv_df_w1) # <0.001 sig
# lm_cycpath2_hh_borrow_lend_sds <- lm(village_wealth_index ~ cycdens2_sds + pathdens2_sds, data = hh_df_w1) # <0.001 sig
# 
# lm_cycpath3_borrow_lend_sds <- lm(village_wealth_index ~ cycdens3_sds + pathdens3_sds, data = indiv_df_w1) # 0.003 sig
# lm_cycpath3_hh_borrow_lend_sds <- lm(village_wealth_index ~ cycdens3_sds + pathdens3_sds, data = hh_df_w1) # <0.001 sig
# 
# lm_cycpath4_borrow_lend_sds <- lm(village_wealth_index ~ cycdens4_sds + pathdens4_sds, data = indiv_df_w1) # 0.004 significance
# lm_cycpath4_hh_borrow_lend_sds <- lm(village_wealth_index ~ cycdens4_sds + pathdens4_sds, data = hh_df_w1)
# 
# lm_cycpath5_borrow_lend_sds <- lm(village_wealth_index ~ cycdens5_sds + pathdens5_sds, data = indiv_df_w1) # 0.08 significance
# lm_cycpath5_hh_borrow_lend_sds <- lm(village_wealth_index ~ cycdens5_sds + pathdens5_sds, data = hh_df_w1)
# 
# # ER graphs
# summary(lm_cycpath2_borrow_lend_er)
# summary(lm_cycpath2_hh_borrow_lend_er)
# 
# summary(lm_cycpath3_borrow_lend_er)
# summary(lm_cycpath3_hh_borrow_lend_er)
# 
# summary(lm_cycpath4_borrow_lend_er)
# summary(lm_cycpath4_hh_borrow_lend_er)
# 
# summary(lm_cycpath5_borrow_lend_er)
# summary(lm_cycpath5_hh_borrow_lend_er)
# 
# # Barabasi-albert graphs
# summary(lm_cycpath2_borrow_lend_ba)
# summary(lm_cycpath2_hh_borrow_lend_ba)
# 
# summary(lm_cycpath3_borrow_lend_ba)
# summary(lm_cycpath3_hh_borrow_lend_ba)
# 
# summary(lm_cycpath4_borrow_lend_ba)
# summary(lm_cycpath4_hh_borrow_lend_ba)
# 
# summary(lm_cycpath5_borrow_lend_ba)
# summary(lm_cycpath5_hh_borrow_lend_ba)
# 
# # Power law distribution
# summary(lm_cycpath2_borrow_lend_pl)
# summary(lm_cycpath2_hh_borrow_lend_pl)
# 
# summary(lm_cycpath3_borrow_lend_pl)
# summary(lm_cycpath3_hh_borrow_lend_pl)
# 
# summary(lm_cycpath4_borrow_lend_pl)
# summary(lm_cycpath4_hh_borrow_lend_pl)
# 
# summary(lm_cycpath5_borrow_lend_pl)
# summary(lm_cycpath5_hh_borrow_lend_pl)
# 
# # Same degree distribution
# summary(lm_cycpath2_borrow_lend_sds)
# summary(lm_cycpath2_hh_borrow_lend_sds)
# 
# summary(lm_cycpath3_borrow_lend_sds)
# summary(lm_cycpath3_hh_borrow_lend_sds)
# 
# summary(lm_cycpath4_borrow_lend_sds)
# summary(lm_cycpath4_hh_borrow_lend_sds)
# 
# summary(lm_cycpath5_borrow_lend_sds)
# summary(lm_cycpath5_hh_borrow_lend_sds)


#------------------- Vill-level correlation btw cycles and wealth -----------------
###### Look at cyclic densities by network type
kstats_df <- kstats_d_bl_k_df %>%
  select(village_code, cycdens2_sds, cycdens3_sds, cycdens4_sds, cycdens5_sds) %>%
  add_column(network_type = "d_bl_k",
             network_name = "Directed, borrow-lend, kin incl") %>%
  bind_rows(kstats_d_bl_nk_df %>%
              select(village_code, cycdens2_sds, cycdens3_sds, cycdens4_sds, cycdens5_sds) %>%
              add_column(network_type = "d_bl_nk",
                         network_name = "Directed, borrow-lend, kin excl")) %>%
  bind_rows(kstats_ud_bl_nk_df %>%
              select(village_code, cycdens2_sds, cycdens3_sds, cycdens4_sds, cycdens5_sds) %>%
              add_column(network_type = "ud_bl_k",
                         network_name = "Undirected, borrow-lend, kin incl")) %>%
  bind_rows(kstats_ud_bl_nk_df %>%
              select(village_code, cycdens2_sds, cycdens3_sds, cycdens4_sds, cycdens5_sds) %>%
              add_column(network_type = "ud_bl_nk",
                         network_name = "Undirected, borrow-lend, kin excl")) %>%
  bind_rows(kstats_d_fr_k_df %>%
              select(village_code, cycdens2_sds, cycdens3_sds, cycdens4_sds, cycdens5_sds) %>%
              add_column(network_type = "d_fr_k",
                         network_name = "Directed, friend, kin incl")) %>%
  bind_rows(kstats_d_fr_nk_df %>%
              select(village_code, cycdens2_sds, cycdens3_sds, cycdens4_sds, cycdens5_sds) %>%
              add_column(network_type = "d_fr_nk",
                         network_name = "Directed, friend, kin excl")) %>%
  bind_rows(kstats_ud_fr_k_df %>%
              select(village_code, cycdens2_sds, cycdens3_sds, cycdens4_sds, cycdens5_sds) %>%
              add_column(network_type = "ud_fr_k",
                         network_name = "Undirected, friend, kin incl")) %>%
  bind_rows(kstats_ud_fr_nk_df %>%
              select(village_code, cycdens2_sds, cycdens3_sds, cycdens4_sds, cycdens5_sds) %>%
              add_column(network_type = "ud_fr_nk",
                         network_name = "Undirected, friend, kin excl"))

# Villages with above expected number of cycles
vill_ids_kcycle2 <- kstats_df %>%
  subset(cycdens2_sds > 0)
table(vill_ids_kcycle2$network_type)

vill_ids_kcycle3 <- kstats_df %>%
  subset(cycdens3_sds > 0)
table(vill_ids_kcycle3$network_type)

vill_ids_kcycle4 <- kstats_df %>%
  subset(cycdens4_sds > 0)
table(vill_ids_kcycle4$network_type)

vill_ids_kcycle5 <- kstats_df %>%
  subset(cycdens5_sds > 0)
table(vill_ids_kcycle5$network_type)


fig_cycdens2_distr <- ggplot(kstats_df, aes(x=cycdens2_sds, fill=network_name)) +
  geom_density(alpha=0.4) +
  geom_vline(xintercept = 0, color="red", linetype="dashed") +
  theme_classic()
fig_cycdens2_distr

fig_cycdens3_distr <- ggplot(kstats_df, aes(x=cycdens3_sds, fill=network_name)) +
  geom_density(alpha=0.4) +
  geom_vline(xintercept = 0, color="red", linetype="dashed") +
  theme_classic()
fig_cycdens3_distr

fig_cycdens4_distr <- ggplot(kstats_df, aes(x=cycdens4_sds, fill=network_name)) +
  geom_density(alpha=0.4) +
  geom_vline(xintercept = 0, color="red", linetype="dashed") +
  theme_classic()
fig_cycdens4_distr

fig_cycdens5_distr <- ggplot(kstats_df, aes(x=cycdens5_sds, fill=network_name)) +
  geom_density(alpha=0.4) +
  geom_vline(xintercept = 0, color="red", linetype="dashed") +
  theme_classic()
fig_cycdens5_distr

ggarrange(fig_cycdens2_distr, fig_cycdens3_distr, fig_cycdens4_distr,
          fig_cycdens5_distr, common.legend=TRUE, ncol = 2, nrow = 2)


# 71 villages have at least 1 kcycle of length 4
vill_ids_kcycle4 <- kstats_df %>%
  subset(cycdens4_sds > 0) #%>%
table(vill_ids_kcycle4$network_type)

# 35 villages have at least 1 kcycle of length 5
vill_ids_kcycle5 <- vill_df %>%
  subset(kcycle5 > 0) %>%
  pull(village_code)

# 33 villages with cyclic density 4 above zero
vill_df %>%
  subset(cycdens_resp4 > 0) %>%
  pull(village_code) %>%
  length()

vill_df_kcycle4 <- vill_df %>%
  subset(kcycle4 > 0)

# # Run linear model to see if log(cyclic density) is correlated with village wealth at w1 or w3
# # NOTE: NaN's produced for linear models including cyclic density for cycles of length 4 and 5
# # (for linear models looking at w1 or w3 village wealth)
# lm_log_cycdens2_wealth_w1 <- lm(village_wealth_index_w1 ~ log1p(cycdens_resp2), data = vill_df) # significant (+)
# lm_log_cycdens3_wealth_w1 <- lm(village_wealth_index_w1 ~ log1p(cycdens_resp3), data = vill_df) # significant (+)
# lm_log_cycdens4_wealth_w1 <- lm(village_wealth_index_w1 ~ log1p(cycdens_resp4), data = vill_df)
# lm_log_cycdens5_wealth_w1 <- lm(village_wealth_index_w1 ~ log1p(cycdens_resp5), data = vill_df)
# 
# lm_log_cycdens2_wealth_w3 <- lm(village_wealth_index_w3 ~ log1p(cycdens_resp2), data = vill_df) # significant (+)
# lm_log_cycdens3_wealth_w3 <- lm(village_wealth_index_w3 ~ log1p(cycdens_resp3), data = vill_df) # significant (+)
# lm_log_cycdens4_wealth_w3 <- lm(village_wealth_index_w3 ~ log1p(cycdens_resp4), data = vill_df)
# lm_log_cycdens5_wealth_w3 <- lm(village_wealth_index_w3 ~ log1p(cycdens_resp5), data = vill_df)
# 
# 
# # Run linear model to see if cyclic density is correlated with village wealth at w1 or w3
# # for subset of villages that have at least one kcycle of length 4 present.
# lm_cycdens2_wealth_w1_kcyc4sub <- lm(village_wealth_index_w1 ~ cycdens_resp2, data = vill_df_kcycle4) # significant (+)
# lm_cycdens3_wealth_w1_kcyc4sub <- lm(village_wealth_index_w1 ~ cycdens_resp3, data = vill_df_kcycle4)
# lm_cycdens4_wealth_w1_kcyc4sub <- lm(village_wealth_index_w1 ~ cycdens_resp4, data = vill_df_kcycle4)
# lm_cycdens5_wealth_w1_kcyc4sub <- lm(village_wealth_index_w1 ~ cycdens_resp5, data = vill_df_kcycle4)
# 
# lm_cycdens2_wealth_w3_kcyc4sub <- lm(village_wealth_index_w3 ~ cycdens_resp2, data = vill_df_kcycle4) # significant (+)
# lm_cycdens3_wealth_w3_kcyc4sub <- lm(village_wealth_index_w3 ~ cycdens_resp3, data = vill_df_kcycle4)
# lm_cycdens4_wealth_w3_kcyc4sub <- lm(village_wealth_index_w3 ~ cycdens_resp4, data = vill_df_kcycle4)
# lm_cycdens5_wealth_w3_kcyc4sub <- lm(village_wealth_index_w3 ~ cycdens_resp5, data = vill_df_kcycle4)
# 
# # Run linear model to see if cyclic density is correlated with village wealth at w1 or w3
# lm_cycdens2_wealth_w1 <- lm(village_wealth_index_w1 ~ cycdens_resp2, data = vill_df) # significant (+)
# lm_cycdens3_wealth_w1 <- lm(village_wealth_index_w1 ~ cycdens_resp3, data = vill_df) # significant (+)
# lm_cycdens4_wealth_w1 <- lm(village_wealth_index_w1 ~ cycdens_resp4, data = vill_df)
# lm_cycdens5_wealth_w1 <- lm(village_wealth_index_w1 ~ cycdens_resp5, data = vill_df)
# 
# lm_cycdens2_wealth_w3 <- lm(village_wealth_index_w3 ~ cycdens_resp2, data = vill_df) # significant
# lm_cycdens3_wealth_w3 <- lm(village_wealth_index_w3 ~ cycdens_resp3, data = vill_df) # significant
# lm_cycdens4_wealth_w3 <- lm(village_wealth_index_w3 ~ cycdens_resp4, data = vill_df)
# lm_cycdens5_wealth_w3 <- lm(village_wealth_index_w3 ~ cycdens_resp5, data = vill_df)
# 
# ggplot(vill_df, aes(x=cycdens_resp5, y=village_wealth_index_w3)) + 
#   geom_point() + 
#   xlab("Cyclic density (kcycle5)") + 
#   ylab("Village wealth w3 (mean)") +
#   theme(legend.position = "none",
#         panel.background = element_blank(), 
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         axis.title.x=element_text(size=20),
#         axis.title.y=element_text(size=20),
#         axis.text.x=element_text(size=15,color="black"),
#         axis.text.y=element_text(size=15,color="black"),
#         axis.line.y = element_line(color="black"),
#         axis.line.x = element_line(color="black"))
# 
# # Run linear model for cyclic density and difference in village wealth between w3 and w1
# lm_cycdens2_wealth_diff_w3_w1 <- lm(wealth_diff_w3_w1 ~ cycdens_resp2, data = vill_df) # significant
# lm_cycdens3_wealth_diff_w3_w1 <- lm(wealth_diff_w3_w1 ~ cycdens_resp3, data = vill_df)
# lm_cycdens4_wealth_diff_w3_w1 <- lm(wealth_diff_w3_w1 ~ cycdens_resp4, data = vill_df)
# lm_cycdens5_wealth_diff_w3_w1 <- lm(wealth_diff_w3_w1 ~ cycdens_resp5, data = vill_df)
# 
# # Check reciprocity, transitivity, and density alone against village wealth and wealth diff
# lm_reciprocity_wealth_w1 <- lm(village_wealth_index_w1 ~ reciprocity.borrow, data = vill_df)
# lm_reciprocity_wealth_w3 <- lm(village_wealth_index_w3 ~ reciprocity.borrow, data = vill_df)
# lm_reciprocity_wealth_diff_w3_w1 <- lm(wealth_diff_w3_w1 ~ reciprocity.borrow, data = vill_df)
# 
# lm_transitivity_wealth_w1 <- lm(village_wealth_index_w1 ~ transitivity.borrow, data = vill_df)
# lm_transitivity_wealth_w3 <- lm(village_wealth_index_w3 ~ transitivity.borrow, data = vill_df)
# lm_transitivity_wealth_diff_w3_w1 <- lm(wealth_diff_w3_w1 ~ transitivity.borrow, data = vill_df)
# 
# lm_density_wealth_w1 <- lm(village_wealth_index_w1 ~ density.borrow, data = vill_df)
# lm_density_wealth_w3 <- lm(village_wealth_index_w3 ~ density.borrow, data = vill_df)
# lm_density_wealth_diff_w3_w1 <- lm(wealth_diff_w3_w1 ~ density.borrow, data = vill_df)
# 
# 
# ###
# 
# # Run linear model to see if presence of kcycle is correlated with village wealth at w1 or w3
# lm_kcycle2_wealth_w1 <- lm(village_wealth_index_w1 ~ kcycle2 + num_nodes.borrow + num_edges.borrow, data = vill_df)
# lm_kcycle3_wealth_w1 <- lm(village_wealth_index_w1 ~ kcycle3 + num_nodes.borrow + num_edges.borrow, data = vill_df)
# lm_kcycle4_wealth_w1 <- lm(village_wealth_index_w1 ~ kcycle4 + num_nodes.borrow + num_edges.borrow, data = vill_df)
# lm_kcycle5_wealth_w1 <- lm(village_wealth_index_w1 ~ kcycle5 + num_nodes.borrow + num_edges.borrow, data = vill_df)
# 
# lm_kcycle2_wealth_w3 <- lm(village_wealth_index_w3 ~ kcycle2 + num_nodes.borrow + num_edges.borrow, data = vill_df)
# lm_kcycle3_wealth_w3 <- lm(village_wealth_index_w3 ~ kcycle3 + num_nodes.borrow + num_edges.borrow, data = vill_df)
# lm_kcycle4_wealth_w3 <- lm(village_wealth_index_w3 ~ kcycle4 + num_nodes.borrow + num_edges.borrow, data = vill_df)
# lm_kcycle5_wealth_w3 <- lm(village_wealth_index_w3 ~ kcycle5 + num_nodes.borrow + num_edges.borrow, data = vill_df)
# 
# # Look at linear regression for kcycle prevalence and difference in village wealth between w3 and w1
# lm_kcycle2_wealth_diff_w3_w1 <- lm(wealth_diff_w3_w1 ~ kcycle2 + num_nodes.borrow + num_edges.borrow, data = vill_df)
# lm_kcycle3_wealth_diff_w3_w1 <- lm(wealth_diff_w3_w1 ~ kcycle3 + num_nodes.borrow + num_edges.borrow, data = vill_df)
# lm_kcycle4_wealth_diff_w3_w1 <- lm(wealth_diff_w3_w1 ~ kcycle4 + num_nodes.borrow + num_edges.borrow, data = vill_df)
# lm_kcycle5_wealth_diff_w3_w1 <- lm(wealth_diff_w3_w1 ~ kcycle5 + num_nodes.borrow + num_edges.borrow, data = vill_df)
# 
# # Run linear model to see if presence of kpath of length 4 is correlated with village wealth at w1 or w3
# lm_kpath4_wealth_w1 <- lm(village_wealth_index_w1 ~ kpath4 + vill_size, data = vill_df)
# 
# lm_kpath5_wealth_w3 <- lm(village_wealth_index_w3 ~ kpath5 + vill_size, data = vill_df)
# lm_kpath4_wealth_w3 <- lm(village_wealth_index_w3 ~ kpath4 + vill_size, data = vill_df)
# lm_kpath3_wealth_w3 <- lm(village_wealth_index_w3 ~ kpath3 + vill_size, data = vill_df)
# lm_kpath2_wealth_w3 <- lm(village_wealth_index_w3 ~ kpath2 + vill_size, data = vill_df)
# 
# # Look at linear regression for kcycle prevalence and difference in village wealth between w3 and w1
# lm_kpath5_wealth_diff_w3_w1 <- lm(wealth_diff_w3_w1 ~ kpath5 + vill_size, data = vill_df)
# lm_kpath4_wealth_diff_w3_w1 <- lm(wealth_diff_w3_w1 ~ kpath4 + vill_size, data = vill_df)
# lm_kpath3_wealth_diff_w3_w1 <- lm(wealth_diff_w3_w1 ~ kpath3 + vill_size, data = vill_df)
# lm_kpath2_wealth_diff_w3_w1 <- lm(wealth_diff_w3_w1 ~ kpath2 + vill_size, data = vill_df)
# 
# 
# # Look at logistic regression for kcycle prevalence and difference in village wealth between w3 and w1
# lgr_kcycle5_wealth_diff.factor <- glm(wealth_diff_w3_w1.factor ~ kcycle5 + vill_size, data = vill_df,
#                                       family = "binomial")
# 
# # Look at linear regression for standard network measures for friend_family name generators and borrowing
# # name generators in predicting village wealth difference between w3 and w1.
# lm_friend_family_transitivity_wealth_diff <- lm(wealth_diff_w3_w1 ~ transitivity.friend_family + 
#                                                   vill_size, data = vill_df)
# lm_friend_family_density_wealth_diff <- lm(wealth_diff_w3_w1 ~ density.friend_family + 
#                                              vill_size, data = vill_df)
# 
# lm_borrow_transitivity_wealth_diff <- lm(wealth_diff_w3_w1 ~ transitivity.borrow + 
#                                            vill_size, data = vill_df)
# lm_borrow_reciprocity_wealth_diff <- lm(wealth_diff_w3_w1 ~ reciprocity.borrow + 
#                                           vill_size, data = vill_df)
# lm_borrow_density_wealth_diff <- lm(wealth_diff_w3_w1 ~ density.borrow + 
#                                       vill_size, data = vill_df)
# 
# # Run linear model to see if presence of kcycle of length 4 are correlated with village health or village
# # mental health at w1
# lm_kcycle4_health_w1 <- lm(avg_health ~ kcycle4 + vill_size, data = vill_df)
# summary(lm_kcycle4_health_w1)
# 
# lm_kcycle4_mental_health_w1 <- lm(avg_mental_health ~ kcycle4 + vill_size, data = vill_df)
# summary(lm_kcycle4_mental_health_w1)
# 
# # Create balloon plot of median hh wealth index per village against kcycle4 levels
# vill_df.balloon <- vill_df %>%
#   select(kcycle2, kcycle3, kcycle4, kcycle5, village_wealth_index_w1, village_wealth_index_w3, 
#          vill_size)
# ggballoonplot(vill_df.balloon)
# 
# # Balloon plot of village_wealth_w1 by kcycle4 prevalence, visualized with village size
# ggplot(vill_df, aes(x = kcycle2, y = wealth_diff_w3_w1, size = vill_size)) +
#   geom_point(shape = 21, colour = "black", fill = "cornsilk") +
#   scale_size_area(max_size = 10) +
#   theme_classic()
# 
# 
# p <- ggplot(vill_df, aes(x = village_code, y = c(kcycle2, kcycle3, kcycle4, kcycle5, village_wealth_index_w1, 
#                                                  village_wealth_index_w3, vill_size))) 
# p+geom_point( aes(size=value))+theme(panel.background=element_blank(), panel.border = element_rect(colour = "blue", fill=NA, size=1))

# # Load cyclic density measure for correlations
# load("~/Desktop/Files/Work/lab/christakis_lab/files/data/borrow/cyc_density_20230911.rda")
# load("~/Desktop/Files/Work/lab/christakis_lab/files/data/borrow/cycdens_borrow_lend_20231010.rda")
# load("~/Desktop/Files/Work/lab/christakis_lab/files/data/borrow/hh_cycdens_20231010.rda")

# # Load path density measure for correlations
# load("~/Desktop/Files/Work/lab/christakis_lab/files/data/borrow/pathdens_ER_20231011.rda")
# load("~/Desktop/Files/Work/lab/christakis_lab/files/data/borrow/pathdens_BA_20231023.rda")
# load("~/Desktop/Files/Work/lab/christakis_lab/files/data/borrow/pathdens_SDS_20231024.rda")
# load("~/Desktop/Files/Work/lab/christakis_lab/files/data/borrow/pathdens_PL_20231023.rda")
# 
# load("~/Desktop/Files/Work/lab/christakis_lab/files/data/borrow/hh_pathdens_ER_20231011.rda")
# load("~/Desktop/Files/Work/lab/christakis_lab/files/data/borrow/hh_pathdens_BA_20231023.rda")
# load("~/Desktop/Files/Work/lab/christakis_lab/files/data/borrow/hh_pathdens_SDS_20231024.rda")
# load("~/Desktop/Files/Work/lab/christakis_lab/files/data/borrow/hh_pathdens_PL_20231023.rda")

# # Kstats for borrow ties only
# kstats_indiv_borrow_df <- vill_df %>%
#   bind_cols(cycdens) %>%
#   rename(cycdens2_er  = kcycle2,
#          cycdens3_er  = kcycle3,
#          cycdens4_er  = kcycle4,
#          cycdens5_er  = kcycle5) %>%
#   bind_cols(cycdens_ba) %>%
#   rename(cycdens2_ba  = kcycle2,
#          cycdens3_ba  = kcycle3,
#          cycdens4_ba  = kcycle4,
#          cycdens5_ba  = kcycle5) %>%
#   bind_cols(cycdens_pl) %>%
#   rename(cycdens2_pl  = kcycle2,
#          cycdens3_pl  = kcycle3,
#          cycdens4_pl  = kcycle4,
#          cycdens5_pl  = kcycle5) %>%
#   bind_cols(cycdens_sds) %>%
#   rename(cycdens2_sds  = kcycle2,
#          cycdens3_sds  = kcycle3,
#          cycdens4_sds  = kcycle4,
#          cycdens5_sds  = kcycle5) %>%
#   bind_cols(pathdens) %>%
#   rename(pathdens2_er = kpath2,
#          pathdens3_er = kpath3,
#          pathdens4_er = kpath4,
#          pathdens5_er = kpath5) %>%
#   bind_cols(pathdens_ba) %>%
#   rename(pathdens2_ba = kpath2,
#          pathdens3_ba = kpath3,
#          pathdens4_ba = kpath4,
#          pathdens5_ba = kpath5) %>%
#   bind_cols(pathdens_pl) %>%
#   rename(pathdens2_pl = kpath2,
#          pathdens3_pl = kpath3,
#          pathdens4_pl = kpath4,
#          pathdens5_pl = kpath5) %>%
#   bind_cols(pathdens_sds) %>%
#   rename(pathdens2_sds = kpath2,
#          pathdens3_sds = kpath3,
#          pathdens4_sds = kpath4,
#          pathdens5_sds = kpath5)
# 
# kstats_hh_borrow_df <- vill_df %>%
#   bind_cols(hh_cycdens) %>%
#   rename(cycdens2_er  = kcycle2,
#          cycdens3_er  = kcycle3,
#          cycdens4_er  = kcycle4,
#          cycdens5_er  = kcycle5) %>%
#   bind_cols(hh_cycdens_ba) %>%
#   rename(cycdens2_ba  = kcycle2,
#          cycdens3_ba  = kcycle3,
#          cycdens4_ba  = kcycle4,
#          cycdens5_ba  = kcycle5) %>%
#   bind_cols(hh_cycdens_pl) %>%
#   rename(cycdens2_pl  = kcycle2,
#          cycdens3_pl  = kcycle3,
#          cycdens4_pl  = kcycle4,
#          cycdens5_pl  = kcycle5) %>%
#   bind_cols(hh_cycdens_sds) %>%
#   rename(cycdens2_sds  = kcycle2,
#          cycdens3_sds  = kcycle3,
#          cycdens4_sds  = kcycle4,
#          cycdens5_sds  = kcycle5) %>%
#   bind_cols(hh_pathdens) %>%
#   rename(pathdens2_er = kpath2,
#          pathdens3_er = kpath3,
#          pathdens4_er = kpath4,
#          pathdens5_er = kpath5) %>%
#   bind_cols(hh_pathdens_ba) %>%
#   rename(pathdens2_ba = kpath2,
#          pathdens3_ba = kpath3,
#          pathdens4_ba = kpath4,
#          pathdens5_ba = kpath5) %>%
#   bind_cols(hh_pathdens_pl) %>%
#   rename(pathdens2_pl = kpath2,
#          pathdens3_pl = kpath3,
#          pathdens4_pl = kpath4,
#          pathdens5_pl = kpath5) %>%
#   bind_cols(hh_pathdens_sds) %>%
#   rename(pathdens2_sds = kpath2,
#          pathdens3_sds = kpath3,
#          pathdens4_sds = kpath4,
#          pathdens5_sds = kpath5)

# # Using erdos-renyi rndom graphs
# lm_cycpath2_borrow_er <- lm(village_wealth_index_w1 ~ cycdens2_er + pathdens2_er, data = kstats_indiv_borrow_df) # <0.001 sig
# lm_cycpath2_hh_borrow_er <- lm(village_wealth_index_w1 ~ cycdens2_er + pathdens2_er, data = kstats_hh_borrow_df) # <0.001 sig
# 
# lm_cycpath3_borrow_er <- lm(village_wealth_index_w1 ~ cycdens3_er + pathdens3_er, data = kstats_indiv_borrow_df) # 0.017 sig
# lm_cycpath3_hh_borrow_er <- lm(village_wealth_index_w1 ~ cycdens3_er + pathdens3_er, data = kstats_hh_borrow_df) # 0.009 sig
# 
# lm_cycpath4_borrow_er <- lm(village_wealth_index_w1 ~ cycdens4_er + pathdens4_er, data = kstats_indiv_borrow_df)
# lm_cycpath4_hh_borrow_er <- lm(village_wealth_index_w1 ~ cycdens4_er + pathdens4_er, data = kstats_hh_borrow_df)
# 
# lm_cycpath5_borrow_er <- lm(village_wealth_index_w1 ~ cycdens5_er + pathdens5_er, data = kstats_indiv_borrow_df)
# lm_cycpath5_hh_borrow_er <- lm(village_wealth_index_w1 ~ cycdens5_er + pathdens5_er, data = kstats_hh_borrow_df)

# # Using barabasi-albert random graphs
# lm_cycpath2_borrow_ba <- lm(village_wealth_index_w1 ~ cycdens2_ba + pathdens2_ba, data = kstats_indiv_borrow_df) # <0.001 sig
# lm_cycpath2_hh_borrow_ba <- lm(village_wealth_index_w1 ~ cycdens2_ba + pathdens2_ba, data = kstats_hh_borrow_df) # <0.001 sig
# 
# lm_cycpath3_borrow_ba <- lm(village_wealth_index_w1 ~ cycdens3_ba + pathdens3_ba, data = kstats_indiv_borrow_df) # 0.017 sig
# lm_cycpath3_hh_borrow_ba <- lm(village_wealth_index_w1 ~ cycdens3_ba + pathdens3_ba, data = kstats_hh_borrow_df) # 0.009 sig
# 
# lm_cycpath4_borrow_ba <- lm(village_wealth_index_w1 ~ cycdens4_ba + pathdens4_ba, data = kstats_indiv_borrow_df)
# lm_cycpath4_hh_borrow_ba <- lm(village_wealth_index_w1 ~ cycdens4_ba + pathdens4_ba, data = kstats_hh_borrow_df)
# 
# lm_cycpath5_borrow_ba <- lm(village_wealth_index_w1 ~ cycdens5_ba + pathdens5_ba, data = kstats_indiv_borrow_df)
# lm_cycpath5_hh_borrow_ba <- lm(village_wealth_index_w1 ~ cycdens5_ba + pathdens5_ba, data = kstats_hh_borrow_df)

# # Using the power law distribution for random graphs
# lm_cycpath2_borrow_pl <- lm(village_wealth_index_w1 ~ cycdens2_pl + pathdens2_pl, data = kstats_indiv_borrow_df) # <0.001 sig
# lm_cycpath2_hh_borrow_pl <- lm(village_wealth_index_w1 ~ cycdens2_pl + pathdens2_pl, data = kstats_hh_borrow_df) # <0.001 sig
# 
# lm_cycpath3_borrow_pl <- lm(village_wealth_index_w1 ~ cycdens3_pl + pathdens3_pl, data = kstats_indiv_borrow_df) # 0.017 sig
# lm_cycpath3_hh_borrow_pl <- lm(village_wealth_index_w1 ~ cycdens3_pl + pathdens3_pl, data = kstats_hh_borrow_df) # 0.009 sig
# 
# lm_cycpath4_borrow_pl <- lm(village_wealth_index_w1 ~ cycdens4_pl + pathdens4_pl, data = kstats_indiv_borrow_df)
# lm_cycpath4_hh_borrow_pl <- lm(village_wealth_index_w1 ~ cycdens4_pl + pathdens4_pl, data = kstats_hh_borrow_df)
# 
# lm_cycpath5_borrow_pl <- lm(village_wealth_index_w1 ~ cycdens5_pl + pathdens5_pl, data = kstats_indiv_borrow_df)
# lm_cycpath5_hh_borrow_pl <- lm(village_wealth_index_w1 ~ cycdens5_pl + pathdens5_pl, data = kstats_hh_borrow_df)

# # Using the same degree distribution as observed for random graphs
# lm_cycpath2_borrow_sds <- lm(village_wealth_index_w1 ~ cycdens2_sds + pathdens2_sds, data = kstats_indiv_borrow_df) # <0.001 sig
# lm_cycpath2_hh_borrow_sds <- lm(village_wealth_index_w1 ~ cycdens2_sds + pathdens2_sds, data = kstats_hh_borrow_df) # <0.001 sig
# 
# lm_cycpath3_borrow_sds <- lm(village_wealth_index_w1 ~ cycdens3_sds + pathdens3_sds, data = kstats_indiv_borrow_df) # 0.017 sig
# lm_cycpath3_hh_borrow_sds <- lm(village_wealth_index_w1 ~ cycdens3_sds + pathdens3_sds, data = kstats_hh_borrow_df) # 0.009 sig
# 
# lm_cycpath4_borrow_sds <- lm(village_wealth_index_w1 ~ cycdens4_sds + pathdens4_sds, data = kstats_indiv_borrow_df)
# lm_cycpath4_hh_borrow_sds <- lm(village_wealth_index_w1 ~ cycdens4_sds + pathdens4_sds, data = kstats_hh_borrow_df)
# 
# lm_cycpath5_borrow_sds <- lm(village_wealth_index_w1 ~ cycdens5_sds + pathdens5_sds, data = kstats_indiv_borrow_df)
# lm_cycpath5_hh_borrow_sds <- lm(village_wealth_index_w1 ~ cycdens5_sds + pathdens5_sds, data = kstats_hh_borrow_df)

#------------------------ Save data ---------------------------
# Saving data to data folder; saved files are ready for figures and tables.

# Save models with MCA wealth as DV
save(lm_cycpath2_d_bl_k_df_mca, lm_cycpath3_d_bl_k_df_mca, lm_cycpath4_d_bl_k_df_mca, lm_cycpath5_d_bl_k_df_mca,
     lm_cycpath2_d_bl_nk_df_mca, lm_cycpath3_d_bl_nk_df_mca, lm_cycpath4_d_bl_nk_df_mca, lm_cycpath5_d_bl_nk_df_mca,
     lm_cycpath2_ud_bl_k_df_mca, lm_cycpath3_ud_bl_k_df_mca, lm_cycpath4_ud_bl_k_df_mca, lm_cycpath5_ud_bl_k_df_mca,
     lm_cycpath2_ud_bl_nk_df_mca, lm_cycpath3_ud_bl_nk_df_mca, lm_cycpath4_ud_bl_nk_df_mca, lm_cycpath5_ud_bl_nk_df_mca,
     lm_cycpath2_d_fr_k_df_mca, lm_cycpath3_d_fr_k_df_mca, lm_cycpath4_d_fr_k_df_mca, lm_cycpath5_d_fr_k_df_mca,
     lm_cycpath2_d_fr_nk_df_mca, lm_cycpath3_d_fr_nk_df_mca, lm_cycpath4_d_fr_nk_df_mca, lm_cycpath5_d_fr_nk_df_mca,
     lm_cycpath2_ud_fr_k_df_mca, lm_cycpath3_ud_fr_k_df_mca, lm_cycpath4_ud_fr_k_df_mca, lm_cycpath5_ud_fr_k_df_mca,
     lm_cycpath2_ud_fr_nk_df_mca, lm_cycpath3_ud_fr_nk_df_mca, lm_cycpath4_ud_fr_nk_df_mca, lm_cycpath5_ud_fr_nk_df_mca,
     file = file.path(data_path, "mca_sds_models_20240613.rda"))

# Save models with vill wealth (avg quintiles) as DV
save(lm_cycpath2_d_bl_k_df_vill_wealth, lm_cycpath3_d_bl_k_df_vill_wealth, lm_cycpath4_d_bl_k_df_vill_wealth, lm_cycpath5_d_bl_k_df_vill_wealth,
     lm_cycpath2_d_bl_nk_df_vill_wealth, lm_cycpath3_d_bl_nk_df_vill_wealth, lm_cycpath4_d_bl_nk_df_vill_wealth, lm_cycpath5_d_bl_nk_df_vill_wealth,
     lm_cycpath2_ud_bl_k_df_vill_wealth, lm_cycpath3_ud_bl_k_df_vill_wealth, lm_cycpath4_ud_bl_k_df_vill_wealth, lm_cycpath5_ud_bl_k_df_vill_wealth,
     lm_cycpath2_ud_bl_nk_df_vill_wealth, lm_cycpath3_ud_bl_nk_df_vill_wealth, lm_cycpath4_ud_bl_nk_df_vill_wealth, lm_cycpath5_ud_bl_nk_df_vill_wealth,
     lm_cycpath2_d_fr_k_df_vill_wealth, lm_cycpath3_d_fr_k_df_vill_wealth, lm_cycpath4_d_fr_k_df_vill_wealth, lm_cycpath5_d_fr_k_df_vill_wealth,
     lm_cycpath2_d_fr_nk_df_vill_wealth, lm_cycpath3_d_fr_nk_df_vill_wealth, lm_cycpath4_d_fr_nk_df_vill_wealth, lm_cycpath5_d_fr_nk_df_vill_wealth,
     lm_cycpath2_ud_fr_k_df_vill_wealth, lm_cycpath3_ud_fr_k_df_vill_wealth, lm_cycpath4_ud_fr_k_df_vill_wealth, lm_cycpath5_ud_fr_k_df_vill_wealth,
     lm_cycpath2_ud_fr_nk_df_vill_wealth, lm_cycpath3_ud_fr_nk_df_vill_wealth, lm_cycpath4_ud_fr_nk_df_vill_wealth, lm_cycpath5_ud_fr_nk_df_vill_wealth,
     file = file.path(data_path, "vill_wealth_sds_models_20240613.rda"))

# Save models with vill gini (MCA coordinates) as DV
save(lm_cycpath2_d_bl_k_df_gini, lm_cycpath3_d_bl_k_df_gini, lm_cycpath4_d_bl_k_df_gini, lm_cycpath5_d_bl_k_df_gini,
     lm_cycpath2_d_bl_nk_df_gini, lm_cycpath3_d_bl_nk_df_gini, lm_cycpath4_d_bl_nk_df_gini, lm_cycpath5_d_bl_nk_df_gini,
     lm_cycpath2_ud_bl_k_df_gini, lm_cycpath3_ud_bl_k_df_gini, lm_cycpath4_ud_bl_k_df_gini, lm_cycpath5_ud_bl_k_df_gini,
     lm_cycpath2_ud_bl_nk_df_gini, lm_cycpath3_ud_bl_nk_df_gini, lm_cycpath4_ud_bl_nk_df_gini, lm_cycpath5_ud_bl_nk_df_gini,
     lm_cycpath2_d_fr_k_df_gini, lm_cycpath3_d_fr_k_df_gini, lm_cycpath4_d_fr_k_df_gini, lm_cycpath5_d_fr_k_df_gini,
     lm_cycpath2_d_fr_nk_df_gini, lm_cycpath3_d_fr_nk_df_gini, lm_cycpath4_d_fr_nk_df_gini, lm_cycpath5_d_fr_nk_df_gini,
     lm_cycpath2_ud_fr_k_df_gini, lm_cycpath3_ud_fr_k_df_gini, lm_cycpath4_ud_fr_k_df_gini, lm_cycpath5_ud_fr_k_df_gini,
     lm_cycpath2_ud_fr_nk_df_gini, lm_cycpath3_ud_fr_nk_df_gini, lm_cycpath4_ud_fr_nk_df_gini, lm_cycpath5_ud_fr_nk_df_gini,
     file = file.path(data_path, "gini_sds_models_20240613.rda"))


# # Save the data frame and models (wave 1)
# save(lm_cycpath2_borrow_lend_er, lm_cycpath2_hh_borrow_lend_er, lm_cycpath3_borrow_lend_er,
#      lm_cycpath3_hh_borrow_lend_er, lm_cycpath4_borrow_lend_er, lm_cycpath4_hh_borrow_lend_er,
#      lm_cycpath5_borrow_lend_er, lm_cycpath5_hh_borrow_lend_er, file = file.path(data_path, "models_w1_er_20231213.rda"))
# 
# save(lm_cycpath2_borrow_lend_ba, lm_cycpath2_hh_borrow_lend_ba, lm_cycpath3_borrow_lend_ba,
#      lm_cycpath3_hh_borrow_lend_ba, lm_cycpath4_borrow_lend_ba, lm_cycpath4_hh_borrow_lend_ba,
#      lm_cycpath5_borrow_lend_ba, lm_cycpath5_hh_borrow_lend_ba, file = file.path(data_path, "models_w1_ba_20231213.rda"))
# 
# save(lm_cycpath2_borrow_lend_pl, lm_cycpath2_hh_borrow_lend_pl, lm_cycpath3_borrow_lend_pl,
#      lm_cycpath3_hh_borrow_lend_pl, lm_cycpath4_borrow_lend_pl, lm_cycpath4_hh_borrow_lend_pl,
#      lm_cycpath5_borrow_lend_pl, lm_cycpath5_hh_borrow_lend_pl, file = file.path(data_path, "models_w1_pl_20231213.rda"))
# 
# save(lm_cycpath2_borrow_lend_sds, lm_cycpath2_hh_borrow_lend_sds, lm_cycpath3_borrow_lend_sds,
#      lm_cycpath3_hh_borrow_lend_sds, lm_cycpath4_borrow_lend_sds, lm_cycpath4_hh_borrow_lend_sds,
#      lm_cycpath5_borrow_lend_sds, lm_cycpath5_hh_borrow_lend_sds, file = file.path(data_path, "models_w1_sds_20231213.rda"))
# 
# # wave 3
# save(lm_cycpath2_borrow_lend_w3_er, lm_cycpath2_hh_borrow_lend_w3_er,
#      lm_cycpath3_borrow_lend_w3_er, lm_cycpath3_hh_borrow_lend_w3_er,
#      lm_cycpath4_borrow_lend_w3_er, lm_cycpath4_hh_borrow_lend_w3_er,
#      lm_cycpath5_borrow_lend_w3_er, lm_cycpath5_hh_borrow_lend_w3_er,
#      file = file.path(data_path, "models_w3_er_20231213.rda"))
# 
# save(lm_cycpath2_borrow_lend_w3_ba, lm_cycpath2_hh_borrow_lend_w3_ba,
#      lm_cycpath3_borrow_lend_w3_ba, lm_cycpath3_hh_borrow_lend_w3_ba,
#      lm_cycpath4_borrow_lend_w3_ba, lm_cycpath4_hh_borrow_lend_w3_ba,
#      lm_cycpath5_borrow_lend_w3_ba, lm_cycpath5_hh_borrow_lend_w3_ba,
#      file = file.path(data_path, "models_w3_ba_20231213.rda"))
# 
# save(lm_cycpath2_borrow_lend_w3_pl, lm_cycpath2_hh_borrow_lend_w3_pl,
#      lm_cycpath3_borrow_lend_w3_pl, lm_cycpath3_hh_borrow_lend_w3_pl,
#      lm_cycpath4_borrow_lend_w3_pl, lm_cycpath4_hh_borrow_lend_w3_pl,
#      lm_cycpath5_borrow_lend_w3_pl, lm_cycpath5_hh_borrow_lend_w3_pl,
#      file = file.path(data_path, "models_w3_pl_20231213.rda"))
# 
# save(lm_cycpath2_hh_borrow_lend_w3_sds, lm_cycpath2_borrow_lend_w3_sds, 
#      lm_cycpath3_hh_borrow_lend_w3_sds, lm_cycpath3_borrow_lend_w3_sds, 
#      lm_cycpath4_hh_borrow_lend_w3_sds, lm_cycpath4_borrow_lend_w3_sds, 
#      lm_cycpath5_hh_borrow_lend_w3_sds, lm_cycpath5_borrow_lend_w3_sds, 
#      file = file.path(data_path, "models_w3_sds_20231213.rda"))
# 
# # combined waves 1 and 3
# save(lm_cycpath2_borrow_lend_w1w3_er, lm_cycpath2_hh_borrow_lend_w1w3_er,
#      lm_cycpath3_borrow_lend_w1w3_er, lm_cycpath3_hh_borrow_lend_w1w3_er,
#      lm_cycpath4_borrow_lend_w1w3_er, lm_cycpath4_hh_borrow_lend_w1w3_er,
#      lm_cycpath5_borrow_lend_w1w3_er, lm_cycpath5_hh_borrow_lend_w1w3_er,
#      file = file.path(data_path, "models_w1w3_er_20231213.rda"))
# 
# save(lm_cycpath2_borrow_lend_w1w3_ba, lm_cycpath2_hh_borrow_lend_w1w3_ba,
#      lm_cycpath3_borrow_lend_w1w3_ba, lm_cycpath3_hh_borrow_lend_w1w3_ba,
#      lm_cycpath4_borrow_lend_w1w3_ba, lm_cycpath4_hh_borrow_lend_w1w3_ba,
#      lm_cycpath5_borrow_lend_w1w3_ba, lm_cycpath5_hh_borrow_lend_w1w3_ba,
#      file = file.path(data_path, "models_w1w3_ba_20231213.rda"))
# 
# save(lm_cycpath2_borrow_lend_w1w3_pl, lm_cycpath2_hh_borrow_lend_w1w3_pl,
#      lm_cycpath3_borrow_lend_w1w3_pl, lm_cycpath3_hh_borrow_lend_w1w3_pl,
#      lm_cycpath4_borrow_lend_w1w3_pl, lm_cycpath4_hh_borrow_lend_w1w3_pl,
#      lm_cycpath5_borrow_lend_w1w3_pl, lm_cycpath5_hh_borrow_lend_w1w3_pl,
#      file = file.path(data_path, "models_w1w3_pl_20231213.rda"))
# 
# save(lm_cycpath2_hh_borrow_lend_w1w3_sds, lm_cycpath2_borrow_lend_w1w3_sds, 
#      lm_cycpath3_hh_borrow_lend_w1w3_sds, lm_cycpath3_borrow_lend_w1w3_sds, 
#      lm_cycpath4_hh_borrow_lend_w1w3_sds, lm_cycpath4_borrow_lend_w1w3_sds, 
#      lm_cycpath5_hh_borrow_lend_w1w3_sds, lm_cycpath5_borrow_lend_w1w3_sds, 
#      file = file.path(data_path, "models_w1w3_sds_20231213.rda"))
