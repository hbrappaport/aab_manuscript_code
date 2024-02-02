
################################################################################
# Look at major stats for VOC data - adapted from https://elifesciences.org/articles/61644
################################################################################

# packages
library(mctoolsr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(scales)
library(plyr)
library(viridis)
library(GGally)
library(matrixStats)
library(pheatmap)
library(dendextend)
library(devtools)
library(tidyr)
library(data.table)
library(dplyr)
library(tidyverse)
library(rstatix)

# functions
find_hull <- function(ord_voc_reps) ord_voc_reps[chull(ord_voc_reps[,1], ord_voc_reps[,2]), ]

## VOC data
voc <- read.csv("raw_data/voc_17_highest_avg_nored.csv", row.names = 1)
voc_reps <- read.csv("raw_data/voc_17_highest_reps_all.csv", row.names = 1)
voc_meta <- read.csv("raw_data/voc_17_meta.csv")
ord_voc_reps_t_subset <- read.csv("raw_data/ord_voc_reps_t_subset.csv")

## clean up 
voc[is.na(voc)] <- 0
sort(colSums(voc))
voc <- convert_to_relative_abundances(voc)

## and all reps
voc_reps[is.na(voc_reps)] <- 0
sort(colSums(voc_reps))
voc_reps <- convert_to_relative_abundances(voc_reps)
voc_reps_abu_t <- pivot_longer(voc_reps, X)

## First look at within vs across reps
dm.voc.avg <- calc_dm(voc)
ord_voc_avg <- calc_ordination(dm.voc.avg, ord_type = 'NMDS') 
ord_voc_avg$sampleID <- rownames(ord_voc_avg)
#ord_voc_avg <- cbind(ord_voc_avg, voc_meta$SampleOr)

dm.voc.reps <- calc_dm(voc_reps)
ord_voc_reps <- calc_ordination(dm.voc.reps, ord_type = 'NMDS')
ord_voc_reps$sampleID <- rownames(ord_voc_reps)
ord_voc_reps <- cbind(ord_voc_reps, voc_meta$simpleOr)
ord_voc_reps <- cbind(ord_voc_reps, voc_meta$simpleID)
ord_voc_reps <- cbind(ord_voc_reps, voc_meta$YL_A)
colnames(ord_voc_reps) <- c("MDS1", "MDS2", "sampleID", "simpleOR", "simpleID", "YL_A")

# find and add hulls
hulls <- ddply(ord_voc_reps, "simpleOR", find_hull)


# And permanova for stats # R2 of sample origin explaining reps
adonis2(formula = dm.voc.reps ~ simpleOR, data = ord_voc_reps, permutations = 999)

#permanova for YL vs. all AAB samples
adonis2(formula = dm.voc.reps ~ YL_A, data = ord_voc_reps, permutations = 999)
             
##dm for the VOCs
voc_t = t(voc_reps)
dm.voc.reps.t <- calc_dm(voc_t)
ord_voc_reps_t <- calc_ordination(dm.voc.reps.t, ord_type = 'NMDS')
ord_voc_reps_t$VOCID <- rownames(ord_voc_reps_t)
  
ggplot() +
  geom_point(data = ord_voc_reps, aes(MDS1, MDS2, fill = simpleOR), shape = 21, size = 3) +
  geom_polygon(data = hulls, alpha = 1, aes(MDS1, MDS2, fill = simpleOR)) +
  geom_text(data = ord_voc_reps, aes(MDS1, MDS2, color = simpleOR, label = simpleID), size =3, vjust = 2) +
  scale_fill_manual(values = c("#4363ef", "#fceab4", "#c6c9ea", "#fceab4", "#fceab4", "#4363ef" ,"#d5936c")) +
  scale_color_manual(values = c("#243580", "#998e6d", "#696b7d", "#998e6d", "#998e6d", "#243580" ,"#7d5741")) +
  geom_point(data = ord_voc_reps_t, aes(-MDS1, -MDS2)) +
  geom_text(data = ord_voc_reps_t_subset, aes(-MDS1, -MDS2, label = VOCID), size = 3, vjust = 2) +
  theme_classic() 
    
# calculate summary stats for most abundant VOCs across dataset
voc_means <- as.data.frame(rowMeans(voc)) 
voc_means$SD <- rowSds(as.matrix(voc))
voc_means$Compound <- rownames(voc_means)
voc_means$Mean <- voc_means$`rowMeans(voc)`
voc_means$`rowMeans(voc)` <- NULL
  
###stats####

voc_reps_groups = read.csv("~/Downloads/voc_reps_t.csv")
voc_reps_groups_abu = convert_to_relative_abundances(voc_reps_groups)

voc_reps_groups %>%
  do(tidy(kruskal.test(x= .$Sample, g = .$YL_A)))

pwc <- voc_reps_groups %>% 
  dunn_test(Palmitic.Acid ~ Group, p.adjust.method = "bonferroni") 
pwc


voc_reps_groups_abu = read.csv("~/Downloads/voc_reps_highest_all_abu_t.csv")

voc_reps_groups_filtered = read.csv("~/Downloads/voc_abu_filtered_one_perc.csv")

## look at all, at the same time

voc_yl_a_filt = voc_reps_groups_filtered %>%
  select(-Group, -Sample) %>%
  pivot_longer(!YL_A) %>%
  group_by(name) %>%
  do(tidy(kruskal.test(x= .$value, g = .$YL_A)))

p_adj = stats::p.adjust(p = voc_yl_a_filt$p.value, method = "fdr")


voc_yl_a_filt$p_adj = p_adj



write.csv(voc_yl_a_filt, "~/Downloads/voc_yl_a_stats_filt.csv")

voc_groups = voc_reps_groups_abu %>%
  select(-YL_A, -Sample) %>%
  pivot_longer(!Group) %>%
  group_by(name) %>%
  do(tidy(kruskal.test(x= .$value, g = .$Group))) %>%
  mutate(p_adj = p.adjust(p.value, method = "bonferroni"))


voc_out = voc_reps_groups_abu %>%
  select(-Group, -Sample) %>%
  pivot_longer(!YL_A) %>%
  group_by(name) %>% 
  do(tidy(kruskal.test(x= .$value, g = .$YL_A))) %>%
  mutate(p_adj = p.adjust(pvalue, method = bonferroni))