# identify genes/processes enriched by SPP or ENV from AAB annotations

# libraries
library(tidyverse)
library(broom)
library(ggpubr)
library(rstatix)

# data
genome_metadata = read.csv("raw_data/aab_genome_metadata.csv")
dram_annot = read.csv("raw_data/aab_genome_dram_annot.csv", check.names=F)
pan_annot = read.csv("raw_data/all_61_pan_gene_clusters_summary.txt",
                       sep= "\t")
anvio_enrich_61 = read.csv("raw_data/anvio_functional_enrichment_61_env.csv")

# format
pan_names = pan_annot %>%
  select(genome_name) %>%
  distinct()

# add cog cats
anvio_enrich_61_sourdough = anvio_enrich_61 %>%
  filter(associated_groups %in% c("Sourdough"),
         adjusted_q_value <= 0.1)
  
# test enrichment by species and environment with carbon utilization
dram_annot_clean = dram_annot %>%
  filter(Category %in% c("carbon utilization",
                         "MISC",
                         "Organic Nitrogen")) %>%
  select(-gene_description, -module, -header, -subheader,
         -Category, -Average) %>%
  distinct() %>%
  column_to_rownames("gene_id") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("dram_genome_id") %>%
  left_join(genome_metadata) %>%
  select(-genus, -species,
         -dram_genome_id, -pangenome_id) %>%
  relocate(environment_2, genus_spp, environment) %>%
  filter(genus_spp %in% c("Gluconobacter oxydans",
                          "Acetobacter malorum",
                          "Acetobacter orientalis"))

# summarize counts
counts = as.data.frame(colSums(dram_annot_clean[,4:length(dram_annot_clean)])) 
colnames(counts) = "count_obs"
counts = counts %>%
  filter(count_obs > 9) %>%
  rownames_to_column("gene_id") 

# filter to count thresh
dram_annot_filt = dram_annot_clean %>%
  select(environment_2, genus_spp, environment, counts$gene_id)

# filter to major pathways
cutil = dram_annot %>% filter(Category == "carbon utilization") %>%
  filter(gene_id %in% colnames(dram_annot_filt))

misc = dram_annot %>% filter(Category == "MISC") %>%
  filter(gene_id %in% colnames(dram_annot_filt))

on = dram_annot %>% filter(Category == "Organic Nitrogen") %>%
  filter(gene_id %in% colnames(dram_annot_filt))

dram_cutil = dram_annot_filt %>%
  select(environment_2, genus_spp, environment, cutil$gene_id)

dram_misc = dram_annot_filt %>%
  select(environment_2, genus_spp, environment, misc$gene_id)

dram_on = dram_annot_filt %>%
  select(environment_2, genus_spp, environment, on$gene_id)

######### KWs
# C utilization 
kw_env_cutil = dram_cutil %>%
  select(-environment, -genus_spp) %>%
  gather(key, value, -environment_2) %>% 
  group_by(key) %>% 
  do(tidy(kruskal.test(x= .$value, g = .$environment_2))) %>%
  filter(!is.na(p.value))

# correction
p.fdr = p.adjust(kw_env_cutil$p.value, method = "fdr")
kw_env_cutil$p.value.fdr = p.fdr

# add annotations
kw_env_cutil_wAnnot = kw_env_cutil %>%
  left_join(dram_annot, by=c("key" = "gene_id"))

# MISC
kw_env_misc = dram_misc %>%
  select(-environment, -genus_spp) %>%
  gather(key, value, -environment_2) %>% 
  group_by(key) %>% 
  do(tidy(kruskal.test(x= .$value, g = .$environment_2))) %>%
  filter(!is.na(p.value))

# correction
p.fdr = p.adjust(kw_env_misc$p.value, method = "fdr")
kw_env_misc$p.value.fdr = p.fdr

# add annotations
kw_env_misc_wAnnot = kw_env_misc %>%
  left_join(dram_annot, by=c("key" = "gene_id"))

# ON
kw_env_on = dram_on %>%
  select(-environment, -genus_spp) %>%
  gather(key, value, -environment_2) %>% 
  group_by(key) %>% 
  do(tidy(kruskal.test(x= .$value, g = .$environment_2))) %>%
  filter(!is.na(p.value))

# correction
p.fdr = p.adjust(kw_env_on$p.value, method = "fdr")
kw_env_on$p.value.fdr = p.fdr

# add annotations
kw_env_on_wAnnot = kw_env_on %>%
  left_join(dram_annot, by=c("key" = "gene_id"))

# put together and write out
kw_all = rbind(kw_env_on_wAnnot, kw_env_misc_wAnnot, kw_env_cutil_wAnnot) %>%
  dplyr::filter(p.value.fdr <= 0.05)

#write.csv(kw_all, "../AAB_amplicon_add/annotation_stats/enrichment_cats_strain.csv")

####################################################

GH13_df = dram_annot_clean %>%
  select(GH13, genus_spp, environment_2) %>%
  mutate(GH13_pa = ifelse(GH13 > 0, "1", "0")) 

# summary stats 
GH13_df %>%
  group_by(genus_spp, environment_2) %>%
  get_summary_stats(GH13, type = "mean_sd")

# plot - Figure 4C
ggplot(GH13_df, aes(x = fct_reorder(genus_spp, GH13, .fun='mean'), y = GH13)) +
  geom_boxplot(outlier.shape = NA, aes(fill=environment_2), lwd=0.3) +
  geom_jitter(aes(fill=environment_2), pch=21, stroke = 0.25, size = 2.5, width = 0.2) +
  theme_classic() +
  scale_fill_manual(values= c("gray", "#FAF0A8")) 

