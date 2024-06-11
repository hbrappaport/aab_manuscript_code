
# packages needed
library(mctoolsr)
library(tidyverse)
library(forcats)

#------------------------------------------------------------
# input, filter, rarefy 16S
#------------------------------------------------------------

# load data
tax_table_its_fp = 'raw_data/sourdough_tab_ITS.txt'
tax_table_16s_fp = 'raw_data/sourdough_tab_16S.txt'
map_fp = 'raw_data/metadata500_starters.txt'

aab_asv_tree_species = read.csv("raw_data/asv_taxonomy_aab_tree.csv")
aab_genome_species = read.csv("raw_data/genome_taxonomy_aab.csv")

input_ITS = load_taxa_table(tax_table_its_fp, map_fp)
input_16S = load_taxa_table(tax_table_16s_fp, map_fp)

#------------------------------------------------------------
# summarize by a few different metrics
#------------------------------------------------------------

# starter overall % AAB (of bacteria)
input_aab = filter_taxa_from_input(input_16S, taxa_to_keep = c("Acetobacterales"))
mean(colSums(input_aab$data_loaded))
sort(colSums(input_aab$data_loaded))

aab_relabu = as.data.frame(sort(colSums(input_aab$data_loaded)))

# AAB detected in 188/500 samples, 38% of starters
# 147 samples at >= 1% cutoff
sum(aab_relabu$`sort(colSums(input_aab$data_loaded))` > 0)
sum(aab_relabu$`sort(colSums(input_aab$data_loaded))` > 0) / 500
sum(aab_relabu$`sort(colSums(input_aab$data_loaded))` > 0.01)
sum(aab_relabu$`sort(colSums(input_aab$data_loaded))` > 0.01) / 500

aab_asvs_tax = input_aab$taxonomy_loaded

# mean when present
input_aab_thresh = filter_samples_by_counts(input_aab, min_seqs = 0.01)
mean(colSums(input_aab_thresh$data_loaded))
sort(colSums(input_aab_thresh$data_loaded))

# look at co-occurrences 
co_occur_tab = input_aab_thresh$data_loaded
co_occur_tab_pa = co_occur_tab
co_occur_tab_pa[co_occur_tab_pa>0] <-1

mean(colSums(co_occur_tab_pa))

co_occur_tab_pa_long = co_occur_tab_pa %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(!sample_id, names_to = "asv_id", values_to = "rel_abun") %>%
  left_join(aab_asv_tree_species) %>%
  left_join(aab_genome_species) %>%
  filter(rel_abun == 1)

ggplot(co_occur_tab_pa_long, aes(asv_id, sample_id)) +
  geom_point()

# ASV level
aab_table = input_aab$data_loaded %>%
  t() %>%
  as.data.frame()

# long format
aab_table_long = aab_table %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(!sample_id, names_to = "asv_id", values_to = "rel_abun") %>%
  left_join(aab_asv_tree_species) %>%
  left_join(aab_genome_species)

## exclude samples without AAB
aab_table_sub = aab_table %>%
  mutate(total_aab = rowSums(.)) %>%
  filter(total_aab > 0.01) %>%
  dplyr::select(-total_aab) %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(!sample_id, names_to = "asv_id", values_to = "rel_abun") %>%
  left_join(aab_asv_tree_species) %>%
  left_join(aab_genome_species)

# plot Fig. 2B here 
ggplot(aab_table_sub, aes(reorder(sample_id, rel_abun), rel_abun)) +
  geom_bar(stat="identity", aes(fill=color), color="black",
           size = 0.15) +
  scale_fill_manual(values=c("#3A62EC", "gold"), na.value = "#A2A9C7") +
  coord_flip() +
  theme_classic() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(y ="Relative abundance of AAB", x= "500 Sourdough starter sample")

# filter to samples where genome was obtained only
sample_ids_genomes = aab_table_long %>%
  filter(genome_obtained == TRUE) %>%
  select(sample_id) %>%
  distinct()

aab_table_genome_subset = aab_table_long %>%
  filter(sample_id %in% sample_ids_genomes$sample_id)

# species comparison plot
count_asvs = aab_table_genome_subset %>%
  select(asv_id) %>%
  distinct()

# unique species - from ASVs
count_asv_species = aab_table_genome_subset %>%
  select(genus, species) %>%
  distinct()

# unique species - from genomes
count_genome_species = aab_table_genome_subset %>%
  select(genome_name) %>%
  distinct()

# asv count
count_asvs = aab_table_genome_subset %>%
  select(asv_id) %>%
  distinct()

# genome count
count_genomes = aab_table_genome_subset %>%
  filter(genome_obtained == TRUE) 

# to plot on tree
count_asv_species = aab_table_sub %>%
  group_by(genus, species, asv_id) %>%
  plyr::summarize(mean_ra = (mean(rel_abun)*100) ) %>%
  filter(mean_ra > 0.01) %>%
  ungroup()


###Plot Fig. 2A
# relative abundance of ASVs (across samples w/any AAB detected)
ggplot(count_asv_species, aes(reorder(species, mean_ra), mean_ra)) +
  geom_bar(stat = "identity", fill="#A2A9C7", color="black",
           size = 0.15) + 
  labs(x="ASV taxon assignment", y = "Mean % relative abundance") +
  coord_flip() +
  theme_classic()


#------------------------------------------------------------
# look at number of LAB / Yeast in samples with >25% AAB in the bacterial community
#------------------------------------------------------------

# only samples with abundant AAB
input_aab_thresh_10 = filter_samples_by_counts(input_aab, min_seqs = 0.25)
aab_dominant_samples = rownames(input_aab_thresh_10$map_loaded)

# LAB in AAB samples
input_lab = filter_taxa_from_input(input_16S, taxa_to_keep = c("Lactobacillales"))
input_lab$map_loaded$sample_id = rownames(input_lab$map_loaded)

input_lab_filt = filter_data(input_lab, "sample_id",
                        keep_vals = rownames(input_aab_thresh_10$map_loaded))

# test enrichment of LAB and yeast in "aab dominant" vs not samples

# add new col to specify
input_lab$map_loaded = input_lab$map_loaded %>%
  mutate(aab_dominant = ifelse(sample_id %in% aab_dominant_samples, "aab_dominant",
                             "not"))

# check result
input_lab$map_loaded$aab_dominant
sumtax_lab = summarize_taxonomy(input_lab, 7, T, T)

# MW enrichment
kw_lab = taxa_summary_by_sample_type(sumtax_lab, 
                                      metadata_map = input_lab$map_loaded,
                                      test_type = 'MW',
                                      type_header = 'aab_dominant',
                                   filter_level = 0.01)

lab_tax = input_lab$taxonomy_loaded %>%
  rownames_to_column("asv_id")

kw_lab = kw_lab %>%
  rownames_to_column("asv_id") %>%
  left_join(lab_tax)

# Yeast in AAB samples
input_yeast = filter_taxa_from_input(input_ITS, taxa_to_keep = c("Saccharomycetales"))
input_yeast$map_loaded$sample_id = rownames(input_yeast$map_loaded)

input_yeast_filt = filter_data(input_yeast, "sample_id",
                             keep_vals = rownames(input_aab_thresh_10$map_loaded))

# test enrichment of yeast and yeast in "aab dominant" vs not samples

# add new col to specify
input_yeast$map_loaded = input_yeast$map_loaded %>%
  mutate(aab_dominant = ifelse(sample_id %in% aab_dominant_samples, "aab_dominant",
                               "not"))

# check result
input_yeast$map_loaded$aab_dominant
sumtax_yeast = summarize_taxonomy(input_yeast, 7, T, T)

# kw enrichment
kw_yeast = taxa_summary_by_sample_type(sumtax_yeast, 
                                     metadata_map = input_yeast$map_loaded,
                                     test_type = 'MW',
                                     type_header = 'aab_dominant',
                                     filter_level = 0.01)

yeast_tax = input_yeast$taxonomy_loaded %>%
  rownames_to_column("asv_id")

kw_yeast = kw_yeast %>%
  rownames_to_column("asv_id") %>%
  left_join(yeast_tax)


#------------------------------------------------------------
write.csv(kw_yeast, "yeast_enriched_AAB_147.csv")
write.csv(kw_lab, "lab_enriched_AAB_147.csv")

write.csv(sumtax_spp_aab, "aab_species_amplicon.csv")
write.csv(aab_table, "aab_asvs_amplicon.csv")

