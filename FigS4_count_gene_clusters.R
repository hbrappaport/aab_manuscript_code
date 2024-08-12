# count basic stats about genes in pangenome etc.

library(tidyverse)

# data
strain = read.csv("../AAB_amplicon_add/annotation_stats/strain_diversity_pan_summary.csv")
pan_annot = read.csv("../AAB_amplicon_add/annotation_stats/all_61_pan_gene_clusters_summary.txt",
                     sep= "\t")

## Number of gene clusters in species genomes

### gene clusters with COG hit

#### iterate to SELECT each of 3 species - first malorum
mal = strain %>%
  filter(Species == "malorum")

mal_clusters = mal %>%
  select(gene_cluster_id) %>%
  distinct()

# core
mal_core = mal %>%
  filter(core_accessory== "core") 

mal_core_clusters = mal_core %>%
  select(gene_cluster_id, num_genes_in_gene_cluster) %>%
  distinct()

mean(mal_core_clusters$num_genes_in_gene_cluster)

# acc
mal_acc = mal %>%
  filter(core_accessory== "accessory") 

mal_acc_clusters = mal_acc %>%
  select(gene_cluster_id, num_genes_in_gene_cluster) %>%
  distinct()

mean(mal_acc_clusters$num_genes_in_gene_cluster)

# oxydans
oxy = strain %>%
  filter(Species == "oxydans")

oxy_clusters = oxy %>%
  select(gene_cluster_id) %>%
  distinct()

# core
oxy_core = oxy %>%
  filter(core_accessory== "core") 

oxy_core_clusters = oxy_core %>%
  select(gene_cluster_id, num_genes_in_gene_cluster) %>%
  distinct()

mean(oxy_core_clusters$num_genes_in_gene_cluster)

# acc
oxy_acc = oxy %>%
  filter(core_accessory== "accessory") 

oxy_acc_clusters = oxy_acc %>%
  select(gene_cluster_id, num_genes_in_gene_cluster) %>%
  distinct()

mean(oxy_acc_clusters$num_genes_in_gene_cluster)

# orientalis
ori = strain %>%
  filter(Species == "orientalis")

ori_clusters = ori %>%
  select(gene_cluster_id) %>%
  distinct()

# core
ori_core = ori %>%
  filter(core_accessory== "core") 

ori_core_clusters = ori_core %>%
  select(gene_cluster_id, num_genes_in_gene_cluster) %>%
  distinct()

mean(ori_core_clusters$num_genes_in_gene_cluster)

# acc
ori_acc = ori %>%
  filter(core_accessory== "accessory") 

ori_acc_clusters = ori_acc %>%
  select(gene_cluster_id, num_genes_in_gene_cluster) %>%
  distinct()

mean(ori_acc_clusters$num_genes_in_gene_cluster)


#############################################################################

#############################################################################

# All AAB

aab_clusters = pan_annot %>%
  select(gene_cluster_id) %>%
  distinct()

# core - 95
aab_95_core = pan_annot %>%
  filter(num_genomes_gene_cluster_has_hits >= 58)

aab_95_core_clusters = aab_95_core %>%
  select(gene_cluster_id, num_genes_in_gene_cluster) %>%
  distinct()

mean(aab_95_core$num_genes_in_gene_cluster)


# acc
aab_95_acc = pan_annot %>%
  filter(!gene_cluster_id %in% aab_95_core_clusters$gene_cluster_id) 

aab_95_acc_clusters = aab_95_acc %>%
  select(gene_cluster_id, num_genes_in_gene_cluster) %>%
  distinct()


# All Acetobacter

ace_annot = pan_annot %>%
  filter(!genome_name %in% c("G_275_3_oxydans_prokka",
                             "G_364_1_oxydans_prokka",
                             "GCF_001580625v1_oxydans_prokka",
                             "GCF_001580635v1_oxydans_prokka",
                             "GCF_001580705v1_oxydans_prokka",
                             "GCF_001580815v1_oxydans_prokka",
                             "GCF_001581045v1_oxydans_prokka",
                             "GCF_024158465v1_oxydans_prokka",
                             "GCF_030450005v1_oxydans_prokka")) 

ace_clusters = ace_annot %>%
  select(gene_cluster_id) %>%
  distinct()

# core - 95
ace_95_core = ace_annot %>%
  filter(num_genomes_gene_cluster_has_hits >= 50)

ace_95_core_clusters = ace_95_core %>%
  select(gene_cluster_id, num_genes_in_gene_cluster) %>%
  distinct()

# acc
ace_95_acc = ace_annot %>%
  filter(!gene_cluster_id %in% ace_95_core_clusters$gene_cluster_id) 

ace_95_acc_clusters = ace_95_acc %>%
  select(gene_cluster_id, num_genes_in_gene_cluster) %>%
  distinct()

############################################################
# subsampling to assess # gene clusters by sample size
############################################################

mal = ori

#### 1 malorum - number of clusters?

mal1 = mal %>%
  dplyr::select(genome_name, gene_cluster_id) %>%
  distinct() %>%
  dplyr::count(genome_name) %>%
  select(-genome_name)

mal1$n_clusters = mal1$n
mal1$n = NULL
mal1$n_group = 1

mal_clusters = mal%>%
  select(genome_name) %>%
  distinct()

#### COMBOS ####

# COMBO2

# assess all combos
comb_2 = as.data.frame(t(combn(mal_clusters$genome_name, 2)))

# initial vector and loop
mal2 = vector(mode = "numeric")

for (i in 1:length(comb_2$V1)) {
  
  clusters = mal %>%
    dplyr::select(genome_name, gene_cluster_id) %>%
    dplyr::filter(mal$genome_name %in%
                    c(comb_2$V1[i],
                      comb_2$V2[i])) %>%
    dplyr::select(gene_cluster_id) %>%
    distinct()
  
  mal2[i] = length(clusters$gene_cluster_id)
  }

# clean up and note n
mal2 = mal2 %>%
  as.data.frame() 
mal2$n_clusters = mal2$.
mal2$. = NULL
mal2$n_group = 2

# COMBO 3
comb_3 = as.data.frame(t(combn(mal_clusters$genome_name, 3)))

# initial vector and loop
mal3 = vector(mode = "numeric")

for (i in 1:length(comb_3$V1)) {
  
  clusters = mal %>%
    dplyr::select(genome_name, gene_cluster_id) %>%
    dplyr::filter(mal$genome_name %in%
                    c(comb_3$V1[i],
                      comb_3$V2[i],
                      comb_3$V3[i])) %>%
    dplyr::select(gene_cluster_id) %>%
    distinct()
  
  mal3[i] = length(clusters$gene_cluster_id)
}

# clean up and note n
mal3 = mal3 %>%
  as.data.frame() 
mal3$n_clusters = mal3$.
mal3$. = NULL
mal3$n_group = 3

# COMBO 4
comb_4 = as.data.frame(t(combn(mal_clusters$genome_name, 4)))

# initial vector and loop
mal4 = vector(mode = "numeric")

for (i in 1:length(comb_4$V1)) {
  
  clusters = mal %>%
    dplyr::select(genome_name, gene_cluster_id) %>%
    dplyr::filter(mal$genome_name %in%
                    c(comb_4$V1[i],
                      comb_4$V2[i],
                      comb_4$V3[i],
                      comb_4$V4[i])) %>%
    dplyr::select(gene_cluster_id) %>%
    distinct()
  
  mal4[i] = length(clusters$gene_cluster_id)
}

# clean up and note n
mal4 = mal4 %>%
  as.data.frame() 
mal4$n_clusters = mal4$.
mal4$. = NULL
mal4$n_group = 4

# COMBO 5
comb_5 = as.data.frame(t(combn(mal_clusters$genome_name, 5)))

# initial vector and loop
mal5 = vector(mode = "numeric")

for (i in 1:length(comb_5$V1)) {
  
  clusters = mal %>%
    dplyr::select(genome_name, gene_cluster_id) %>%
    dplyr::filter(mal$genome_name %in%
                    c(comb_5$V1[i],
                      comb_5$V2[i],
                      comb_5$V3[i],
                      comb_5$V4[i],
                      comb_5$V5[i])) %>%
    dplyr::select(gene_cluster_id) %>%
    distinct()
  
  mal5[i] = length(clusters$gene_cluster_id)
}

# clean up and note n
mal5 = mal5 %>%
  as.data.frame() 
mal5$n_clusters = mal5$.
mal5$. = NULL
mal5$n_group = 5

# COMBO 6
comb_6 = as.data.frame(t(combn(mal_clusters$genome_name, 6)))

# initial vector and loop
mal6 = vector(mode = "numeric")

for (i in 1:length(comb_6$V1)) {
  
  clusters = mal %>%
    dplyr::select(genome_name, gene_cluster_id) %>%
    dplyr::filter(mal$genome_name %in%
                    c(comb_6$V1[i],
                      comb_6$V2[i],
                      comb_6$V3[i],
                      comb_6$V4[i],
                      comb_6$V5[i],
                      comb_6$V6[i])) %>%
    dplyr::select(gene_cluster_id) %>%
    distinct()
  
  mal6[i] = length(clusters$gene_cluster_id)
}

# clean up and note n
mal6 = mal6 %>%
  as.data.frame() 
mal6$n_clusters = mal6$.
mal6$. = NULL
mal6$n_group = 6

# COMBO 7
comb_7 = as.data.frame(t(combn(mal_clusters$genome_name, 7)))

# initial vector and loop
mal7 = vector(mode = "numeric")

for (i in 1:length(comb_7$V1)) {
  
  clusters = mal %>%
    dplyr::select(genome_name, gene_cluster_id) %>%
    dplyr::filter(mal$genome_name %in%
                    c(comb_7$V1[i],
                      comb_7$V2[i],
                      comb_7$V3[i],
                      comb_7$V4[i],
                      comb_7$V5[i],
                      comb_7$V6[i],
                      comb_7$V7[i])) %>%
    dplyr::select(gene_cluster_id) %>%
    distinct()
  
  mal7[i] = length(clusters$gene_cluster_id)
}

# clean up and note n
mal7 = mal7 %>%
  as.data.frame() 
mal7$n_clusters = mal7$.
mal7$. = NULL
mal7$n_group = 7

# COMBO 8
comb_8 = as.data.frame(t(combn(mal_clusters$genome_name, 8)))

# initial vector and loop
mal8 = vector(mode = "numeric")

for (i in 1:length(comb_8$V1)) {
  
  clusters = mal %>%
    dplyr::select(genome_name, gene_cluster_id) %>%
    dplyr::filter(mal$genome_name %in%
                    c(comb_8$V1[i],
                      comb_8$V2[i],
                      comb_8$V3[i],
                      comb_8$V4[i],
                      comb_8$V5[i],
                      comb_8$V6[i],
                      comb_8$V7[i],
                      comb_8$V8[i])) %>%
    dplyr::select(gene_cluster_id) %>%
    distinct()
  
  mal8[i] = length(clusters$gene_cluster_id)
}

# clean up and note n
mal8 = mal8 %>%
  as.data.frame() 
mal8$n_clusters = mal8$.
mal8$. = NULL
mal8$n_group = 8

# COMBO 9
comb_9 = as.data.frame(t(combn(mal_clusters$genome_name, 9)))

# initial vector and loop
mal9 = vector(mode = "numeric")

for (i in 1:length(comb_9$V1)) {
  
  clusters = mal %>%
    dplyr::select(genome_name, gene_cluster_id) %>%
    dplyr::filter(mal$genome_name %in%
                    c(comb_9$V1[i],
                      comb_9$V2[i],
                      comb_9$V3[i],
                      comb_9$V4[i],
                      comb_9$V5[i],
                      comb_9$V6[i],
                      comb_9$V7[i],
                      comb_9$V8[i],
                      comb_9$V9[i])) %>%
    dplyr::select(gene_cluster_id) %>%
    distinct()
  
  mal9[i] = length(clusters$gene_cluster_id)
}

# clean up and note n
mal9 = mal9 %>%
  as.data.frame() 
mal9$n_clusters = mal9$.
mal9$. = NULL
mal9$n_group = 9

#### put all together and plot

mal_subsample = rbind(mal1, mal2, mal3, mal4, mal5, mal6, mal7, mal8, mal9)

ggplot(mal_subsample, aes(n_group, n_clusters)) +
  geom_point(pch=21, fill="darkgray", stroke=0.1, size=2) +
  geom_smooth(color="black") +
  theme_classic() + 
  labs(x= "Number of genomes", y="Number of gene clusters",
       title = "A. malorum")
