
# Co-occurrence networks
install.packages("pacman")
pacman::p_load(reshape2, dplyr, mctoolsr, ggplot2, vegan, ggraph, igraph, Hmisc)

#load data
tax_table_its_fp = "raw_data/sourdough_tab_ITS.txt"
tax_table_16s_fp = "raw_data/sourdough_tab_16S.txt"
map_fp = "raw_data/metadata500_starters.txt"
input_ITS = load_taxa_table(tax_table_its_fp, map_fp)
input_16S = load_taxa_table(tax_table_16s_fp, map_fp)

# unique names for euks and bacteria
rownames(input_ITS$data_loaded) <- paste("yeast", rownames(input_ITS$data_loaded), sep = "_")
rownames(input_16S$data_loaded) <- paste("bac", rownames(input_16S$data_loaded), sep = "_")
rownames(input_ITS$taxonomy_loaded) <- paste("yeast", rownames(input_ITS$taxonomy_loaded), sep = "_")
rownames(input_16S$taxonomy_loaded) <- paste("bac", rownames(input_16S$taxonomy_loaded), sep = "_")
full_table <- rbind.data.frame(input_ITS$data_loaded, input_16S$data_loaded)

# taxonomy
bac_tax <- input_16S$taxonomy_loaded
yeast_tax <- input_ITS$taxonomy_loaded
yeast_tax$taxonomy8 <- NULL
bac_tax$taxonomy8 <- NULL
taxonomy <- rbind.data.frame(bac_tax, yeast_tax)
taxonomy$ASVID <- rownames(taxonomy)

# filter by abundance
sort(colSums(full_table))
full_table_filt <- filter_taxa_from_table(full_table, filter_thresh = 0.001)

# and only include if found in more than 10% of sites
full_table_filt_PA <- decostand(full_table_filt, method = "pa") #same number here

# set number of samples threshold (>= 1 site)
site.num <- ncol(full_table_filt_PA)
no_sites <- sort(rowSums(full_table_filt_PA))
no_sites_index <- (no_sites >= 3)

# must be found in > 10% of all sites
no_sites <- as.data.frame(no_sites[no_sites_index])
no_sites$ASV <- rownames(no_sites)

# filter table by ASVs
full_table_filt$ASVID <- rownames(full_table_filt)
full_table_filt <- filter(full_table_filt, ASVID %in% no_sites$ASV)

# check table
rownames(full_table_filt) <- full_table_filt$ASVID
full_table_filt$ASVID <- NULL
sort(colSums(full_table_filt))

# correlate
full_table_filt <- as.matrix(t(full_table_filt))
all <- rcorr(full_table_filt, type="spearman")

# wrangle p values and rho values for ASV and genus
all_pvalues <- all$P
all_rho <- all$r
all_pvalues <- melt(all_pvalues)
all_pvalues$pvalues <- all_pvalues$value
all_pvalues$value <- NULL
all_rho <- melt(all_rho)
all_rho$Var1 <- NULL
all_rho$Var2 <- NULL
all_edges <- cbind(all_pvalues, all_rho)

# filter by strength and significance ## current is 0.5 with p < 0.00001
all_edges <- filter(all_edges, all_edges$value < -0.2 | all_edges$value > 0.2)
all_edges <- filter(all_edges, all_edges$pvalues < 0.001) # filter by pval
all_nodes <- as.vector(unique(all_edges$Var1))

 # filter to only include yeast-bac correlations for now
 yeast_bac_edges_filt <- all_edges
 
####
 yeast_bac_edges_filt$ASVID <- yeast_bac_edges_filt$Var1
 yeast_bac_edges_filt_wTax <- merge(yeast_bac_edges_filt, taxonomy, by  ="ASVID")
 yeast_bac_edges_filt_wTax$ASVID <- yeast_bac_edges_filt_wTax$Var2
 yeast_bac_edges_filt_wTax <- merge(yeast_bac_edges_filt_wTax, taxonomy, by  ="ASVID")
# 
 aab_edges_filt_wTax_RED <- filter(yeast_bac_edges_filt_wTax,taxonomy5.y == "Acetobacteraceae" | taxonomy5.x == "Acetobacteraceae")
 unique(aab_edges_filt_wTax_RED$Var1)

# and clean up nodes
all_nodes.a <- as.data.frame(unique(aab_edges_filt_wTax_RED$Var1))
all_nodes.a$ASVID <- all_nodes.a$`unique(aab_edges_filt_wTax_RED$Var1)`
all_nodes.a$`unique(aab_edges_filt_wTax_RED$Var1)` <- NULL
all_nodes.b <- as.data.frame(unique(aab_edges_filt_wTax_RED$Var2))
all_nodes.b$ASVID <- all_nodes.b$`unique(aab_edges_filt_wTax_RED$Var2)`
all_nodes.b$`unique(aab_edges_filt_wTax_RED$Var2)` <- NULL
all_nodes <- rbind(all_nodes.a, all_nodes.b)

# merge w tax
all_nodes_wTax <- merge(all_nodes, taxonomy, by = "ASVID")

# add abundances to nodes
abu <- rowMeans(full_table)
ASVID <- rownames(full_table)
abu <- cbind(abu, ASVID)
all_nodes_wTax <- merge(all_nodes_wTax, abu, by = "ASVID")

