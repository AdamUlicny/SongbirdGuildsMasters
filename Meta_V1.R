############### Analysing large scale patterns of Guild membership and phylogenetic signal ###############
library(tidyverse)
library(readxl)
library(rstudioapi)
library(vegan)
library(ape)
library(dendextend)
library(RColorBrewer)
library(factoextra)
library(treemap)
library(patchwork)
library(gridExtra)

# Load the data
sp_list_meta <- read.csv("data/sp_list_meta.csv")
method_substrate_meta <- read.csv("data/method_substrate_meta.csv")

# filtering out non-passerines
method_substrate_meta<-method_substrate_meta%>%
  inner_join(sp_list_meta%>%filter(Passeriformes=="PASSERIFORMES"), by="Sp_BirdLife")

# summarising dataset
method_substrate_subset<-method_substrate_meta%>%
  group_by(Sp_BirdLife)%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  summarise(across(where(is.numeric), sum, na.rm = F))%>%
  filter(N_BEH >= 10  & N_SUB >= 10)%>%
  select(-X)

# make list of species in the subset
list_sp <- method_substrate_subset%>%
  pull(Sp_BirdLife)

write.csv(list_sp, "resources/list_sp_meta_v2.csv")

# Load the new phylogeny matching tips  
phylo_meta <- read.nexus("resources/phylo/phylo_meta_v2/output.nex")

# filtering list for drop.tip
filter_list<-sp_list_meta%>%
  filter(!(Passeriformes=="PASSERIFORMES"))%>%
  pull(Sp_BirdLife)


# Species count per continent
counts_per_continent <- method_substrate_meta%>%
  group_by(continent)%>%
  summarise(n_distinct(Sp_BirdLife))

# matrix preparation
matrix_meta_full <- method_substrate_subset%>%
  select(-c(N_BEH,N_SUB))%>%
  remove_rownames%>%
  column_to_rownames(var="Sp_BirdLife")

# Bray-Curtis distance calculation
matrix_meta_full <- matrix_meta_full%>%
  vegdist(method = "bray")%>%
  as.matrix()%>%
  as.dist(matrix_meta_full[order(rownames(matrix_meta_full)),order(colnames(matrix_meta_full))])

dendro_meta_bray <- matrix_meta_full %>% 
  hclust(method = "ward.D2") %>%
  as.dendrogram

############### Phylogeny ##################
phylo_meta <- phylo_meta[[10]]
phylo_meta <- ape::drop.tip(phylo_meta, filter_list)
phylo_meta <- cophenetic(phylo_meta)
phylo_meta <- as.matrix(phylo_meta)
phylo_meta <- as.dist(phylo_meta[order(rownames(phylo_meta)),order(colnames(phylo_meta))])

dendro_meta_phylo <- hclust(phylo_meta) %>%
  as.dendrogram

# Plotting the dendrogram
par(mar=c(5,1,1,12))
plot(dendro_meta_bray, main = "Foraging guilds Meta ", type = "rectangle", horiz = T)

# Phylogeny - Guilds Co-dendrogram
set.seed(12345)
dendlist(dendro_meta_bray, dendro_meta_phylo)%>%
  dendextend::untangle(method="random", R=100)%>%####### Crucial step to produce human readable codendrograms! Use lower R on slower machines.
  dendextend::untangle(method="step2side")%>%############## For troubleshooting, use "%>% entanglement()" to assess entanglement (lower is better)
  tanglegram(common_subtrees_color_lines = TRUE, # Do NOT include "sort=T" argument if using untangle before (sort overrides it)
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="behavior",
             main_right="phylogeny",
             hang=F)
