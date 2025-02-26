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

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the data
sp_list_meta <- read.csv("data/sp_list_meta.csv")
method_substrate_meta <- read.csv("data/method_substrate_meta.csv")

##### Data preparation ###########
# filtering out non-passerines
method_substrate_meta<-method_substrate_meta%>%
  inner_join(sp_list_meta%>%filter(Passeriformes=="PASSERIFORMES"), by="Sp_BirdLife")

# summarising dataset
method_substrate_subset<-method_substrate_meta%>%
  group_by(Sp_BirdLife)%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  summarise(across(where(is.numeric), sum, na.rm = F))%>%
  filter(N_BEH >= 10  & N_SUB >= 10) # current arbitrary criteria for inclusion in the analysis


method_substrate_continents<-method_substrate_meta%>%
  group_by(continent,Sp_BirdLife)%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  summarise(across(where(is.numeric), sum, na.rm = F))%>%
  filter(N_BEH >= 10  & N_SUB >= 10) # current arbitrary criteria for inclusion in the analysis


# make list of species in  global subset
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
  filter(N_BEH >= 10  & N_SUB >= 10)%>%
  summarise(n_distinct(Sp_BirdLife))

# matrix preparation
matrix_meta_full <- method_substrate_subset%>%
  select(-c(N_BEH,N_SUB))%>%
  remove_rownames%>%
  column_to_rownames(var="Sp_BirdLife")

############## Bray-Curtis distance calculation ##################
dist_Bray_Global <- matrix_meta_full%>%
  vegdist(method = "bray")
dist_Bray_Global<-as.matrix(dist_Bray_Global)
dist_Bray_Global<-as.dist(dist_Bray_Global[order(rownames(dist_Bray_Global)),order(colnames(dist_Bray_Global))])

dendro_meta_bray <- dist_Bray_Global %>% 
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

###########  Plotting the guilds dendrogram ###########
par(mar=c(5,1,1,12))
plot(dendro_meta_bray, main = "Foraging guilds Meta ", type = "rectangle", horiz = T)

############# Phylogeny - Global Guilds Co-dendrogram ####################
set.seed(12345)
dendlist(dendro_meta_bray, dendro_meta_phylo)%>%
  dendextend::untangle(method="random", R=50)%>%####### Crucial step to produce human readable codendrograms! Use lower R on slower machines.
  dendextend::untangle(method="step2side")%>%
  tanglegram(common_subtrees_color_lines = TRUE, # Do NOT include "sort=T" argument if using untangle before (sort overrides it)
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="foraging guilds",
             main_right="phylogeny",
             hang=F)%>%
  entanglement()# lower entanglement = better readability
mantel_Global <- mantel(dist_Bray_Global, phylo_meta, method = "spearman", permutations = 999)
print(mantel_Global)



######################### Function for separate continents###########################
#loop that produces list_sp per continent
# Get the continents
unique_continents <- unique(method_substrate_continents$continent)

# Loop through each continent
for (current_continent in unique_continents) {
  continent_species <- method_substrate_continents %>%
    filter(continent == current_continent) %>%
    pull(Sp_BirdLife)
  filtered_species<-setdiff(list_sp, continent_species)
  assign(paste0("filter_not_in_", current_continent), filtered_species)
}

# create subsets per continen
per_continent <- function(method_substrate_meta) {
  continents <- unique(method_substrate_meta$continent)  # Extract unique continents
  distance_list_continents <- list()  
  
  for (cont in continents) {
    # Subset and process the data
    filtered_data <- method_substrate_meta %>%
      filter(continent == cont) %>%
      group_by(Sp_BirdLife) %>%
      mutate_all(~replace(., is.na(.), 0)) %>%
      summarise(across(where(is.numeric), sum, na.rm = FALSE)) %>%
      filter(N_BEH >= 10 & N_SUB >= 10) %>%
      select(-c(N_BEH, N_SUB)) %>%
      remove_rownames() %>%
      column_to_rownames(var = "Sp_BirdLife")
    
    # Compute Bray-Curtis dissimilarity
    dist_matrix <- vegdist(filtered_data, method = "bray")
    
    # Store distance matrix in the list
    distance_list_continents[[cont]] <- dist_matrix
  }
  
  return(distance_list_continents)  # Return list of distance matrices
}

distance_list_continents<-per_continent(method_substrate_meta)

phylo_meta <- read.nexus("resources/phylo/phylo_meta_v2/output.nex")
phylo_meta <- phylo_meta[[10]]

# Asia
phylo_Asia <- ape::drop.tip(phylo_meta, filter_not_in_Asia)
phylo_Asia <- cophenetic(phylo_Asia)
phylo_Asia <- as.matrix(phylo_Asia)
phylo_Asia <- as.dist(phylo_Asia[order(rownames(phylo_Asia)),order(colnames(phylo_Asia))])
dendro_Asia_phylo <- hclust(phylo_Asia) %>%
  as.dendrogram

#Australia
phylo_Australia <- ape::drop.tip(phylo_meta, filter_not_in_Australia)
phylo_Australia <- cophenetic(phylo_Australia)
phylo_Australia <- as.matrix(phylo_Australia)
phylo_Australia <- as.dist(phylo_Australia[order(rownames(phylo_Australia)),order(colnames(phylo_Australia))])
dendro_Australia_phylo <- hclust(phylo_Australia) %>%
  as.dendrogram


# Europe
phylo_Europe <- ape::drop.tip(phylo_meta, filter_not_in_Europe)
phylo_Europe <- cophenetic(phylo_Europe)
phylo_Europe <- as.matrix(phylo_Europe)
phylo_Europe <- as.dist(phylo_Europe[order(rownames(phylo_Europe)),order(colnames(phylo_Europe))])
dendro_Europe_phylo <- hclust(phylo_Europe) %>%
  as.dendrogram

# North America
filter_not_in_North_America <- `filter_not_in_North America`
rm(`filter_not_in_North America`)
phylo_North_America <- ape::drop.tip(phylo_meta, filter_not_in_North_America)
phylo_North_America <- cophenetic(phylo_North_America)
phylo_North_America <- as.matrix(phylo_North_America)
phylo_North_America <- as.dist(phylo_North_America[order(rownames(phylo_North_America)),order(colnames(phylo_North_America))])
dendro_North_America_phylo <- hclust(phylo_North_America) %>%
  as.dendrogram

# order distance matrices alphabetically and create dendrograms
dist_Bray_Asia <- distance_list_continents[["Asia"]]
dist_Bray_Asia<- as.matrix (dist_Bray_Asia)
dist_Bray_Asia<-as.dist(dist_Bray_Asia[order(rownames(dist_Bray_Asia)),order(colnames(dist_Bray_Asia))])
dendro_Asia_bray <- dist_Bray_Asia%>% 
  hclust(method = "ward.D2") %>%
  as.dendrogram

dist_Bray_Australia <- distance_list_continents[[ "Australia"]]
dist_Bray_Australia<- as.matrix (dist_Bray_Australia)
dist_Bray_Australia<-as.dist(dist_Bray_Australia[order(rownames(dist_Bray_Australia)),order(colnames(dist_Bray_Australia))])
dendro_Australia_bray <- dist_Bray_Australia%>% 
  hclust(method = "ward.D2") %>%
  as.dendrogram

dist_Bray_Europe <- distance_list_continents[["Europe"]]
dist_Bray_Europe<- as.matrix (dist_Bray_Europe)
dist_Bray_Europe<-as.dist(dist_Bray_Europe[order(rownames(dist_Bray_Europe)),order(colnames(dist_Bray_Europe))])
dendro_Europe_bray <- dist_Bray_Europe%>% 
  hclust(method = "ward.D2") %>%
  as.dendrogram

dist_Bray_North_America <- distance_list_continents[["North America"]]
dist_Bray_North_America<- as.matrix (dist_Bray_North_America)
dist_Bray_North_America<-as.dist(dist_Bray_North_America[order(rownames(dist_Bray_North_America)),order(colnames(dist_Bray_North_America))])
dendro_North_America_bray <- dist_Bray_North_America%>% 
  hclust(method = "ward.D2") %>%
  as.dendrogram

### Dendrograms per continent
dendlist(dendro_Europe_bray, dendro_Europe_phylo)%>%
  dendextend::untangle(method="random", R=50)%>%####### Crucial step to produce human readable codendrograms! Use lower R on slower machines.
  dendextend::untangle(method="step2side")%>%
  tanglegram(common_subtrees_color_lines = TRUE, # Do NOT include "sort=T" argument if using untangle before (sort overrides it)
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="foraging guilds",
             main_right="phylogeny",
             main="Europe",
             hang=F)%>%
  entanglement() # lower entanglement = better readability
mantel_Europe <- mantel(dist_Bray_Europe, phylo_Europe, method = "spearman", permutations = 999)
print(mantel_Europe)

dendlist(dendro_Asia_bray, dendro_Asia_phylo)%>%
  dendextend::untangle(method="random", R=50)%>%####### Crucial step to produce human readable codendrograms! Use lower R on slower machines.
  dendextend::untangle(method="step2side")%>%
  tanglegram(common_subtrees_color_lines = TRUE, # Do NOT include "sort=T" argument if using untangle before (sort overrides it)
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="foraging guilds",
             main_right="phylogeny",
             main="Asia",
             hang=F)%>%
  entanglement()# lower entanglement = better readability
mantel_Asia <- mantel(dist_Bray_Asia, phylo_Asia, method = "spearman", permutations = 999)
print(mantel_Asia)

dendlist(dendro_North_America_bray, dendro_North_America_phylo)%>%
  dendextend::untangle(method="random", R=50)%>%####### Crucial step to produce human readable codendrograms! Use lower R on slower machines.
  dendextend::untangle(method="step2side")%>%
  tanglegram(common_subtrees_color_lines = TRUE, # Do NOT include "sort=T" argument if using untangle before (sort overrides it)
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="foraging guilds",
             main_right="phylogeny",
             main="North_America",
             hang=F)%>%
  entanglement()# lower entanglement = better readability
mantel_North_America <- mantel(dist_Bray_North_America, phylo_North_America, method = "spearman", permutations = 999)
print(mantel_North_America)

dendlist(dendro_Australia_bray, dendro_Australia_phylo)%>%
  dendextend::untangle(method="random", R=50)%>%####### Crucial step to produce human readable codendrograms! Use lower R on slower machines.
  dendextend::untangle(method="step2side")%>%
  tanglegram(common_subtrees_color_lines = TRUE, # Do NOT include "sort=T" argument if using untangle before
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="foraging guilds",
             main_right="phylogeny",
             main="Australia",
             hang=F)%>%
  entanglement()# lower entanglement = better readability
mantel_Australia <- mantel(dist_Bray_Australia, phylo_Australia, method = "spearman", permutations = 999)
print(mantel_Australia)
