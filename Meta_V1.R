############### Analysing large scale patterns of Guild membership and phylogenetic signal ###############

############### Load libraries ###############
library(tidyverse)
library(readxl)
library(rstudioapi)
library(vegan)
library(ape)
library(dendextend)
library(RColorBrewer)
library(factoextra)
library(treemap)
library(clootl)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the data
method_substrate_meta <- read.csv("data/method_substrate_meta.csv") # behavioral data

morphology <- read_xlsx("resources/AVONET/AVONET Supplementary dataset 1.xlsx", sheet = "AVONET2_eBird") # AVONET morphology


##### Data preparation ###########
# version without OTHER substrate, recalculated N_BEH and N_SUB
method_substrate_meta <- method_substrate_meta %>%
  select(-(15:17)) %>%
  mutate_all(~replace(., is.na(.), 0))%>%
  rowwise() %>%
  mutate(N_BEH = sum(c_across(4:9)),
         N_SUB = sum(c_across(10:14)))

# filtering out non-passerines
method_substrate_meta<-method_substrate_meta%>%
  filter(Order=="PASSERIFORMES")

# summarising dataset
method_substrate_subset<-method_substrate_meta%>%
  group_by(Sp_eBird)%>%
  summarise(across(where(is.numeric), sum, na.rm = F))%>%
  filter(N_BEH >= 30  & N_SUB >= 30) # current arbitrary criteria for inclusion in the analysis


method_substrate_continents<-method_substrate_meta%>%
  group_by(continent,Sp_eBird)%>%
  summarise(across(where(is.numeric), sum, na.rm = F))%>%
  filter(N_BEH >= 30  & N_SUB >= 30) # current arbitrary criteria for inclusion in the analysis


############################ Phylogeny###################
sp<-method_substrate_subset%>%
  pull(Sp_eBird)

phylo_meta <- extractTree(species=sp, output.type="scientific", taxonomy.year=2021, version="current")


########################## Morphology #################
morphology<-morphology%>%
  rename(Sp_eBird=Species2)

morphology_Global<-morphology%>%
  inner_join(method_substrate_subset, by="Sp_eBird")%>%
  select(Sp_eBird,Beak.Length_Culmen, Beak.Width, 
         Beak.Depth,Tarsus.Length, Wing.Length, `Hand-Wing.Index`, Tail.Length)%>%
  mutate_at(2:8, log)

morphology_Global_matrix <- morphology_Global%>%
  remove_rownames%>%
  column_to_rownames(var="Sp_eBird")

### distance matrix
dist_morpho_Global<-vegdist(morphology_Global_matrix, method="euclidean")
## alphabetical sorting
dist_morpho_Global <- as.matrix(dist_morpho_Global)
dist_morpho_Global <- as.dist(dist_morpho_Global[order(rownames(dist_morpho_Global)),order(colnames(dist_morpho_Global))])

## cluster 
dendro_Global_morpho<-dist_morpho_Global %>%
  hclust(method="ward.D2")%>%
  as.dendrogram()

#################### Species count per continent ##################
counts_per_continent <- method_substrate_meta%>%
  group_by(continent)%>%
  filter(N_BEH >= 30  & N_SUB >= 30)%>%
  summarise(n_distinct(Sp_eBird))

#################### Matrix preparation ##########################
matrix_meta_full <- method_substrate_subset%>%
  select(-c(N_BEH,N_SUB))%>%
  remove_rownames%>%
  column_to_rownames(var="Sp_eBird")

############## Bray-Curtis dissimilarity calculation ##################
dist_Bray_Global <- matrix_meta_full%>%
  vegdist(method = "bray")
dist_Bray_Global<-as.matrix(dist_Bray_Global)
dist_Bray_Global<-as.dist(dist_Bray_Global[order(rownames(dist_Bray_Global)),order(colnames(dist_Bray_Global))])

dendro_meta_bray <- dist_Bray_Global %>% # clustering using ward.D2
  hclust(method = "ward.D2") %>%
  as.dendrogram

########## Gower distance calculation ################
dist_Gower_Global <- matrix_meta_full%>%
  vegdist(method = "altGower")
dist_Gower_Global<-as.matrix(dist_Gower_Global)
dist_Gower_Global<-as.dist(dist_Gower_Global[order(rownames(dist_Gower_Global)),order(colnames(dist_Gower_Global))])

dendro_meta_gower <- dist_Gower_Global %>% # clustering using ward.D2
  hclust(method = "ward.D2") %>%
  as.dendrogram


############### Phylogeny ##################
phylo_meta_mantel <- as.dist(
  as.matrix(cophenetic(phylo_meta))[order(rownames(cophenetic(phylo_meta))), 
                                    order(colnames(cophenetic(phylo_meta)))]
)

# dendro_meta_phylo <- hclust(phylo_meta, method = "ward.D2") %>%as.dendrogram

phylo_meta_ultra <- chronos(phylo_meta)
phylo_meta_ultra$tip.label <- gsub("_", " ", phylo_meta_ultra$tip.label)
dendro_meta_phylo <- as.dendrogram(phylo_meta_ultra)

###########  Plotting the guilds dendrogram ###########
par(mar=c(5,1,1,12))
plot(dendro_meta_bray, main = "Foraging guilds Bray-Curtis", type = "rectangle", horiz = T)
plot(dendro_meta_gower, main = "Foraging guilds Gower", type = "rectangle", horiz = T)


######################### Separate continents ###########################
#loop that produces list_sp per continent
# Get the continents
unique_continents <- unique(method_substrate_continents$continent)

# Loop through each continent
for (current_continent in unique_continents) {
  continent_species <- method_substrate_continents %>%
    filter(continent == current_continent) %>%
    pull(Sp_eBird)
  assign(paste0("sp_", current_continent), continent_species)
}

# create subsets per continen
per_continent <- function(method_substrate_meta) {
  continents <- unique(method_substrate_meta$continent)  # Extract unique continents
  distance_list_continents <- list()  
  
  for (cont in continents) {
    # Subset and process the data
    filtered_data <- method_substrate_meta %>%
      filter(continent == cont) %>%
      group_by(Sp_eBird) %>%
      mutate_all(~replace(., is.na(.), 0)) %>%
      summarise(across(where(is.numeric), sum, na.rm = FALSE)) %>%
      filter(N_BEH >= 30 & N_SUB >= 30) %>%
      select(-c(N_BEH, N_SUB)) %>%
      remove_rownames() %>%
      column_to_rownames(var = "Sp_eBird")
    
    # Compute Bray-Curtis dissimilarity
    dist_matrix <- vegdist(filtered_data, method = "bray")
    
    # Store distance matrix in the list
    distance_list_continents[[cont]] <- dist_matrix
  }
  
  return(distance_list_continents)  # Return list of distance matrices
}

distance_list_continents<-per_continent(method_substrate_meta)

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

dist_Bray_North_America <- distance_list_continents[["North_America"]]
dist_Bray_North_America<- as.matrix (dist_Bray_North_America)
dist_Bray_North_America<-as.dist(dist_Bray_North_America[order(rownames(dist_Bray_North_America)),order(colnames(dist_Bray_North_America))])
dendro_North_America_bray <- dist_Bray_North_America%>% 
  hclust(method = "ward.D2") %>%
  as.dendrogram

############### Phylogeny per continents ##################
# load phylogeny for each continent
phylo_Asia <- extractTree(species=sp_Asia, output.type="scientific", taxonomy.year=2021, version="current")
phylo_Australia <- extractTree(species=sp_Australia, output.type="scientific", taxonomy.year=2021, version="current")
phylo_Europe <- extractTree(species=sp_Europe, output.type="scientific", taxonomy.year=2021, version="current")
phylo_North_America <- extractTree(species=sp_North_America, output.type="scientific", taxonomy.year=2021, version="current")

#matrix for mantel
phylo_Asia_mantel <- as.dist(
  as.matrix(cophenetic(phylo_Asia))[order(rownames(cophenetic(phylo_Asia))), 
                                    order(colnames(cophenetic(phylo_Asia)))])

phylo_Australia_mantel <- as.dist(
  as.matrix(cophenetic(phylo_Australia))[order(rownames(cophenetic(phylo_Australia))), 
                                    order(colnames(cophenetic(phylo_Australia)))])

phylo_Europe_mantel <- as.dist(
  as.matrix(cophenetic(phylo_Europe))[order(rownames(cophenetic(phylo_Europe))), 
                                    order(colnames(cophenetic(phylo_Europe)))])

phylo_North_America_mantel <- as.dist(
  as.matrix(cophenetic(phylo_North_America))[order(rownames(cophenetic(phylo_North_America))), 
                                    order(colnames(cophenetic(phylo_North_America)))])

# dendro for each continent
phylo_Asia_ultra <- chronos(phylo_Asia)
phylo_Asia_ultra$tip.label <- gsub("_", " ", phylo_Asia_ultra$tip.label)
dendro_Asia_phylo <- as.dendrogram(phylo_Asia_ultra)

phylo_Australia_ultra <- chronos(phylo_Australia)
phylo_Australia_ultra$tip.label <- gsub("_", " ", phylo_Australia_ultra$tip.label)
dendro_Australia_phylo <- as.dendrogram(phylo_Australia_ultra)

phylo_Europe_ultra <- chronos(phylo_Europe)
phylo_Europe_ultra$tip.label <- gsub("_", " ", phylo_Europe_ultra$tip.label)
dendro_Europe_phylo <- as.dendrogram(phylo_Europe_ultra)

phylo_North_America_ultra <- chronos(phylo_North_America)
phylo_North_America_ultra$tip.label <- gsub("_", " ", phylo_North_America_ultra$tip.label)
dendro_North_America_phylo <- as.dendrogram(phylo_North_America_ultra)

#################### Morphology per continent ##################


############## Asia
morphology_Asia<-morphology%>%
  filter(Sp_eBird %in% sp_Asia)%>%
  select(Sp_eBird,Beak.Length_Culmen, Beak.Width, 
         Beak.Depth,Tarsus.Length, Wing.Length, `Hand-Wing.Index`, Tail.Length)%>%
  mutate_at(2:8, log)

morphology_Asia_matrix <- morphology_Asia%>%
  remove_rownames%>%
  column_to_rownames(var="Sp_eBird")

# distance matrix, sorted alphabetically
dist_morpho_Asia<-vegdist(morphology_Asia_matrix, method="euclidean")
dist_morpho_Asia <- as.matrix(dist_morpho_Asia)
dist_morpho_Asia <- as.dist(dist_morpho_Asia[order(rownames(dist_morpho_Asia)),order(colnames(dist_morpho_Asia))])

# clustering 
dendro_Asia_morpho<-dist_morpho_Asia %>%
  hclust(method="ward.D2")%>%
  as.dendrogram()

################ Australia
morphology_Australia<-morphology%>%
  filter(Sp_eBird %in% sp_Australia)%>%
  select(Sp_eBird,Beak.Length_Culmen, Beak.Width, 
         Beak.Depth,Tarsus.Length, Wing.Length, `Hand-Wing.Index`, Tail.Length)%>%
  mutate_at(2:8, log)

morphology_Australia_matrix <- morphology_Australia%>%
  remove_rownames%>%
  column_to_rownames(var="Sp_eBird")

### distance matrix
dist_morpho_Australia<-vegdist(morphology_Australia_matrix, method="euclidean")
## alphabetical sorting
dist_morpho_Australia <- as.matrix(dist_morpho_Australia)
dist_morpho_Australia <- as.dist(dist_morpho_Australia[order(rownames(dist_morpho_Australia)),order(colnames(dist_morpho_Australia))])

## cluster 
dendro_Australia_morpho<-dist_morpho_Australia %>%
  hclust(method="ward.D2")%>%
  as.dendrogram()

################ Europe
morphology_Europe<-morphology%>%
  filter(Sp_eBird %in% sp_Europe)%>%
  select(Sp_eBird,Beak.Length_Culmen, Beak.Width, 
         Beak.Depth,Tarsus.Length, Wing.Length, `Hand-Wing.Index`, Tail.Length)%>%
  mutate_at(2:8, log)

morphology_Europe_matrix <- morphology_Europe%>%
  remove_rownames%>%
  column_to_rownames(var="Sp_eBird")

### distance matrix
dist_morpho_Europe<-vegdist(morphology_Europe_matrix, method="euclidean")
## alphabetical sorting
dist_morpho_Europe <- as.matrix(dist_morpho_Europe)
dist_morpho_Europe <- as.dist(dist_morpho_Europe[order(rownames(dist_morpho_Europe)),order(colnames(dist_morpho_Europe))])

## cluster 
dendro_Europe_morpho<-dist_morpho_Europe %>%
  hclust(method="ward.D2")%>%
  as.dendrogram()


######## North America

morphology_North_America<-morphology%>%
  filter(Sp_eBird %in% sp_North_America)%>%
  select(Sp_eBird,Beak.Length_Culmen, Beak.Width, 
         Beak.Depth,Tarsus.Length, Wing.Length, `Hand-Wing.Index`, Tail.Length)%>%
  mutate_at(2:8, log)

morphology_North_America_matrix <- morphology_North_America%>%
  remove_rownames%>%
  column_to_rownames(var="Sp_eBird")

### distance matrix
dist_morpho_North_America<-vegdist(morphology_North_America_matrix, method="euclidean")
## alphabetical sorting
dist_morpho_North_America <- as.matrix(dist_morpho_North_America)
dist_morpho_North_America <- as.dist(dist_morpho_North_America[order(rownames(dist_morpho_North_America)),order(colnames(dist_morpho_North_America))])

## cluster 
dendro_North_America_morpho<-dist_morpho_North_America %>%
  hclust(method="ward.D2")%>%
  as.dendrogram()



############# Global Co-dendrogram ####################
set.seed(12345)
dendlist(dendro_meta_phylo, dendro_meta_bray)%>%
  dendextend::untangle(method="ladderize")%>%
  tanglegram(common_subtrees_color_lines = TRUE, # Do NOT include "sort=T" argument if using untangle before (sort overrides it)
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="Phylogeny",
             main_right="Bray-Curtis",
             hang=F)%>%
  entanglement()# lower entanglement = better readability

mantel_Global <- mantel(dist_Bray_Global, phylo_meta_mantel, method = "spearman", permutations = 999)
print(mantel_Global)

dendlist(dendro_meta_phylo, dendro_meta_gower)%>%
  dendextend::untangle(method="ladderize")%>%
  tanglegram(common_subtrees_color_lines = TRUE, 
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="Phylogeny",
             main_right="Gower",
             hang=F)%>%
  entanglement()# lower entanglement = better read

mantel_Global2 <- mantel(dist_Gower_Global, phylo_meta_mantel, method = "spearman", permutations = 999)
print(mantel_Global2)

############# Global phylogeny-morphology codendrogram ########################
dendlist(dendro_meta_phylo, dendro_Global_morpho)%>%
  dendextend::untangle(method="ladderize")%>%
  tanglegram(common_subtrees_color_lines = TRUE, 
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="Phylogeny",
             main_right="Morphology",
             hang=F)%>%
  entanglement()# lower entanglement = better read

mantel_Global3 <- mantel(dist_morpho_Global, phylo_meta_mantel, method = "spearman", permutations = 999)
print(mantel_Global3)

############# Global guilds-morphology codendrogram (Bray-Curtis) ########################
dendlist(dendro_meta_bray, dendro_Global_morpho)%>%
  dendextend::untangle(method="ladderize")%>%
  tanglegram(common_subtrees_color_lines = TRUE, 
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="Bray-Curtis guilds",
             main_right="Morphology",
             hang=F)%>%
  entanglement()# lower entanglement = better read

mantel_Global4 <- mantel(dist_morpho_Global, dist_Bray_Global, method = "spearman", permutations = 999)
print(mantel_Global4)



######################## Per-Continent Codendrograms ################################
dendlist(dendro_Asia_phylo, dendro_Asia_bray)%>%
  dendextend::untangle(method="ladderize")%>%
  tanglegram(common_subtrees_color_lines = TRUE,
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="phylogeny",
             main_right="Bray-Curtis",
             main="Asia",
             hang=F)%>%
  entanglement() # lower entanglement = better readability
mantel_Asia <- mantel(dist_Bray_Asia, phylo_Asia_mantel, method = "spearman", permutations = 999)
print(mantel_Asia)

dendlist(dendro_Australia_phylo, dendro_Australia_bray)%>%
  dendextend::untangle(method="ladderize")%>%
  tanglegram(common_subtrees_color_lines = TRUE, # Do NOT include "sort=T" argument if using untangle before (sort overrides it)
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="phylogeny",
             main_right="Bray-Curtis",
             main="Australia",
             hang=F)%>%
  entanglement() # lower entanglement = better readability
mantel_Australia <- mantel(dist_Bray_Australia, phylo_Australia_mantel, method = "spearman", permutations = 999)
print(mantel_Australia)

dendlist(dendro_Europe_phylo, dendro_Europe_bray)%>%
  dendextend::untangle(method="ladderize")%>%
  tanglegram(common_subtrees_color_lines = TRUE, 
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="phylogeny",
             main_right="Bray-Curtis",
             main="Europe",
             hang=F)%>%
  entanglement() # lower entanglement = better readability
mantel_Europe <- mantel(dist_Bray_Europe, phylo_Europe_mantel, method = "spearman", permutations = 999)
print(mantel_Europe)

dendlist(dendro_North_America_phylo, dendro_North_America_bray)%>%
  dendextend::untangle(method="ladderize")%>%
  tanglegram(common_subtrees_color_lines = TRUE, 
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="phylogeny",
             main_right="Bray-Curtis",
             main="North_America",
             hang=F)%>%
  entanglement() # lower entanglement = better readability
mantel_North_America <- mantel(dist_Bray_North_America, phylo_North_America_mantel, method = "spearman", permutations = 999)
print(mantel_North_America)
