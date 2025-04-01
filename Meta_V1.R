############### Analysing large scale patterns of Guild membership and phylogenetic signal ###############
# Test test
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
library(ggraph)
library(igraph)
library(phylobase)
library(phylosignal)
library(adephylo)

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

phylo_meta <- clootl::extractTree(species=sp, output.type="scientific", taxonomy.year=2021, version="current")


########################## Morphology #################
morphology<-morphology%>%
  rename(Sp_eBird=Species2)

morphology_Global<-morphology%>%
  inner_join(method_substrate_subset, by="Sp_eBird")%>%
  select(Sp_eBird,Beak.Length_Culmen, Beak.Width, 
         Beak.Depth,Tarsus.Length, Wing.Length, `Hand-Wing.Index`, Tail.Length)%>%
  mutate_at(2:8, log10)

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

behavior_columns <- colnames(matrix_meta_full)[1:6]
substrate_columns <- colnames(matrix_meta_full)[7:11]

matrix_Global_prop <- matrix_meta_full %>%
  mutate(
    behavior_sum = rowSums(across(all_of(behavior_columns))),
    substrate_sum = rowSums(across(all_of(substrate_columns)))
  ) %>%
  mutate(across(all_of(behavior_columns), ~ ifelse(behavior_sum == 0, 0, . / behavior_sum))) %>%
  mutate(across(all_of(substrate_columns), ~ ifelse(substrate_sum == 0, 0, . / substrate_sum))) %>%
  select(-behavior_sum, -substrate_sum)


matrix_Asia_prop <- method_substrate_continents%>%
  filter(continent=="Asia")%>%
  select(-c(N_BEH,N_SUB))%>%
  remove_rownames%>%
  column_to_rownames(var="Sp_eBird")%>%
  mutate(
    behavior_sum = rowSums(across(all_of(behavior_columns))),
    substrate_sum = rowSums(across(all_of(substrate_columns)))
  ) %>%
  mutate(across(all_of(behavior_columns), ~ ifelse(behavior_sum == 0, 0, . / behavior_sum))) %>%
  mutate(across(all_of(substrate_columns), ~ ifelse(substrate_sum == 0, 0, . / substrate_sum))) %>%
  select(-behavior_sum, -substrate_sum,-continent)

matrix_Australia_prop <-method_substrate_continents%>%
  filter(continent=="Australia")%>%
  select(-c(N_BEH,N_SUB))%>%
  remove_rownames%>%
  column_to_rownames(var="Sp_eBird")%>%
  mutate(
    behavior_sum = rowSums(across(all_of(behavior_columns))),
    substrate_sum = rowSums(across(all_of(substrate_columns)))
  ) %>%
  mutate(across(all_of(behavior_columns), ~ ifelse(behavior_sum == 0, 0, . / behavior_sum))) %>%
  mutate(across(all_of(substrate_columns), ~ ifelse(substrate_sum == 0, 0, . / substrate_sum))) %>%
  select(-behavior_sum, -substrate_sum,-continent)
  
matrix_Europe_prop <-method_substrate_continents%>%
  filter(continent=="Europe")%>%
  select(-c(N_BEH,N_SUB))%>%
  remove_rownames%>%
  column_to_rownames(var="Sp_eBird")%>%
  mutate(
    behavior_sum = rowSums(across(all_of(behavior_columns))),
    substrate_sum = rowSums(across(all_of(substrate_columns)))
  ) %>%
  mutate(across(all_of(behavior_columns), ~ ifelse(behavior_sum == 0, 0, . / behavior_sum))) %>%
  mutate(across(all_of(substrate_columns), ~ ifelse(substrate_sum == 0, 0, . / substrate_sum))) %>%
  select(-behavior_sum, -substrate_sum,-continent)

matrix_North_America_prop <-method_substrate_continents%>%
  filter(continent=="North_America")%>%
  select(-c(N_BEH,N_SUB))%>%
  remove_rownames%>%
  column_to_rownames(var="Sp_eBird")%>%
  mutate(
    behavior_sum = rowSums(across(all_of(behavior_columns))),
    substrate_sum = rowSums(across(all_of(substrate_columns)))
  ) %>%
  mutate(across(all_of(behavior_columns), ~ ifelse(behavior_sum == 0, 0, . / behavior_sum))) %>%
  mutate(across(all_of(substrate_columns), ~ ifelse(substrate_sum == 0, 0, . / substrate_sum))) %>%
  select(-behavior_sum, -substrate_sum,-continent)


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
  vegdist(method = "gower")
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
phylo_Asia <- clootl::extractTree(species=sp_Asia, output.type="scientific", taxonomy.year=2021, version="current")
phylo_Australia <- clootl::extractTree(species=sp_Australia, output.type="scientific", taxonomy.year=2021, version="current")
phylo_Europe <- clootl::extractTree(species=sp_Europe, output.type="scientific", taxonomy.year=2021, version="current")
phylo_North_America <- clootl::extractTree(species=sp_North_America, output.type="scientific", taxonomy.year=2021, version="current")

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
  mutate_at(2:8, log10)

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
  mutate_at(2:8, log10)

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
  mutate_at(2:8, log10)

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
  mutate_at(2:8, log10)

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

############### XXXXXXXXXXXXXXXXX ##########################
############## Graphs ######################################
############# XXXXXXXXXXXXXXXXXX ##########################
########### Coloring Guilds ##############################
passeriformes_meta <- method_substrate_meta %>%
  filter(Order == "PASSERIFORMES") %>%
  select(Sp_eBird, continent) %>%
  group_by(Sp_eBird) %>%
  summarize(continent = ifelse(n_distinct(continent) > 1, "Multiple", first(continent)))  # Mark multiple-continent species
species_to_continent <- setNames(passeriformes_meta$continent, passeriformes_meta$Sp_eBird)
continent_colors<-c("North_America" = "darkblue",
                    "Europe" = "darkgreen",
                    "Asia" = "darkred",
                    "Australia" = "orange",
                    "Multiple" = "black")

# Vector colorCodes where each species from species_to_continent is assigned a color from continent_colors
colorCodes <- sapply(species_to_continent, function(x) continent_colors[x])
# remove .Asia etc
names(colorCodes) <- gsub("\\..*", "", names(colorCodes))



############# Global Co-dendrogram ####################
set.seed(12345)


dendro_meta_phylo<-ladderize(dendro_meta_phylo)# ladderized phylogeny (ape)
dendro_meta_bray<-untangle(dendro_meta_bray, dendro_meta_phylo, method = "step1side")# first dendrogram is untangled, second dendrogram is fixed. Produces
dendro_meta_bray<-dendro_meta_bray[[1]]# saving untangled dendrogram back to dend object from dendlist
# these 2 steps allows us to plot ladderized phylogeny on the right and untangled dendrogram on the left
# without this workflow, we could still do ladderization+step1side, but the phylogeny would be on the left


# extract labels/leaves from dendrogram
labels_phylo_Global<-order.dendrogram(dendro_meta_phylo)
labels_phylo_Global<-colorCodes[labels_phylo_Global]# assign colors to labels

labels_bray_Global<-dendro_meta_bray%>%labels%>%rev()
labels_bray_Global<-colorCodes[labels_bray_Global]# assign colors to labels

labels_colors(dendro_meta_phylo) <- colorCodes[labels(dendro_meta_phylo)][order.dendrogram(dendro_meta_phylo)]
labels_colors(dendro_meta_bray) <- colorCodes[labels(dendro_meta_bray)][order.dendrogram(dendro_meta_bray)]

ordered_labels_phylo <- labels(dendro_meta_phylo)[order.dendrogram(dendro_meta_phylo)]
labels_colors_phylo<-colorCodes[ordered_labels_phylo]
data.frame(Species = ordered_labels_phylo, Color = labels_colors_phylo)

labels_colors(dendro_meta_phylo)<-labels_phylo_Global
labels_colors(dendro_meta_bray)<-labels_bray_Global
colorCodes[groupCodes][order.dendrogram(dend)]



plot(dendro_meta_phylo, horiz=T)

tanglegram(dendro_meta_phylo,dendro_meta_bray,
           common_subtrees_color_lines = TRUE,
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=12.5,
             margin_outer=3,
             lwd=3,
           columns_width=c(5,2,5),
             main_left="Phylogeny",
             main_right="Bray-Curtis",
             hang=F)%>%
  entanglement()# lower entanglement = better readability


mantel_Global_P_G1 <- mantel(dist_Bray_Global, phylo_meta_mantel, method = "spearman", permutations = 999)
print("Results of Mantel test between Bray-Curtis Guilds and Phylogeny")
print(mantel_Global_P_G1)



############# Global phylogeny-morphology codendrogram ########################
dendro_Global_morpho_step1side<-untangle(dendro_Global_morpho, dendro_meta_phylo, method = "step1side")
# second dendrogram is fixed
# this steps allows us to plot ladderized phylogeny and untangled left dendrogram 

tanglegram(dendro_meta_phylo, dendro_Global_morpho_step1side[[1]],
             common_subtrees_color_lines = TRUE, 
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="Phylogeny",
             main_right="Morphology",
             hang=F)%>%
  entanglement()# lower entanglement = better readability

mantel_Global_P_M <- mantel(dist_morpho_Global, phylo_meta_mantel, method = "spearman", permutations = 999)
print("Results of Mantel test between Morphology and Phylogeny")
print(mantel_Global_P_M)

############# Global guilds-morphology codendrogram (Bray-Curtis) ########################
dendlist(dendro_meta_bray, dendro_Global_morpho)%>%
  dendextend::untangle(method="step2side")%>%# untangles both sides for best readability
  tanglegram(common_subtrees_color_lines = TRUE, 
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=,
             lwd=3,
             main_left="Bray-Curtis guilds",
             main_right="Morphology",
             hang=F)%>%
  entanglement()# lower entanglement = better readability

mantel_Global_M_G1 <- mantel(dist_morpho_Global, dist_Bray_Global, method = "spearman", permutations = 999)
print("Results of Mantel test between Bray-Curtis Guilds and Morphology")
print(mantel_Global_M_G1)


mantel_Global_M_G2 <- mantel(dist_morpho_Global, dist_Gower_Global, method = "spearman", permutations = 999)
print("Results of Mantel test between Gower Guilds and Morphology")
print(mantel_Global_M_G2)


######################## Per-Continent Codendrograms ################################
# prepare the dendrograms
dendro_Asia_phylo<-ladderize(dendro_Asia_phylo)# ladderized phylogeny (ape)
dendro_Asia_bray_step1side<-untangle(dendro_Asia_bray, dendro_Asia_phylo, method = "step1side")# untangle bray guilds
dendro_Asia_morpho_step1side<-untangle(dendro_Asia_morpho, dendro_Asia_phylo, method = "step1side")# untangle morphology

dendro_Australia_phylo<-ladderize(dendro_Australia_phylo)# ladderized phylogeny (ape)
dendro_Australia_bray_step1side<-untangle(dendro_Australia_bray, dendro_Australia_phylo, method = "step1side")# untangle bray guilds
dendro_Australia_morpho_step1side<-untangle(dendro_Australia_morpho, dendro_Australia_phylo, method = "step1side")# untangle morphology

dendro_Europe_phylo<-ladderize(dendro_Europe_phylo)# ladderized phylogeny (ape)
dendro_Europe_bray_step1side<-untangle(dendro_Europe_bray, dendro_Europe_phylo, method = "step1side")# untangle bray guilds
dendro_Europe_morpho_step1side<-untangle(dendro_Europe_morpho, dendro_Europe_phylo, method = "step1side")# untangle

dendro_North_America_phylo<-ladderize(dendro_North_America_phylo)# ladderized phylogeny (ape)
dendro_North_America_bray_step1side<-untangle(dendro_North_America_bray, dendro_North_America_phylo, method = "step1side")# untangle bray guilds
dendro_North_America_morpho_step1side<-untangle(dendro_North_America_morpho, dendro_North_America_phylo, method = "step1side")# untangle morphology

#### Asia
tanglegram(dendro_Asia_phylo, dendro_Asia_bray_step1side[[1]],
           common_subtrees_color_lines = TRUE,
           center=T,
           highlight_distinct_edges = FALSE,
           highlight_branches_lwd=FALSE,
           margin_inner=12.5,
           margin_outer=3,
           columns_width=c(5,1,5),
           lwd=2.5,
           main_left="Phylogeny",
           main_right="Bray-Curtis",
           main="Asia",
           hang=F)%>%
  entanglement()
mantel_Asia_P_G1 <- mantel(dist_Bray_Asia, phylo_Asia_mantel, method = "spearman", permutations = 999)
print(mantel_Asia_P_G1)
mantel_Asia_P_M <- mantel(dist_morpho_Asia, phylo_Asia_mantel, method = "spearman", permutations = 999)
mantel_Asia_M_G1 <- mantel(dist_morpho_Asia, dist_Bray_Asia, method = "spearman", permutations = 999)

#### Australia
tanglegram(dendro_Australia_phylo, dendro_Australia_bray_step1side[[1]],
           common_subtrees_color_lines = TRUE,
           center=T,
           highlight_distinct_edges = FALSE,
           highlight_branches_lwd=FALSE,
           margin_inner=12.5,
           margin_outer=5,
           columns_width=c(5,1,5),
           lwd=2.5,
           main_left="Phylogeny",
           main_right="Bray-Curtis",
           main="Asia",
           hang=F)%>%
  entanglement() # lower entanglement = better readability
mantel_Australia_P_G1 <- mantel(dist_Bray_Australia, phylo_Australia_mantel, method = "spearman", permutations = 999)
print(mantel_Australia_P_G1)
mantel_Australia_P_M <- mantel(dist_morpho_Australia, phylo_Australia_mantel, method = "spearman", permutations = 999)
mantel_Australia_M_G1 <- mantel(dist_morpho_Australia, dist_Bray_Australia, method = "spearman", permutations = 999)

#### Europe
tanglegram(dendro_Europe_phylo, dendro_Europe_bray_step1side[[1]],
           common_subtrees_color_lines = TRUE,
           center=T,
           highlight_distinct_edges = FALSE,
           highlight_branches_lwd=FALSE,
           margin_inner=12.5,
           margin_outer=5,
           columns_width=c(5,1,5),
           lwd=2.5,
           main_left="Phylogeny",
           main_right="Bray-Curtis",
           main="Asia",
           hang=F)%>%
  entanglement()# lower entanglement = better readability
mantel_Europe_P_G1 <- mantel(dist_Bray_Europe, phylo_Europe_mantel, method = "spearman", permutations = 999)
print(mantel_Europe_P_G1)
mantel_Europe_P_M <- mantel(dist_morpho_Europe, phylo_Europe_mantel, method = "spearman", permutations = 999)
mantel_Europe_M_G1 <- mantel(dist_morpho_Europe, dist_Bray_Europe, method = "spearman", permutations = 999)

#### North America
tanglegram(dendro_North_America_phylo, dendro_North_America_bray_step1side[[1]],
           common_subtrees_color_lines = TRUE,
           center=T,
           highlight_distinct_edges = FALSE,
           highlight_branches_lwd=FALSE,
           margin_inner=12.5,
           margin_outer=5,
           columns_width=c(5,1,5),
           lwd=2.5,
           main_left="Phylogeny",
           main_right="Bray-Curtis",
           main="Asia",
           hang=F)%>%
  entanglement() # lower entanglement = better readability
mantel_North_America_P_G1 <- mantel(dist_Bray_North_America, phylo_North_America_mantel, method = "spearman", permutations = 999)
print(mantel_North_America_P_G1)
mantel_North_America_P_M <- mantel(dist_morpho_North_America, phylo_North_America_mantel, method = "spearman", permutations = 999)
mantel_North_America_M_G1 <- mantel(dist_morpho_North_America, dist_Bray_North_America, method = "spearman", permutations = 999)

############### Network graphs #########################

###################### Global Mantel Graph ####################################
nodes <- data.frame(
  id = c("Phylogeny", "Morphology", "Guilds Bray"),
  x = c(0, -0.5, 0.5),  # X-coordinates (adjust to keep symmetry)
  y = c(0.5, -0.5, -0.5)  # Y-coordinates (equilateral triangle)
)

nodes2 <- data.frame(
  id = c("Phylogeny", "Morphology", "Guilds Gower"),
  x = c(0, -0.5, 0.5),  # X-coordinates (adjust to keep symmetry)
  y = c(0.5, -0.5, -0.5)  # Y-coordinates (equilateral triangle)
)

edges <- data.frame(
  from = c("Phylogeny", "Morphology", "Guilds Bray"),
  to = c("Guilds Bray", "Phylogeny", "Morphology"),
  weight = c(mantel_Global_P_G1[["statistic"]], mantel_Global_P_M[["statistic"]], mantel_Global_M_G1[["statistic"]])
)

triangle_graph_Global1 <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

ggraph(triangle_graph_Global1, layout = "manual", x = nodes$x, y = nodes$y) + 
  geom_edge_link(aes(width = weight, label = round(weight, 3)), 
                 color = "darkkhaki", edge_alpha = 0.8, 
                 label_size = 5, show.legend = FALSE) +
  geom_node_point(size = 40, color = "cyan4") +  # Bigger nodes
  geom_node_text(aes(label = name), size = 5, color = "black") +  # Labels inside nodes
  xlim(-1, 1) + ylim(-1, 1) + # Adjust edge width scaling
  theme_void() +  # Remove background
  ggtitle("Network of Correlations between Global Guilds, Phylogeny and Morphology")




##################### Asia Mantel Graph #############################
edges_Asia <- data.frame(
  from = c("Phylogeny", "Morphology", "Guilds Bray"),
  to = c("Guilds Bray", "Phylogeny", "Morphology"),
  weight = c(mantel_Asia_P_G1[["statistic"]], mantel_Asia_P_M[["statistic"]], mantel_Asia_M_G1[["statistic"]])
)

triangle_graph_Asia <- graph_from_data_frame(edges_Asia, vertices = nodes, directed = FALSE)

ggraph(triangle_graph_Asia, layout = "manual", x = nodes$x, y = nodes$y) + 
  geom_edge_link(aes(width = weight, label = round(weight, 3)), 
                 color = "darkkhaki", edge_alpha = 0.8, 
                 label_size = 5, show.legend = FALSE) +
  geom_node_point(size = 40, color = "cyan4") +  # Bigger nodes
  geom_node_text(aes(label = name), size = 5, color = "black") +  # Labels inside nodes
  xlim(-1, 1) + ylim(-1, 1) +
  theme_void() +  # Remove background
  ggtitle("Asia")

######################## Australia Mantel Graph #############################

edges_Australia <- data.frame(
  from = c("Phylogeny", "Morphology", "Guilds Bray"),
  to = c("Guilds Bray", "Phylogeny", "Morphology"),
  weight = c(mantel_Australia_P_G1[["statistic"]], mantel_Australia_P_M[["statistic"]], mantel_Australia_M_G1[["statistic"]])
)

triangle_graph_Australia <- graph_from_data_frame(edges_Australia, vertices = nodes, directed = FALSE)

ggraph(triangle_graph_Australia, layout = "manual", x = nodes$x, y = nodes$y) + 
  geom_edge_link(aes(width = weight, label = round(weight, 3)), 
                 color = "darkkhaki", edge_alpha = 0.8, 
                 label_size = 5, show.legend = FALSE) +
  geom_node_point(size = 40, color = "cyan4") +  # Bigger nodes
  geom_node_text(aes(label = name), size = 5, color = "black") +  # Labels inside nodes
  scale_edge_width(range = c(1, 10)) +
  xlim(-1, 1) + ylim(-1, 1) + # Adjust edge width scaling
  theme_void() +  # Remove background
  ggtitle("Australia")

#################### Europe Mantel Graph #############################

edges_Europe <- data.frame(
  from = c("Phylogeny", "Morphology", "Guilds Bray"),
  to = c("Guilds Bray", "Phylogeny", "Morphology"),
  weight = c(mantel_Europe_P_G1[["statistic"]], mantel_Europe_P_M[["statistic"]], mantel_Europe_M_G1[["statistic"]])
)

triangle_graph_Europe <- graph_from_data_frame(edges_Europe, vertices = nodes, directed = FALSE)

ggraph(triangle_graph_Europe, layout = "manual", x = nodes$x, y = nodes$y) + 
  geom_edge_link(aes(width = weight, label = round(weight, 3)), 
                 color = "darkkhaki", edge_alpha = 0.8, 
                 label_size = 5, show.legend = FALSE) +
  geom_node_point(size = 40, color = "cyan4") +  # Bigger nodes
  geom_node_text(aes(label = name), size = 5, color = "black") +  # Labels inside nodes
  scale_edge_width(range = c(1, 20)) +
  xlim(-1, 1) + ylim(-1, 1) + # Adjust edge width scaling
  theme_void() +  # Remove background
  ggtitle("Europe")

####################### North America Mantel Graph ############################

edges_North_America <- data.frame(
  from = c("Phylogeny", "Morphology", "Guilds Bray"),
  to = c("Guilds Bray", "Phylogeny", "Morphology"),
  weight = c(mantel_North_America_P_G1[["statistic"]], mantel_North_America_P_M[["statistic"]], mantel_North_America_M_G1[["statistic"]])
)

triangle_graph_North_America <- graph_from_data_frame(edges_North_America, vertices = nodes, directed = FALSE)

ggraph(triangle_graph_North_America, layout = "manual", x = nodes$x, y = nodes$y) + 
  geom_edge_link(aes(width = weight, label = round(weight, 3)), 
                 color = "darkkhaki", edge_alpha = 0.8, 
                 label_size = 5, show.legend = FALSE) +
  geom_node_point(size = 40, color = "cyan4") +  # Node size
  geom_node_text(aes(label = name), size = 5, color = "black") +  # Labels inside nodes
  scale_edge_width(range = c(1, 10)) +
  xlim(-1, 1) + ylim(-1, 1) + 
  theme_void() +
  ggtitle("North_America")

##################################### Guild visualization Global #######################################################

bray_tree_Global  <- as.phylo(dendro_meta_bray)
bray_tree_Asia  <- as.phylo(dendro_Asia_bray)
bray_tree_Australia  <- as.phylo(dendro_Australia_bray)
bray_tree_Europe  <- as.phylo(dendro_Europe_bray)
bray_tree_North_America  <- as.phylo(dendro_North_America_bray)

trait_labels<-c("Flycatch", "Glean", "Hover", "Pounce", "Probe", "Snatch", "Air", "Bark", "Flower", "Ground", "Leaf")



traits_Global <- phylo4d( x=bray_tree_Global, tip.data=matrix_Global_prop )
traits_Asia <- phylo4d( x=bray_tree_Asia, tip.data=matrix_Asia_prop )
traits_Australia <- phylo4d( x=bray_tree_Australia, tip.data=matrix_Australia_prop )
traits_Europe <- phylo4d( x=bray_tree_Europe, tip.data=matrix_Europe_prop )
traits_North_America <- phylo4d( x=bray_tree_North_America, tip.data=matrix_North_America_prop )


par(oma=c(6,2,2,0))
table.phylo4d(traits_Global, treetype="phylogram", symbol="circles", ratio.tree=0.2, center=F, scale=F, legend=F, grid=T, box=F, cex.symbol=0.3, cex.label=0.6, cex.legend=0.8, var.label=trait_labels, main="Guilds Global")

table.phylo4d(traits_Asia, treetype="phylogram", symbol="circles", ratio.tree=0.2, center=F, scale=F, legend=F, grid=T, box=F, cex.symbol=0.3, cex.label=0.6, cex.legend=0.8, var.label=trait_labels, main="Guilds Asia")

gridplot.phylo4d(traits_Global, tree.ladderize=T, center=F, scale=F, tree.type="phylogram", tree.ratio=0.15, trait.bg.col = "white", show.box = T, trait.labels = trait_labels, main="Guilds Global", cex.main=1.2, cell.col = white2red(200))

################## Per-Continent Guild visualization ############################
gridplot.phylo4d(traits_Asia, tree.ladderize=T, center=F, scale=F, tree.type="phylogram", tree.ratio=0.15, trait.bg.col = "white", show.box = T, trait.labels = trait_labels, main="Guilds Asia", cex.main=1.2, cell.col = white2red(200))

gridplot.phylo4d(traits_Australia, tree.ladderize=T, center=F, scale=F, tree.type="phylogram", tree.ratio=0.15, trait.bg.col = "white", show.box = T, trait.labels = trait_labels, main="Guilds Australia", cex.main=1.2, cell.col = white2red(200))

gridplot.phylo4d(traits_Europe, tree.ladderize=T, center=F, scale=F, tree.type="phylogram", tree.ratio=0.15, trait.bg.col = "white", show.box = T, trait.labels = trait_labels, main="Guilds Europe", cex.main=1.2, cell.col = white2red(200))

gridplot.phylo4d(traits_North_America, tree.ladderize=T, center=F, scale=F, tree.type="phylogram", tree.ratio=0.15, trait.bg.col = "white", show.box = T, trait.labels = trait_labels, main="Guilds North America", cex.main=1.2, cell.col = white2red(200))




