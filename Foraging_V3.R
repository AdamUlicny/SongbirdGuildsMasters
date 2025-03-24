######### Foraging Guilds and Specialization of songbirds in a Czech Lowland Deciduous forest ############
################################# Libraries ##############################################################
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
library(igraph)
library(ggraph)
library(ggalluvial)

################################### Load data ###########################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data_23 <- read_excel("./data/behav_data_23.xlsx")
data_24 <- read_excel("./data/behav_data_24.xlsx")
data_bodovka <- read_excel("./data/bodovka_data_23.xlsx")
phylo_cz <- "./resources/phylo/phylo_cz_v1/output.nex"
phylo_cz <- ape::read.nexus(phylo_cz)
################################### Data preparation#####################
# in data bodovka unite columns genus and species, separator _, name sp_orig
data_bodovka <- data_bodovka %>%
  unite(col = "sp_orig", Genus, Species, sep = "_", remove = TRUE)

data_bodovka$sp_orig<-gsub("_", " ", data_bodovka$sp_orig)
data_bodovka$sp_orig <- gsub("Carduelis chloris", "Chloris chloris", data_bodovka$sp_orig)
data_bodovka$sp_orig <- gsub("Phoenicorus phoenicorus", "Phoenicurus phoenicurus", data_bodovka$sp_orig)
data_bodovka$sp_orig <- gsub("Parus palustris", "Poecile palustris", data_bodovka$sp_orig)

# In data_23 replace line column with year 2023, then combine ID with year, separate 
data_23 <- data_23 %>%
  mutate(line = 2023) %>%
  unite(col = "ID", ID, line, sep = "_", remove = TRUE)

# same for data_24
data_24 <- data_24 %>%
  mutate(line = 2024) %>%
  unite(col = "ID", ID, line, sep = "_", remove = TRUE)

# combine data_23 and data_24
data_cz <- bind_rows(data_23, data_24)

# Method - remove NA, unite genus and species, pivot longer
data_method_cz <- data_cz %>%
  filter(!is.na(behavior_1)) %>%
  unite(col = "sp_orig", genus, species, sep = "_", remove = TRUE) %>%
  pivot_longer(cols = starts_with("behavior"), 
               names_to = "x", values_to = "behav") %>%
  select(sp_orig, behav, ID ) %>%
  na.omit()

# Substrate - remove NA, unite genus and species, pivot longer
data_substrate_cz <- data_cz %>%
  filter(!is.na(behavior_1)) %>%
  unite(col = "sp_orig", genus, species, sep = "_", remove = TRUE) %>%
  pivot_longer(cols = starts_with("substrate_main"),
               names_to = "x", values_to = "substrate") %>%
  select(sp_orig, substrate, ID) %>%
  na.omit(F)

# Merge data_method_cz and data_substrate_cz
data_cz_long <- bind_cols(data_method_cz, data_substrate_cz)%>%
  mutate(line=1)

data_cz_long <- data_cz_long %>%
  select(ID...3, sp_orig...1, behav, substrate,line)%>%
  rename(ID = ID...3, sp_orig = sp_orig...1)

# remove "_" and replace with " " in sp_orig
data_cz_long$sp_orig <- gsub("_", " ", data_cz_long$sp_orig)

# replace all "Parus caeruleus" with "Cyanistes caeruleus"
data_cz_long$sp_orig <- gsub("Parus caeruleus", "Cyanistes caeruleus", data_cz_long$sp_orig)
data_cz_long$sp_orig <- gsub("Carduelis chloris", "Chloris chloris", data_cz_long$sp_orig)
data_cz_long$sp_orig <- gsub("Phoenicorus phoenicorus", "Phoenicurus phoenicurus", data_cz_long$sp_orig)
data_cz_long$sp_orig <- gsub("Parus palustris", "Poecile palustris", data_cz_long$sp_orig)


# count numbers of actions and individuals per species
count_actions_sp <- data_cz_long%>%
  group_by(sp_orig) %>%
  summarise(n = n())

count_individuals_sp <- data_cz_long %>%
  group_by(sp_orig) %>%
  summarise(n = n_distinct(ID))

# merge count_actions_sp and count_individuals_sp
counts_sp <- left_join(count_actions_sp, count_individuals_sp, by = "sp_orig")
  
counts_sp <- counts_sp%>%
  rename(actions = n.x, individuals = n.y)

# create filtering list for easy species removal actions  and individuals < 3 + remove Piciformes 
filter_list <- counts_sp %>%
  filter(actions < 10 | individuals < 3) %>%
  pull(sp_orig)%>%
  c("Dendrocopos_major", "Dendrocopos_minor", "Dryocopus_martius")

# Filter out species from filter_list and provide final list of passerines (n=19)
data_cz_long <- data_cz_long %>%
  filter(!sp_orig %in% filter_list)

# Table with species and percentage of all observations
frequency_cz <- data_cz_long %>%
  group_by(sp_orig) %>%
  summarise(n = n()) %>%
  mutate(percentage = n/sum(n)*100) %>%
  arrange(desc(percentage))

# abundance curve histogram per species
ggplot(data = frequency_cz, aes(x = reorder(sp_orig, -percentage), y = percentage)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Species", y = "Percentage of all observations (%)") +
  theme_minimal()

# abundance curve histogram per species in data_bodovka
frequency_bodovka <- data_bodovka %>%
  group_by(sp_orig) %>%
  summarise(n = n()) %>%
  mutate(percentage = n/sum(n)*100) %>%
  arrange(desc(percentage))

# compare histograms of frequency per species in data_cz and data_bodovka showing only species shared by both tables
ggplot(data = frequency_cz %>% filter(sp_orig %in% frequency_bodovka$sp_orig), aes(x = reorder(sp_orig, -percentage), y = percentage)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Species", y = "Percentage of all observations (%)") +
  theme_minimal() +
  geom_bar(data = frequency_bodovka %>% filter(sp_orig %in% frequency_cz$sp_orig), aes(x = reorder(sp_orig, -percentage), y = percentage), stat = "identity", fill = "red", alpha = 0.5)

# species in data_bodovka, not present in passerines_cz 
missed_species <- data_bodovka %>%
  filter(!sp_orig %in% data_cz_long$sp_orig) %>%
  select(sp_orig) %>%
  distinct()

# sp_orig in data_cz_long to csv
write.csv(counts_sp$sp_orig, file = "./resources/sp_orig.csv")

######################### Basic Graphs #####################################
# stacked bar chart of method, substrate, foliage cover and bird height
graph_method<-ggplot(data_cz_long, aes(x = line, fill = behav)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title = element_blank(), axis.text.x=element_blank(), legend.background = element_rect(fill='transparent'),
        axis.ticks.x=element_blank())+
  geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(
    values = c(flycatach = "#0cf0e8",
               glean = "#06c24b",
               hang_glean = "#8B7500",
               hover_snatch = "#9606c2",
               manipulation = "#0677c2",
               pounce = "#2F4F4F",
               probe = "#c2b906",
               snatch = "#ed8105"))+
  labs(x="", y="", title = "Foraging method")

graph_substrate<-ggplot(data_cz_long, aes(x = line, fill = substrate)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank(),axis.text.x=element_blank(),legend.background = element_rect(fill='transparent'), 
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank())+
  geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(
    values = c(air = "#AEEEEE",
               bark = "#F5DEB3",
               ground = "#8B7500",
               leaf = "#556B2F",
               other = "#2F4F4F"))+
  labs(x="", y="", title = "Foraging substrate")

graph_method+graph_substrate




graph_foliage<- data_cz%>% 
  filter(!is.na(dist_stem))%>%
  drop_na(foliage_dens)%>%
  ggplot(aes(x = plot, fill = factor(foliage_dens, levels=c("low", "medium", "high")))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank())+
  geom_bar(position="fill", color= "black")+
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_brewer(palette = "Greens")+
  labs(x="", y="", title = "Foliage density")

graph_distance<-data_cz%>% 
  drop_na(dist_stem)%>%
  ggplot(aes(x = plot, fill = factor(dist_stem, levels=c("edge", "outer", "inner", "stem")))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank())+
  geom_bar(position="fill", color= "black")+
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_brewer(palette = "Oranges")+
  labs(x="", y="", title = "Bird distance from stem")

graph_foliage+graph_distance


grid.arrange(graph_method, graph_substrate, graph_foliage, graph_distance, ncol = 2)


###### Connections graph between method-substrate
edges<- data_cz_long%>%
  count(behav, substrate, name="weight")
edges <- edges %>%
  mutate(weight = weight / sum(weight) * 100)
substrate_colors <- c(air = "#AEEEEE", bark = "#F5DEB3", ground = "#8B7500", 
                      leaf = "#556B2F", other = "#2F4F4F")

behav_colors <- c(flycatach = "#0cf0e8", glean = "#06c24b", hang_glean = "#8B7500",
                  hover_snatch = "#9606c2", manipulation = "#0677c2", pounce = "#2F4F4F", 
                  probe = "#c2b906", snatch = "#ed8105")


graph_connections <- graph_from_data_frame(edges, directed = FALSE)

node_type <- ifelse(V(graph_connections)$name %in% edges$behav, "Behavior", "Substrate")
node_colors <- ifelse(node_type == "Behavior", "#F8766D", "#00BFC4")  # Red for behaviors, Blue for substrates

ggraph(graph_connections, layout = "fr") + 
  geom_edge_link(aes(edge_alpha = weight, edge_width = weight), 
                 edge_color = "gray70", curvature = 0.3) + 
  geom_node_point(aes(color = node_type), size = 6) +  
  geom_node_text(aes(label = name), repel = TRUE, size = 5) +  
  scale_edge_alpha(range = c(0.3, 1)) +  
  scale_edge_width(range = c(0.5, 2)) +  
  scale_color_manual(values = c("Behavior" = "#F8766D", "Substrate" = "#00BFC4")) +  
  theme_void() +  
  theme(legend.position = "none")


ggplot(edges, aes(axis1 = behav, axis2 = substrate, y = weight)) +
  geom_alluvium(aes(fill = behav), width = 1/8, alpha = 0.9, knot.pos = 0.3) + 
  geom_stratum(aes(fill = behav), width = 1/6) +  
  geom_stratum(aes(fill = substrate), width = 1/6) +  
  scale_fill_manual(values = c(behav_colors, substrate_colors), name = "Category") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) + # Convert y-axis to percentages
  theme_minimal() +  
  theme(panel.grid = element_blank(),  # Remove gridlines
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(size = 14),
        axis.ticks = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(title = "Behavior-Substrate Co-occurrence",
       x = NULL, y = "Percentage (%)")


############# Levins Specialization index ##########
# 𝑩= 𝟏Σ𝒑𝒊𝟐
# 𝑩a=𝟏−(𝑩−𝟏)/(𝒏−𝟏) (where n=number of available categories to specialize in)
calculate_index <- function(data, select_column, n_categories = 5) {
  data %>%
    group_by(sp_orig) %>%
    count({{ select_column }}) %>%
    mutate(prop_category = prop.table(n)) %>%
    mutate(pi2 = prop_category^2) %>%
    summarise(B = 1 / sum(pi2), .groups = 'drop') %>%
    mutate(Ba = 1 - (B - 1) / (n_categories - 1))
}


# Calculate specialization index for foraging method
levins_method <- calculate_index(data_cz_long, behav, 7)%>%
  rename(Ba_method = Ba, B_method = B)

levins_substrate <- calculate_index(data_cz_long, substrate, 5)%>%
  rename(Ba_substrate = Ba, B_substrate = B)

# merge specialization indexes
levins_method_substrate <- left_join(levins_method, levins_substrate, by = "sp_orig")

# Scatterplot comparing method/substrate specialization
par(mar = c(5, 5, 5, 5))
plot(Ba_substrate ~ Ba_method, xlim = c(0.3, 1), ylim = c(0.3, 1), ylab="Substrate specialization",xlab="Method specialization", data=levins_method_substrate, cex.lab=2, pch=19)
abline(lm(Ba_substrate ~ Ba_method, levins_method_substrate), lw=1.3)
abline(c(0,1), lty=2, col="red")

# ggplot
graph_specialization <- ggplot(levins_method_substrate, aes(x = Ba_method, y = Ba_substrate)) +
  geom_point(size = 3) +
  labs(x = "Method specialization", y = "Substrate specialization") +
  xlim(0.3, 1) +
  ylim(0.3, 1) +
  geom_text(aes(label = sp_orig), vjust = -0.5, hjust = 0.5, size = 3) +  # Add species labels
  geom_smooth(method = "lm", color = "black", se = F, size = 1) +  # Line of best fit
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Diagonal reference line
  theme_classic() +
  theme(axis.title = element_text(size = 16))
plot(graph_specialization)

# linear regression
summary(lm(Ba_substrate ~ Ba_method, levins_method_substrate))

################## Dissimilarity matrix ############################
# prepare data by pivoting wider
data_cz_wide <- data_cz_long %>%
  group_by(sp_orig) %>% 
  count(sp_orig,behav,substrate, sort=TRUE)%>% 
  unite(col = "behav_substrate", behav, substrate, sep = "-", remove = TRUE) %>%
  pivot_wider(names_from="behav_substrate",values_from="n")%>%
  replace(is.na(.), 0)%>%
  remove_rownames%>% 
  column_to_rownames(var="sp_orig")

# calculate dissimilarity/distance matrix
bray_dis_matrix <- vegdist(data_cz_wide, method = "bray")
bray_dis_matrix <- as.matrix(bray_dis_matrix)
bray_dis_matrix <- as.dist(bray_dis_matrix[order(rownames(bray_dis_matrix)),order(colnames(bray_dis_matrix))])

jaccard_dis_matrix <- vegdist(data_cz_wide, method = "jaccard")
jaccard_dis_matrix <- as.matrix(jaccard_dis_matrix)
jaccard_dis_matrix <- as.dist(jaccard_dis_matrix[order(rownames(jaccard_dis_matrix)),order(colnames(jaccard_dis_matrix))])

# prepare phylogenetic distances
phylo_cz <- phylo_cz[[10]]
phylo_cz <- ape::drop.tip(phylo_cz, filter_list)
phylo_cz <- cophenetic(phylo_cz)
phylo_cz <- as.matrix(phylo_cz)
phylo_cz <- as.dist(phylo_cz[order(rownames(phylo_cz)),order(colnames(phylo_cz))])



# Prepare dendrograms
dendro_bray <- bray_dis_matrix %>% 
  hclust(method = "ward.D2") %>%
  as.dendrogram

dendro_jaccard <- jaccard_dis_matrix %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram

dendro_phylo <- hclust(phylo_cz) %>%
  as.dendrogram

############# Plotting dendrograms ###############
# Bray dissimilarity dendrogram
plot(dendro_bray,horiz=T, main="Dendrogram Bray-Curtis dissimilarity", )
plot(dendro_jaccard,horiz=T, main="Dendrogram Jaccard distance", )
plot(dendro_phylo,horiz=T, main="Phylogeny", )

######### Colored guilds ##############

# cut dendrogram into 5 guilds
clusMember = cutree(dendro_bray, 5)

# vector of colors
labelColors = c("#D46B37", "#020C45","#036564","#7A67EE", "#D43B22")
labelLegend = c("bark", "probe", "leaf", "air", "ground")

# guild color function
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}

# dendrapply function to guilds
dendro_bray_aes = dendrapply(dendro_bray, colLab)

### plot guilds, colors and legend
par(mar=c(5,1,1,12))
plot(dendro_bray_aes, main = "Foraging guilds", type = "rectangle", horiz = T)
legend("topleft", 
       legend = labelLegend, 
       col = labelColors, 
       pch = c(20,20,20,20), bty = "y",  pt.cex = 1.5, cex = 0.8 , 
       text.col = "black", horiz = FALSE,
       title="Guilds",
       inset = c(0, 0.1))

#################### Tanglegram ############################

set.seed(12345)
dendlist(dendro_bray, dendro_phylo)%>%
  dendextend::untangle(method="random", R=100)%>%####### Crucial step to produce human readable codendrograms! Use lower R on slower machines.
  dendextend::untangle(method="step2side")%>%
  tanglegram(common_subtrees_color_lines = TRUE, # Do NOT include "sort=T" argument if using untangle before (sort overrides it)
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="behavior",
             main_right="phylogeny",
             hang=F)
mantel_CZ <- mantel(dist_Bray_North_America, phylo_North_America, method = "spearman", permutations = 9999)
print(mantel_North_America)
# Using Jaccard

dendlist(dendro_bray, dendro_jaccard)%>%
  dendextend::untangle(method="random", R=100)%>%####### Crucial step to produce human readable codendrograms! Use lower R on slower machines.
  dendextend::untangle(method="step2side")%>%
  tanglegram(common_subtrees_color_lines = TRUE, # Do NOT include "sort=T" argument if using untangle before (sort overrides it)
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=10,
             margin_outer=7,
             lwd=3,
             main_left="Bray-Curtis",
             main_right="Jaccard",
             hang=F)
