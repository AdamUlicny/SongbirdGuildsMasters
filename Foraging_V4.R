# Foraging Guilds and Specialization of Songbirds in a Czech Lowland Deciduous Forest

################################################################################
# 1. Load Libraries
################################################################################

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
library(ggrepel)
library(gt)
library(scales)

################################################################################
# 2. Load Data
################################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data_23 <- read_excel("./data/behav_data_23.xlsx") # Behavioral data 2023
data_24 <- read_excel("./data/behav_data_24.xlsx") # Behavioral data 2024
data_bodovka <- read_excel("./data/bodovka_data_23.xlsx") # Point transect data 2023

################################################################################
# 3. Data Preparation and Cleaning
################################################################################

# 3.1. Clean 'data_bodovka'
data_bodovka <- data_bodovka %>%
  unite(col = "sp_orig", Genus, Species, sep = "_", remove = TRUE) # Combine Genus and Species

data_bodovka$sp_orig <- gsub("_", " ", data_bodovka$sp_orig) # Replace underscores with spaces
data_bodovka$sp_orig <- gsub("Carduelis chloris", "Chloris chloris", data_bodovka$sp_orig) # Correct species names
data_bodovka$sp_orig <- gsub("Phoenicorus phoenicorus", "Phoenicurus phoenicurus", data_bodovka$sp_orig)
data_bodovka$sp_orig <- gsub("Parus palustris", "Poecile palustris", data_bodovka$sp_orig)

# 3.2. Combine and clean 'data_23' and 'data_24'
data_23 <- data_23 %>%
  mutate(line = 2023) %>%
  unite(col = "ID", ID, line, sep = "_", remove = TRUE) # Add year and combine ID

data_24 <- data_24 %>%
  mutate(line = 2024) %>%
  unite(col = "ID", ID, line, sep = "_", remove = TRUE) # Add year and combine ID

data_cz <- bind_rows(data_23, data_24) %>%
  unite(col = "sp_orig", genus, species, sep = "_", remove = TRUE) # Combine data and species names

# 3.3. Prepare long format data for method and substrate analysis
data_method_cz <- data_cz %>%
  filter(!is.na(behavior_1)) %>%
  pivot_longer(cols = starts_with("behavior"), names_to = "x", values_to = "behav") %>%
  select(sp_orig, behav, ID) %>%
  na.omit()

data_substrate_cz <- data_cz %>%
  filter(!is.na(behavior_1)) %>%
  pivot_longer(cols = starts_with("substrate_main"), names_to = "x", values_to = "substrate") %>%
  select(sp_orig, substrate, ID) %>%
  na.omit()

data_cz_long <- bind_cols(data_method_cz, data_substrate_cz)%>%
  mutate(line=1)

data_cz_long <- data_cz_long %>%
  select(ID...3, sp_orig...1, behav, substrate,line)%>%
  rename(ID = ID...3, sp_orig = sp_orig...1)

data_cz_long$sp_orig <- gsub("_", " ", data_cz_long$sp_orig) # Clean species names
data_cz_long$sp_orig <- gsub("Parus caeruleus", "Cyanistes caeruleus", data_cz_long$sp_orig)
data_cz_long$sp_orig <- gsub("Carduelis chloris", "Chloris chloris", data_cz_long$sp_orig)
data_cz_long$sp_orig <- gsub("Phoenicorus phoenicorus", "Phoenicurus phoenicurus", data_cz_long$sp_orig)
data_cz_long$sp_orig <- gsub("Parus palustris", "Poecile palustris", data_cz_long$sp_orig)

# 3.4. Filter species based on counts and remove unwanted taxa
count_actions_sp <- data_cz_long %>% group_by(sp_orig) %>% summarise(n = n())
count_individuals_sp <- data_cz_long %>% group_by(sp_orig) %>% summarise(n = n_distinct(ID))
counts_sp <- left_join(count_actions_sp, count_individuals_sp, by = "sp_orig") %>% rename(actions = n.x, individuals = n.y)

filter_list <- counts_sp %>%
  filter(actions < 20 | individuals < 5) %>%
  pull(sp_orig) %>%
  c("Dendrocopos major", "Dendrocopos minor", "Dryocopus martius", "Columba palumbus", "Columba oenas", "Cuculus canorus", "Buteo buteo", "Corvus corax", "Garrulus glandarius")

data_cz_long <- data_cz_long %>% filter(!sp_orig %in% filter_list)
data_cz <- data_cz %>% filter(!sp_orig %in% filter_list)

################################################################################
# 4. Descriptive Statistics and Visualization
################################################################################

# 4.1. Species Frequency
frequency_cz <- data_cz_long %>%
  group_by(sp_orig) %>%
  summarise(n = n()) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  arrange(desc(percentage))

frequency_bodovka <- data_bodovka %>%
  group_by(sp_orig) %>%
  filter(!sp_orig %in% filter_list) %>%
  summarise(n = n()) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  arrange(desc(percentage))

# 4.2. Frequency Plots
ggplot(data = frequency_cz, aes(x = reorder(sp_orig, -percentage), y = percentage)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Species", y = "Percentage of all observations (%)") +
  theme_minimal()

ggplot(data = frequency_cz %>% filter(sp_orig %in% frequency_bodovka$sp_orig), aes(x = reorder(sp_orig, -percentage), y = percentage, fill = "Behavioral observations")) +
  geom_bar(stat = "identity") +
  geom_bar(data = frequency_bodovka %>% filter(sp_orig %in% frequency_cz$sp_orig), aes(x = reorder(sp_orig, -percentage), y = percentage, fill = "Point transect"), stat = "identity", alpha = 0.5) +
  ylim(0, 20) +
  coord_flip() +
  labs(x = "", y = "Frequency of observations in %", fill = "Dataset type") +
  scale_fill_manual(values = c("Behavioral observations" = "steelblue", "Point transect" = "red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = c(0.6, 0.9), legend.justification = c(0, 1))


# 4.3. Species not present in behavioral data
missed_species <- data_bodovka %>%
  filter(!sp_orig %in% data_cz_long$sp_orig) %>%
  filter(!sp_orig %in% filter_list) %>%
  select(sp_orig) %>%
  distinct()

counts_bodovka <- data_bodovka %>%
  group_by(sp_orig) %>%
  summarise(n = n()) %>%
  filter(!sp_orig %in% filter_list)

write.csv(counts_sp$sp_orig, file = "./resources/sp_orig.csv")

# 4.4 Stacked Bar Charts for Method and Substrate
graph_method <- ggplot(data_cz_long, aes(x = line, fill = behav)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title = element_blank(), axis.text.x = element_blank(), legend.background = element_rect(fill = 'transparent'), axis.ticks.x = element_blank()) +
  geom_bar(position = "fill", color = "black") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(values = c(flycatch = "#0cf0e8", glean = "#06c24b", hang_glean = "#8B7500", hover_snatch = "#9606c2", manipulation = "#0677c2", pounce = "#2F4F4F", probe = "#c2b906", snatch = "#ed8105"), labels = c("Flycatch", "Glean", "Hang-Glean", "Hover-Snatch", "Manipulation", "Pounce", "Probe", "Snatch")) +
  labs(x = "Method", y = "Method/Substrate usage in %", title = "")

graph_substrate <- ggplot(data_cz_long, aes(x = line, fill = substrate)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title = element_blank(), axis.text.x = element_blank(), legend.background = element_rect(fill = 'transparent'), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_bar(position = "fill", color = "black") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(values = c(air = "#AEEEEE", bark = "#F5DEB3", ground = "#8B7500", leaf = "#556B2F", other = "#2F5F8F"), labels = c("Air", "Bark", "Ground", "Leaf", "Other")) +
  labs(x = "Substrate", y = "", title = "")

graph_foliage <- data_cz %>%
  filter(!is.na(dist_stem)) %>%
  drop_na(foliage_dens) %>%
  ggplot(aes(x = plot, fill = factor(foliage_dens, levels = c("low", "medium", "high")))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  geom_bar(position = "fill", color = "black") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Greens"), labels = c("Low", "Medium", "High")) +
  labs(x = "Foliage density", y = "Foliage preference/vegetation position in %", title = "")

graph_distance <- data_cz %>%
  drop_na(dist_stem) %>%
  ggplot(aes(x = plot, fill = factor(dist_stem, levels = c("edge", "outer", "inner", "stem")))) +
  theme(element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  geom_bar(position = "fill", color = "black") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Oranges"), labels = c("Edge", "Outer", "Inner", "Stem")) +
  labs(x = "Vegetation position", y = "", title = "")

graph_method + graph_substrate + graph_foliage + graph_distance



# 4.5. Counts and percentages of method and substrate
method_counts <- data_cz_long %>% count(behav) %>% mutate(percentage = n / sum(n) * 100) %>% arrange(desc(percentage))
substrate_counts <- data_cz_long %>% count(substrate) %>% mutate(percentage = n / sum(n) * 100) %>% arrange(desc(percentage))
veget_counts <- data_cz_long %>% filter(substrate == "bark" | substrate == "leaf") %>% count(substrate) %>% mutate(percentage = n / sum(n) * 100) %>% arrange(desc(percentage))
distance_counts <- data_cz %>% drop_na(dist_stem) %>% count(dist_stem) %>% mutate(percentage = n / sum(n) * 100) %>% arrange(desc(percentage))
foliage_counts <- data_cz %>% drop_na(foliage_dens) %>% count(foliage_dens) %>% mutate(percentage = n / sum(n) * 100) %>% arrange(desc(percentage))

data_cz_long %>% count(behav, substrate) %>% mutate(percentage = n / sum(n) * 100) %>% arrange(desc(percentage))

# 4.6. Connections Graph
edges <- data_cz_long %>% count(behav, substrate, name = "weight") %>% mutate(weight = weight / sum(weight) * 100)

substrate_colors <- c(air = "#AEEEEE", bark = "#F5DEB3", ground = "#8B7500", leaf = "#556B2F", other = "#2F4F4F")
behav_colors <- c(flycatach = "#0cf0e8", glean = "#06c24b", hang_glean = "#8B7500", hover_snatch = "#9606c2", manipulation = "#0677c2", pounce = "#2F4F4F", probe = "#c2b906", snatch = "#ed8105")

graph_connections <- graph_from_data_frame(edges, directed = FALSE)
node_type <- ifelse(V(graph_connections)$name %in% edges$behav, "Behavior", "Substrate")
node_colors <- ifelse(node_type == "Behavior", "#F8766D", "#00BFC4")


ggplot(edges, aes(axis1 = behav, axis2 = substrate, y = weight)) +
  geom_alluvium(aes(), width = 1/8, alpha = 0.9, knot.pos = 0.3) +
  geom_stratum(aes(fill = behav), width = 1/6) +
  geom_stratum(aes(fill = substrate), width = 1/6) +
  scale_fill_manual(values = c(behav_colors, substrate_colors), name = "Category") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 14), axis.ticks = element_blank(), legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  labs(title = "Behavior-Substrate Co-occurrence", x = NULL, y = "Percentage (%)")

################################################################################
# 5. Levins Specialization Index
################################################################################

calculate_index <- function(data, select_column, n_categories = 5) {
  data %>%
    group_by(sp_orig) %>%
    count({{ select_column }}) %>%
    mutate(prop_category = prop.table(n)) %>%
    mutate(pi2 = prop_category^2) %>%
    summarise(B = 1 / sum(pi2), .groups = 'drop') %>%
    mutate(Ba = 1 - (B - 1) / (n_categories - 1))
}

levins_method <- calculate_index(data_cz_long, behav, 8) %>% rename(Ba_method = Ba, B_method = B)
levins_substrate <- calculate_index(data_cz_long, substrate, 5) %>% rename(Ba_substrate = Ba, B_substrate = B)
levins_method_substrate <- left_join(levins_method, levins_substrate, by = "sp_orig")

graph_specialization <- ggplot(levins_method_substrate, aes(x = Ba_method, y = Ba_substrate)) +
  geom_point(size = 3) +
  labs(x = "Method specialization", y = "Substrate specialization") +
  geom_smooth(method = "lm", color = "black", se = FALSE, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = sp_orig), size = 5, force = 1000, segment.size = 0.5, segment.color = "grey50", max.overlaps = Inf) +
  theme_classic() +
  scale_x_continuous(labels = label_comma(decimal.mark = ",")) +
  scale_y_continuous(labels = label_comma(decimal.mark = ",")) +
  theme(axis.title = element_text(size = 20))

plot(graph_specialization)

summary(lm(Ba_substrate ~ Ba_method, levins_method_substrate))
cor.test(levins_method_substrate$Ba_method, levins_method_substrate$Ba_substrate, method = "pearson")

table_specialization <- levins_method_substrate %>%
  select(sp_orig, Ba_substrate, Ba_method) %>%
  rename("Latin species name" = sp_orig, "Substrate specialization" = Ba_substrate, "Method specialization" = Ba_method) %>%
  gt() %>%
  tab_header(title = "") %>%
  fmt_number(columns = c("Substrate specialization", "Method specialization"), decimals = 2) %>%
  cols_align(align = "center") %>%
  tab_options(table.font.size = "medium")

table_specialization

################################################################################
# 6. Dissimilarity Matrix and Dendrograms
################################################################################

data_cz_wide <- data_cz_long %>%
  group_by(sp_orig) %>%
  count(sp_orig, behav, substrate, sort = TRUE) %>%
  unite(col = "behav_substrate", behav, substrate, sep = "-", remove = TRUE) %>%
  pivot_wider(names_from = "behav_substrate", values_from = "n") %>%
  replace(is.na(.), 0) %>%
  remove_rownames() %>%
  column_to_rownames(var = "sp_orig")

bray_dis_matrix <- vegdist(data_cz_wide, method = "bray") %>% as.matrix() %>% as.dist()
jaccard_dis_matrix <- vegdist(data_cz_wide, method = "jaccard") %>% as.matrix() %>% as.dist()

phylo_cz <- phylo_cz[[10]] %>% ape::drop.tip(filter_list) %>% cophenetic() %>% as.matrix() %>% as.dist()

dendro_bray <- hclust(bray_dis_matrix, method = "ward.D2") %>% as.dendrogram()
dendro_jaccard <- hclust(jaccard_dis_matrix, method = "ward.D2") %>% as.dendrogram()
dendro_phylo <- hclust(phylo_cz) %>% as.dendrogram()

plot(dendro_bray, horiz = TRUE, main = "Dendrogram Bray-Curtis dissimilarity")
plot(dendro_jaccard, horiz = TRUE, main = "Dendrogram Jaccard distance")
plot(dendro_phylo, horiz = TRUE, main = "Phylogeny")

clusMember <- cutree(dendro_bray, 5)
labelColors <- c("#D46B37", "#020C45", "#036564", "#7A67EE", "#D43B22")
labelLegend <- c("Bark", "Probe", "Leaf", "Air", "Ground")

colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}

dendro_bray_aes <- dendrapply(dendro_bray, colLab)

par(mar = c(5, 1, 1, 12))
plot(dendro_bray_aes, main = "", type = "rectangle", horiz = TRUE, xlab = "Bray-Curtis distance")
legend("topleft", legend = labelLegend, col = labelColors, pch = c(20, 20, 20, 20), bty = "y", pt.cex = 1.5, cex = 1.2, text.col = "black", horiz = FALSE, title = "Guilds", inset = c(0, 0.05))

################################################################################
# 7. Tanglegram
################################################################################

set.seed(12345)
dendlist(dendro_bray, dendro_phylo) %>%
  dendextend::untangle(method = "random", R = 100) %>%
  dendextend::untangle(method = "step2side") %>%
  tanglegram(common_subtrees_color_lines = TRUE, highlight_distinct_edges = FALSE, highlight_branches_lwd = FALSE, margin_inner = 10, margin_outer = 7, lwd = 3, main_left = "behavior", main_right = "phylogeny", hang = FALSE)

dendlist(dendro_bray, dendro_jaccard) %>%
  dendextend::untangle(method = "random", R = 100) %>%
  dendextend::untangle(method = "step2side") %>%
  tanglegram(common_subtrees_color_lines = TRUE, highlight_distinct_edges = FALSE, highlight_branches_lwd = FALSE, margin_inner = 10, margin_outer = 7, lwd = 3, main_left = "Bray-Curtis", main_right = "Jaccard", hang = FALSE)

################################################################################
# 8. Phylogenetic Signal
################################################################################

Bray_ <- as.phylo(dendro_Europe_bray)
data_cz_prop <- data_cz_wide %>% mutate_all(~ . / sum(.))
trait_labels <- c("Flycatch", "Glean", "Hover", "Pounce", "Probe", "Snatch", "Air", "Bark", "Flower", "Ground", "Leaf")
traits_Europe <- phylo4d(x = bray_tree_Europe, tip.data = matrix_Europe_prop)

gridplot.phylo4d(traits_Europe, tree.ladderize = TRUE, center = FALSE, scale = FALSE, tree.type = "phylogram", tree.ratio = 0.15, trait.bg.col = "white", show.box = TRUE, trait.labels = trait_labels, main = "Guilds Europe", cex.main = 1.2, cell.col = white2red(200))

