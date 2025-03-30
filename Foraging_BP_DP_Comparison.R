##### Comparing 2021 results with 2023&2024 ##################
################################### Load data ###########################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data_21 <- read_excel("./data/behav_data_21.xlsx") 
data_23 <- read_excel("./data/behav_data_23.xlsx")
data_24 <- read_excel("./data/behav_data_24.xlsx")

################################### Data preparation ###########################
data_21 <- data_21 %>%
  filter(line == "2") %>%
  unite(col = "sp_orig", genus, species, sep = "_", remove = TRUE)


data_method_21 <- data_21 %>%
  filter(!is.na(behavior1)) %>%
  pivot_longer(cols = starts_with("behavior"), 
               names_to = "x", values_to = "behav") %>%
  select(sp_orig, behav) %>%
  na.omit()

data_substrate_21 <- data_21 %>%
  filter(!is.na(behavior1)) %>%
  pivot_longer(cols = starts_with("substrate_main"),
               names_to = "x", values_to = "substrate") %>%
  select(sp_orig, substrate) %>%
  na.omit(F)

# combine data_method_21 and data_substrate_21 without duplicating the sp_orig columns using bind_cols
data_long_21 <- bind_cols(data_method_21, data_substrate_21)

data_long_21 <- data_long_21%>%
  select(-c(sp_orig...3))%>%
  rename(sp_orig=sp_orig...1)

data_long_21$sp_orig<-gsub("_", " ", data_long_21$sp_orig)
data_long_21$sp_orig<-gsub("Phylloscopus sp.", "Phylloscopus collybita", data_long_21$sp_orig)
data_long_21$sp_orig<-gsub("Certhia sp.", "Certhia familiaris", data_long_21$sp_orig)
data_long_21$sp_orig<-gsub("Phylloscopus sp.", "Phylloscopus collybita", data_long_21$sp_orig)
data_long_21$sp_orig<-gsub("Erithacus rubecola", "Erithacus rubecula", data_long_21$sp_orig)
data_long_21$sp_orig<-gsub("Phylloscopus sp.", "Phylloscopus collybita", data_long_21$sp_orig)

data_long_21<-data_long_21%>%
  filter(!sp_orig %in% c("Dendrocopos major", "Dendrocopos minor", "Dryocopus martius", "Dendrocopos medius"))

data_wide_21 <- data_long_21 %>%
  group_by(sp_orig) %>% 
  count(sp_orig,behav,substrate, sort=TRUE)%>% 
  unite(col = "behav_substrate", behav, substrate, sep = "-", remove = TRUE) %>%
  pivot_wider(names_from="behav_substrate",values_from="n")%>%
  replace(is.na(.), 0)%>%
  remove_rownames%>% 
  column_to_rownames(var="sp_orig")

# calculate dissimilarity/distance matrix
bray_dis_matrix_21 <- vegdist(data_wide_21, method = "bray")
bray_dis_matrix_21 <- as.matrix(bray_dis_matrix_21)
bray_dis_matrix_21 <- as.dist(bray_dis_matrix_21[order(rownames(bray_dis_matrix_21)),order(colnames(bray_dis_matrix_21))])


data_long_23_24 <- data_cz_long%>%
  filter(sp_orig %in% data_long_21$sp_orig)

data_wide_23_24 <- data_long_23_24 %>%
  group_by(sp_orig) %>% 
  count(sp_orig,behav,substrate, sort=TRUE)%>% 
  unite(col = "behav_substrate", behav, substrate, sep = "-", remove = TRUE) %>%
  pivot_wider(names_from="behav_substrate",values_from="n")%>%
  replace(is.na(.), 0)%>%
  remove_rownames%>% 
  column_to_rownames(var="sp_orig")


bray_dis_matrix_23_24 <- vegdist(data_wide_23_24, method = "bray")
bray_dis_matrix_23_24 <- as.matrix(bray_dis_matrix_23_24)
bray_dis_matrix_23_24 <- as.dist(bray_dis_matrix_23_24[order(rownames(bray_dis_matrix_23_24)),order(colnames(bray_dis_matrix_23_24))])

mantel(bray_dis_matrix_21, bray_dis_matrix_23_24, method = "pearson", permutations = 999)

dendro_bray_21 <- bray_dis_matrix_21 %>% 
  hclust(method = "ward.D2") %>%
  as.dendrogram

dendro_bray_23_24 <- bray_dis_matrix_23_24 %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram

dendlist(dendro_bray_23_24, dendro_bray_21)%>%
  dendextend::untangle(method="random", R=100)%>%
  dendextend::untangle(method="step2side")%>%
  tanglegram(common_subtrees_color_lines = TRUE, # Do NOT include "sort=T" argument if using untangle before (sort overrides it)
             highlight_distinct_edges  = FALSE,
             highlight_branches_lwd=FALSE,
             margin_inner=13,
             margin_outer=7,
             lwd=3,
             main_left="Metodika 2023 & 2024",
             main_right="Metodika 2021",
             hang=F)

