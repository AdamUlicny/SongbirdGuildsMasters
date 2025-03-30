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
  filter(!sp_orig %in% c("Dendrocopos major", "Dendrocopos minor", "Dryocopus martius"))
