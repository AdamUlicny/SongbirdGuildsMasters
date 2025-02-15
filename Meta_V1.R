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

list_sp<-sp_list_meta%>%filter(Passeriformes=="PASSERIFORMES")%>%pull(Sp_BirdLife)

# Load the phylogeny
phylo_meta <- read.nexus("resources/phylo/phylo_meta_v1/output.nex")

# filtering out non-passerines
method_substrate_meta<-method_substrate_meta%>%
  inner_join(sp_list_meta%>%filter(Passeriformes=="PASSERIFORMES"), by="Sp_BirdLife")

# in method_substrate_meta, deselect columns Study, Location and zdroj_dat, then sum up repeating Sp_BirdLife names and leave others untouched
method_substrate_meta<-method_substrate_meta%>%
  select(-c(Study, location, zdroj_dat))%>%
  group_by(Sp_BirdLife)%>%
  summarise(across(where(is.numeric), sum, na.rm = F))
