################### R skript DP ########################
setwd("E:/SYNC/DP")
require(tidyverse)
require(readxl)
data_all <- read_excel("behav_data_23.xlsx")
data_bodovka <- read_excel("bodovka_data_23.xlsx")

##################################################################

#odstranit radky bez forag action
data_forag <- data_all %>%
  filter(!is.na(behavior_1))

#sloucit sloupce genus a species
data_sp <- data_forag %>%
  unite(col = "druh", genus, species, sep = "_", remove = TRUE)

data_bodovka <- data_bodovka %>%
  unite(col = "druh", Genus, Species, sep= "_", remove =TRUE)

#pivot tabulky pro foraging
data_behav <- data_sp %>%
  pivot_longer(cols = c("behavior_1", "behavior_2", "behavior_3", "behavior_4", "behavior_5"), names_to = "x", values_to = "behav") %>%
  select(druh, behav, line ) %>%
  na.omit()

##pivot tabulky pro substrate
data_substrate <- data_sp %>%
  pivot_longer(cols = c("substrate_main_1", "substrate_main_2", "substrate_main_3", "substrate_main_4", "substrate_main_5"), names_to = "x", values_to = "substrate") %>%
  select(druh, line, substrate) %>%
  na.omit(F)

## pivot tabulky pro substrate_fine
data_substrate_fine <- data_sp %>%
  pivot_longer(cols = c("substrate_fine_1", "substrate_fine_2", "substrate_fine_3", "substrate_fine_4", "substrate_fine_5"), names_to = "x", values_to = "substrate_fine", values_drop_na = T) %>%
  select(druh, line, substrate_fine) 

## pivot tabulky pro druhy bodovka
data_bodovka_sp <- data_bodovka %>%
  group_by(Datum, druh)%>%
  filter(!(druh=="Columba_oenas"| druh=="Columba_palumbus"| druh=="Buteo_buteo"| druh=="Cuculus_canorus"| druh=="Corvus_corax"| druh=="Dryocopus_martius"| druh=="Dendrocopos_medius"| druh=="Dendrocopos_major"| druh=="Garrulus_glandarius"))%>%
  count()

#######################################################################

## square plot pro presence/absence z bodovky
ggplot(data_bodovka_sp) +
  aes(x = Datum, y = druh, fill = n) +
  geom_tile(size = 1.2) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous(labels = as.Date(data_bodovka_sp$Datum), breaks = data_bodovka_sp$Datum)

## bubble plot

ggplot(data_bodovka_sp, aes(Datum, druh)) +
  geom_point(aes(size = n), colour = "blue", fill = "blue", shape = 21, ) +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous(labels = as.Date(data_bodovka_sp$Datum), breaks = data_bodovka_sp$Datum)

#############################################



#tabulka s počty foraging actions
tabulka_behav <- as.data.frame.matrix( table(data_behav$line, data_behav$behav) )

# histogram výšky pozorování foraging actions
ggplot(data = data_sp, aes(x = bird_height)) +
  geom_histogram(binwidth = 1)

# abundancni matice pro linie:

tabulka_abundance_behav <- as.data.frame.matrix( table(data_sp$line, data_sp$druh) )
tabulka_abundance_bodovka <- as.data.frame.matrix(table(data_bodovka$druh, data_bodovka$Datum))

#zbytecne slozity zpusob výpočtu SR pro behav
species_richness <- data_sp %>%  
  group_by(line) %>%        
  summarise(n_distinct(druh)) 
#n=25 

### graf linie - foraging method
library(ggplot2)
sbp1<-ggplot(data_behav, aes(x = line, fill = behav)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title = element_blank())+
  geom_bar(position="fill", color= "black")+
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_brewer(palette = "Spectral")+
  labs(x="linie", y="", title = "Foraging method")

## graf frekvenci vyuziti substratu

sbp2<-ggplot(data_substrate, aes(x = line, fill = substrate)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank())+
  geom_bar(position="fill", color= "black")+
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_brewer(palette = "Spectral")+
  labs(x="line", y="", title = "Foraging substrate")
library(patchwork)
sbp1 + sbp2 ### plot side by side

# preference for foliage density and distance from stem
## foliage density 
sbp3<- data_sp%>% 
  filter(!is.na(dist_stem))%>%
  drop_na(foliage_dens)%>%
  ggplot(aes(x = line, fill = factor(foliage_dens, levels=c("low", "medium", "high")))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank())+
  geom_bar(position="fill", color= "black")+
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_brewer(palette = "Greens")+
  labs(x="linie", y="", title = "Foliage density")

## distance from stem
sbp4<-data_sp%>% 
  drop_na(dist_stem)%>%
  ggplot(aes(x = line, fill = factor(dist_stem, levels=c("edge", "outer", "inner", "stem")))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank())+
  geom_bar(position="fill", color= "black")+
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_brewer(palette = "Oranges")+
  labs(x="linie", y="", title = "Distance from stem")

sbp3 + sbp4 ### plot side by side

#################### Uprava dat pro vypocet index ##############################

behav_substrate <- bind_cols(data_behav,data_substrate)%>%
  select(line...3, druh...1, behav, substrate)%>%
  rename(line = line...3, druh=druh...1) 

behav_substrate_fine <- bind_cols(data_behav,data_substrate, data_substrate_fine)%>%
  select(line...3, druh...1, behav, substrate, substrate_fine,)%>%
  rename(line = line...3, druh=druh...1) 

#vypocet proporci vyuziti substratu a metody
prop_substrate<- behav_substrate %>%
  group_by(druh) %>% 
  count(substrate) %>%
  mutate(prop_substrate = prop.table(n))

prop_method <- behav_substrate %>%
  group_by(druh) %>%
  count(behav) %>%
  mutate(prop_method = prop.table(n))

##kod pro vypocet indexu

index_substrate <- behav_substrate %>% #vypocet indexu B a Ba pro substrat
  group_by(druh) %>% 
  count(substrate) %>%
  mutate(prop_substrate = prop.table(n)) %>%
  mutate(pi2=prop_substrate^2)%>%
  select(druh, pi2) %>%
  group_by(druh)%>%
  summarise(B=1/sum(pi2))%>%
  select(druh, B)%>%
  group_by(druh,B)%>%
  summarise(Ba=1-(B-1)/(6-1))

index_method <- behav_substrate %>% #vypocet indexu B a Ba pro metodu
  group_by(druh) %>% 
  count(behav) %>%
  mutate(prop_behav = prop.table(n)) %>%
  mutate(pi2=prop_behav^2)%>%
  select(druh, pi2) %>%
  group_by(druh)%>%
  summarise(B=1/sum(pi2))%>%
  select(druh, B)%>%
  group_by(druh,B)%>%
  summarise(Ba=1-(B-1)/(7-1))

index<- bind_cols(index_method, index_substrate) %>%  #do jednoho df
  select(druh...1, B...2, Ba...3, B...5, Ba...6) %>%
  rename(druh = druh...1, BM=B...2,BaM=Ba...3, BS=B...5, BaS=Ba...6)

### Scatterplot specializace metoda/substrat
plot(BaS ~ BaM, xlim = c(0.4, 1), ylim = c(0.4, 1), data=index, xlab="Specializace na substrat",ylab="Specializace na metodu", pch=19) 
par(mar = c(5, 5, 5, 5))
abline(lm(BaS ~ BaM, index), lw=1.3)
abline(c(0,1), lty=2, col="red")




################## D I S T A N C N I     M A T I C E ###################

dist_mat<-behav_substrate %>%
  group_by(druh) %>% 
  count(druh,behav,substrate, sort=TRUE)%>% 
  unite(col = "behav_substrate", behav, substrate, sep = "-", remove = TRUE) %>% #distancni matice kde behav-substrate kombinace tvori sloupce
  pivot_wider(names_from="behav_substrate",values_from="n")%>%
  replace(is.na(.), 0)

dist_mat_fine <-behav_substrate_fine %>%
  group_by(druh) %>% 
  count(druh, behav ,substrate, substrate_fine, sort=TRUE)%>% 
  unite(col = "behav_substrate_fine", behav, substrate, substrate_fine, sep = "-", remove = TRUE) %>% #distancni matice kde behav-substrate kombinace tvori sloupce
  pivot_wider(names_from="behav_substrate_fine",values_from="n")%>%
  replace(is.na(.), 0)

#odstraneni druhu s n_actions<5 a n_pozorovani<2
dist_mat_best <- dist_mat[!(dist_mat$druh=="Aegithalos_caudatus"|dist_mat$druh=="Oriolus_oriolus"| dist_mat$druh=="Phoenicorus_phoenicorus"| dist_mat$druh=="Sturnus_vulgaris"| dist_mat$druh=="Sylvia_borin"| dist_mat$druh=="Certhia_brachydactyla"| dist_mat$druh=="Turdus_viscivorus"),] %>%
  remove_rownames%>% 
  column_to_rownames(var="druh")

dist_mat_best_fine <- dist_mat_fine[!(dist_mat_fine$druh=="Aegithalos_caudatus"| dist_mat_fine$druh=="Oriolus_oriolus"| dist_mat_fine$druh=="Phoenicorus_phoenicorus"| dist_mat_fine$druh=="Sturnus_vulgaris"| dist_mat_fine$druh=="Sylvia_borin"| dist_mat_fine$druh=="Certhia_brachydactyla"| dist_mat_fine$druh=="Turdus_viscivorus"),] %>%
  remove_rownames%>% 
  column_to_rownames(var="druh")

library(vegan)

vzdalenost <- vegdist(dist_mat_best, method = "bray",na.rm=TRUE)#bray-curtis distance
vzdalenost_fine <- vegdist(dist_mat_best_fine, method = "bray",na.rm=TRUE)#bray-curtis distance

###################### D E N D R O ###############################
#dendrogram
dendro <- vzdalenost %>%
  hclust (method="ward.D2") %>%
  as.dendrogram()


dendro_fine <- vzdalenost_fine %>%
  hclust (method="ward.D2") %>%
  as.dendrogram()

dendro_phylo_best <- phylo_best %>%
  as.dendrogram()

library(ape)
par(mar=c(3,1,1,12))
plot(dendro,horiz=T, main="dendro behav_substrate", )

plot(as.dendrogram(dendro_fine),horiz=T, main = "dendro behav_substrate_fine")



library(dendextend)
dendro_srovnani <- dendlist(dendro, dendro_fine)

tanglegram(dendro_srovnani, 
           common_subtrees_color_lines = TRUE, highlight_distinct_edges  = FALSE, highlight_branches_lwd=FALSE,
           margin_inner=14,
           lwd=2,
           main_left="dendro",
           main_right="dendro_fine"
)

dendro_phylo <- dendlist(dendro, dendro_phylo_best)
tangle_phylo <- tanglegram(dendro_phylo_best, dendro, 
           common_subtrees_color_lines = TRUE,
           highlight_distinct_edges  = FALSE,
           sort=T,
           highlight_branches_lwd=FALSE,
           margin_inner=14,
           lwd=2,
           main_left="phylogeny",
           main_right="behavior"
)



op = par(bg = "#DDE3CA")
plot(tangle_phylo, col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071",
     col.axis = "#F38630", lwd = 2, lty = 3, sub = '', hang = -1, axes = FALSE)


# vector of colors
labelColors = c("#CDB380", "#036564", "#EB6841", "#D46B37")

# cut dendrogram in 4 clusters
clusMember = cutree(dendro, 4)

# function to get color labels
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}

# using dendrapply
clusDendro = dendrapply(dendro, colLab)

# make plot
plot(clusDendro, main = "Behavior dendro", type = "rectangle", horiz = T)






############################################################################
#vypocty
#procentualni vyuziti substratu
data_substrate %>% 
  count(substrate) %>% 
  mutate(percent = n/sum(n)) %>% 
  select(-n) %>% 
  spread(substrate, percent)

#procentualni vyuziti metody
data_behav %>%
  count(behav) %>% 
  mutate(percent = n/sum(n)) %>% 
  select(-n) %>% 
  spread(behav, percent)

#pocet druhu
data_sp %>% 
  summarise(n_distinct(druh))


#preference pozice na strome
data_sp %>%   
  filter(!is.na(dist_stem))  %>%
  count(dist_stem) %>% 
  mutate(percent = n/sum(n)) %>% 
  select(-n) %>% 
  spread(dist_stem, percent)

#preference miry olisteni
data_sp %>% 
  filter(!is.na(foliage_dens))  %>%
  count(foliage_dens) %>% 
  mutate(percent = n/sum(n)) %>% 
  select(-n) %>% 
  spread(foliage_dens, percent)


#korelace Pearson
cor(index$BaS, index$BaM)

model<-lm(index$BaS ~ index$BaM)
summary (model)


########################################################

druhy<-as.data.frame(unique(data_sp$druh))
library("writexl")
write_xlsx(druhy,"E:/SYNC/DP/druhy.xlsx")

#####################################################
## Fylo dendro

phylo_data <- "E:/SYNC/DP/phylo/output.nex"
phylo_data<-ape::read.nexus(phylo_data)
phylo_best<-phylo_data[[1]]
ape::plot.phylo(phylo_data[[92]])
plot.phylo(phylo_best)

