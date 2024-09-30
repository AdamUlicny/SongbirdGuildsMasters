################################### Guild structure of a passerine assemblage in a Czech lowland deciduous forest ##############################
################## Libraries #################

#packages_to_install <- c("tidyverse", "ape", "vegan", "png","dendextend", "readxl", "patchwork", "treemap", "Rtapas")
#install.packages(packages_to_install)
library(tidyverse)
library(readxl)
library(vegan)
library(dendextend)
library(ape)
library(patchwork)
library(rstudioapi)
library(treemap)
library(ggrepel)
library(Rtapas)

################ Importing data ################

setwd(dirname(getActiveDocumentContext()$path))
data_23 <- read_excel("./data/behav_data_23.xlsx")
data_24 <- read_excel("./data/behav_data_24.xlsx")
data_bodovka <- read_excel("./data/bodovka_data_23.xlsx")
data_lit<-read_excel("./data/literature_data.xlsx")

############# Data wrangling ##################

### behavior, delete empty observations, unite columns into species, pivot to long format
data_behav_23 <- data_23 %>%
  filter(!is.na(behavior_1)) %>%
  unite(col = "druh", genus, species, sep = "_", remove = TRUE) %>%
  pivot_longer(cols = c("behavior_1", "behavior_2", "behavior_3", "behavior_4", "behavior_5"), 
               names_to = "x", values_to = "behav") %>%
  select(druh, behav, line ) %>%
  na.omit()

data_behav_24 <- data_24 %>%
  filter(!is.na(behavior_1)) %>%
  unite(col = "druh", genus, species, sep = "_", remove = TRUE) %>%
  pivot_longer(cols = c("behavior_1", "behavior_2", "behavior_3", "behavior_4", "behavior_5"), 
               names_to = "x", values_to = "behav") %>%
  select(druh, behav, line ) %>%
  na.omit()



### substrate, delete empty observations, unite columns into species, pivot to long format
data_substrate_23 <- data_23 %>%
  filter(!is.na(behavior_1)) %>%
  unite(col = "druh", genus, species, sep = "_", remove = TRUE) %>%
  pivot_longer(cols = c("substrate_main_1", "substrate_main_2", "substrate_main_3", 
                        "substrate_main_4", "substrate_main_5"), names_to = "x", values_to = "substrate") %>%
  select(druh, line, substrate) %>%
  na.omit(F)

data_substrate_24 <- data_24 %>%
  filter(!is.na(behavior_1)) %>%
  unite(col = "druh", genus, species, sep = "_", remove = TRUE) %>%
  pivot_longer(cols = c("substrate_main_1", "substrate_main_2", "substrate_main_3", 
                        "substrate_main_4", "substrate_main_5"), names_to = "x", values_to = "substrate") %>%
  select(druh, line, substrate) %>%
  na.omit(F)

### fine substrate, delete empty observations, unite columns into species, pivot to long format
data_substrate_fine <-data_all %>%
  filter(!is.na(behavior_1)) %>%
  unite(col = "druh", genus, species, sep = "_", remove = TRUE) %>%
  pivot_longer(cols = c("substrate_fine_1", "substrate_fine_2", "substrate_fine_3", 
                        "substrate_fine_4", "substrate_fine_5"), names_to = "x", 
               values_to = "substrate_fine", values_drop_na = T) %>%
  select(druh, line, substrate_fine) 

### literature data

data_lit_behav <- data_lit %>%
  select(species_orig,flycatch,glean,hover_snatch,snatch,pounce, probe)%>%
  pivot_longer(cols=c("flycatch","glean","hover_snatch","snatch","pounce", "probe"), names_to="behav",values_to="x")

data_lit_substrate <- data_lit %>%
  select(species_orig, air, bark, flower, ground, leaf, other)%>%
  pivot_longer(cols=c("air","bark","flower","ground","leaf","other"), names_to="substrate",values_to="x")

### point transect, unite columns, remove non-oscines, count observations
data_bodovka_sp <- data_bodovka %>%
  unite(col = "druh", Genus, Species, sep= "_", remove =TRUE) %>%
  group_by(Datum, druh)%>%
  filter(!(druh=="Columba_oenas"| druh=="Columba_palumbus"| druh=="Buteo_buteo"| druh=="Cuculus_canorus"| druh=="Corvus_corax"| druh=="Dryocopus_martius"| druh=="Dendrocopos_medius"| druh=="Dendrocopos_major"| druh=="Garrulus_glandarius"))%>%
  count()

############### Point transect results #######################

## square plot pro presence/absence z bodovky
ggplot(data_bodovka_sp) +
  aes(x = Datum, y = druh, fill = n) +
  geom_tile(linewidth = 1.2) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous(labels = as.Date(data_bodovka_sp$Datum), breaks = data_bodovka_sp$Datum)

## bubble plot

ggplot(data_bodovka_sp, aes(Datum, druh)) +
  geom_point(aes(size = n), colour = "blue", fill = "blue", shape = 21, ) +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_x_continuous(labels = as.Date(data_bodovka_sp$Datum), breaks = data_bodovka_sp$Datum)

## transformace dat pro porovnani

data_freq_bodovka<- data_bodovka %>%
  unite(col = "druh", Genus, Species, sep= "_", remove =TRUE) %>%
  group_by(druh)%>%
  filter(!(druh=="Columba_oenas"| druh=="Columba_palumbus"| druh=="Buteo_buteo"| druh=="Cuculus_canorus"| druh=="Corvus_corax"| druh=="Dryocopus_martius"| druh=="Dendrocopos_medius"| druh=="Dendrocopos_major"| druh=="Garrulus_glandarius"))%>%
  count()%>%
  summarise(n)%>%
  mutate(prop_bodovka = proportions(n))

data_freq_behav <- data_behav_23%>%
  group_by(druh)%>%
  count()%>%
  summarise(n)%>%
  mutate(prop_behav=proportions(n))

## treemap zastoupení druhu bodovka vs pozorovani

treemap_bodovka<-treemap(data_freq_bodovka,
                         index="druh",
                         vSize="n",
                         type="index",
                         title="freq point transect"
)

treemap_behav<-treemap(data_freq_behav,
                       index="druh",
                       vSize="n",
                       type="index",
                       title="freq obs behav"
)

## frekvenční tabulka

data_freq_compare<-merge(data_freq_bodovka, data_freq_behav, by="druh", all=TRUE)

## proporce druhu na bodovce vs pozorování
ggplot(data_freq_compare) +
  aes(x = prop_behav, y = prop_bodovka) +
  geom_point(
    shape = "asterisk",
    size = 1.45,
    colour = "#0C350C"
  ) +
  labs(
    x = "Proportion of behavioral observations",
    y = "Proportion of observations in point transect"
  ) +
  theme_bw() +
  xlim(0, 0.2) +
  ylim(0, 0.2) + geom_abline()+ geom_text(label=data_freq_compare$druh)

############ Graphs for behavior and substrate #############

#tabulka s počty foraging actions
tabulka_behav <- as.data.frame.matrix( table(data_behav$line, data_behav$behav) )

# histogram výšky pozorování foraging actions
ggplot(data = data_sp, aes(x = bird_height)) +
  geom_histogram(binwidth = 1)

# abundancni matice pro linie:

tabulka_abundance_behav <- as.data.frame.matrix( table(data_sp$line, data_sp$druh) )
tabulka_abundance_bodovka <- as.data.frame.matrix(table(data_bodovka_sp$druh, data_bodovka_sp$Datum))


### graf linie - foraging method
library(ggplot2)
sbp1<-ggplot(data_behav, aes(x = line, fill = behav)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title = element_blank(), axis.text.x=element_blank(), legend.background = element_rect(fill='transparent'),
        axis.ticks.x=element_blank())+
  geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_brewer(palette = "Pastel1")+
  labs(x="", y="", title = "Foraging method")

## graf frekvenci vyuziti substratu

sbp2<-ggplot(data_substrate, aes(x = line, fill = substrate)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank(),axis.text.x=element_blank(),legend.background = element_rect(fill='transparent'), 
        axis.ticks.x=element_blank())+
  geom_bar(position="fill")+
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +scale_fill_manual(
    values = c(air = "#AEEEEE",
               bark = "#F5DEB3",
               ground = "#8B7500",
               leaf = "#556B2F",
               other = "#2F4F4F"))+
  labs(x="", y="", title = "Foraging substrate")

sbp1 + sbp2

## old colors
#values = c(air = "#5AB2A8",
#           bark = "#CFA154",
#           ground = "#543005",
#           leaf = "#008B45",
#           other = "#003D88"))

##scale_fill_manual(
#values = c(flycatch = "#5AB2A8",
#           glean = "#CFA154",
#           hang_glean = "#FFC17A", 
#           hover_snatch = "#003C30",
#           manipulation = "#543005",
#           pounce = "goldenrod",
#           probe = "#8B8B00",
#           snatch =  "#003D88"))

transparent_substrate <- ggplot(data_substrate) +
  aes(x = line, fill = substrate) +
  geom_histogram(bins = 30L) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'), legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(
    values = c(air = "#5AB2A8",
               bark = "#CFA154",
               ground = "#543005",
               leaf = "#003C30",
               other = "#003D88"))

transparent_substrate


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

################# Specialisation index ############

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

index_method_subset <-index_method[!(index_method$druh=="Aegithalos_caudatus"|index_method$druh=="Oriolus_oriolus"| index_method$druh=="Phoenicorus_phoenicorus"| index_method$druh=="Sturnus_vulgaris"| index_method$druh=="Sylvia_borin"| index_method$druh=="Certhia_brachydactyla"| index_method$druh=="Turdus_viscivorus"),]

index_substrate_subset <-index_substrate[!(index_substrate$druh=="Aegithalos_caudatus"|index_substrate$druh=="Oriolus_oriolus"| index_substrate$druh=="Phoenicorus_phoenicorus"| index_substrate$druh=="Sturnus_vulgaris"| index_substrate$druh=="Sylvia_borin"| index_substrate$druh=="Certhia_brachydactyla"| index_substrate$druh=="Turdus_viscivorus"),]

index<- bind_cols(index_method, index_substrate) %>%  #do jednoho df
  select(druh...1, B...2, Ba...3, B...5, Ba...6) %>%
  rename(druh = druh...1, BM=B...2,BaM=Ba...3, BS=B...5, BaS=Ba...6)

index_subset<- bind_cols(index_method_subset, index_substrate_subset) %>%  #do jednoho df
  select(druh...1, B...2, Ba...3, B...5, Ba...6) %>%
  rename(druh = druh...1, BM=B...2,BaM=Ba...3, BS=B...5, BaS=Ba...6)



### Scatterplot specializace metoda/substrat
par(mar = c(5, 5, 5, 5))
plot(BaS ~ BaM, xlim = c(0.4, 1), ylim = c(0.4, 1), data=index, xlab="Specializace na substrat",ylab="Specializace na metodu", pch=19) 
abline(lm(BaS ~ BaM, index), lw=1.3)
abline(c(0,1), lty=2, col="red")

par(mar = c(5, 5, 5, 5))
plot(BaS ~ BaM, xlim = c(0.4, 1), ylim = c(0.4, 1), ylab="Substrate specialization",xlab="Method specialization", data=index_subset, cex.lab=2, pch=19)
abline(lm(BaS ~ BaM, index_subset), lw=1.3)
abline(c(0,1), lty=2, col="red")

summary(lm(BaS ~ BaM, index_subset))


################# Disimilarity matrices ###############################

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


#################### Literature ###############
dist_mat_lit_behav <- data_lit_behav%>%
  mutate(across(where(is.numeric), round, 0))%>%
  replace(is.na(.), 0)%>%
  group_by(species_orig, behav)%>%
  summarise(x=sum(x))%>%
  pivot_wider(names_from="behav", values_from="x")%>%
  remove_rownames%>% 
  column_to_rownames(var="species_orig")


dist_mat_lit_substrate <- data_lit_substrate%>%
  mutate(across(where(is.numeric), round, 0))%>%
  replace(is.na(.), 0)%>%
  group_by(species_orig, substrate)%>%
  summarise(x=sum(x))%>%
  pivot_wider(names_from="substrate", values_from="x")%>%
  remove_rownames%>% 
  column_to_rownames(var="species_orig")


#odstraneni druhu s n_actions<5 a n_pozorovani<2
dist_mat_best <- dist_mat[!(dist_mat$druh=="Aegithalos_caudatus"|dist_mat$druh=="Oriolus_oriolus"| dist_mat$druh=="Phoenicorus_phoenicorus"| dist_mat$druh=="Sturnus_vulgaris"| dist_mat$druh=="Sylvia_borin"| dist_mat$druh=="Certhia_brachydactyla"| dist_mat$druh=="Turdus_viscivorus"),] %>%
  remove_rownames%>% 
  column_to_rownames(var="druh")

dist_mat_best_fine <- dist_mat_fine[!(dist_mat_fine$druh=="Aegithalos_caudatus"| dist_mat_fine$druh=="Oriolus_oriolus"| dist_mat_fine$druh=="Phoenicorus_phoenicorus"| dist_mat_fine$druh=="Sturnus_vulgaris"| dist_mat_fine$druh=="Sylvia_borin"| dist_mat_fine$druh=="Certhia_brachydactyla"| dist_mat_fine$druh=="Turdus_viscivorus"),] %>%
  remove_rownames%>% 
  column_to_rownames(var="druh")

vzdalenost <- vegdist(dist_mat_best, method = "bray",na.rm=TRUE)#bray-curtis distance
vzdalenost <- as.matrix(vzdalenost)
vzdalenost <- as.dist(vzdalenost[order(rownames(vzdalenost)),order(colnames(vzdalenost))])

vzdalenost_fine <- vegdist(dist_mat_best_fine, method = "bray",na.rm=TRUE)#bray-curtis distance

dist_lit_behav <- vegdist(dist_mat_lit_behav, method = "bray",na.rm=TRUE)
dist_lit_behav <- as.matrix(dist_lit_behav)
dist_lit_behav <- as.dist(dist_lit_behav[order(rownames(dist_lit_behav)),order(colnames(dist_lit_behav))])

dist_lit_substrate <-vegdist(dist_mat_lit_substrate, method = "bray",na.rm=TRUE)
dist_lit_substrate <- as.matrix(dist_lit_substrate)
dist_lit_substrate <- as.dist(dist_lit_substrate[order(rownames(dist_lit_substrate)),order(colnames(dist_lit_substrate))])

#################### Morpho data ###############################
data_morfo <- read_xlsx("./data/morfo_data.xlsx")

### drop unused columns
data_morfo_matrix <- data_morfo %>%
  select(species_birdlife,Beak.Length_Culmen, Beak.Width, 
         Beak.Depth,Tarsus.Length, Wing.Length, `Hand-Wing.Index`, Tail.Length)

### log transform data
data_morfo_matrix_log <- data_morfo_matrix%>%
  mutate_at(2:8, log)%>%
  remove_rownames%>% 
  column_to_rownames(var="species_birdlife")


### distance matrix
dist_morfo<-vegdist(data_morfo_matrix_log, method="euclidean")
## alphabetical sorting
dist_morfo <- as.matrix(dist_morfo)
dist_morfo <- as.dist(dist_morfo[order(rownames(dist_morfo)),order(colnames(dist_morfo))])

## cluster 
dendro_morfo<-dist_morfo %>%
  hclust(method="ward.D2")%>%
  as.dendrogram()

#################### Phylogenetic data ############################

#### Exporting species names

#druhy<-as.data.frame(unique(data_sp$druh))
#library("writexl")
#write_xlsx(druhy,"druhy.xlsx")

### Importing trees

phylo_data <- "./resources/output.nex"
phylo_data<-ape::read.nexus(phylo_data)

### Selecting best tree
ape::plot.phylo(phylo_data[[92]])
phylo_best<-phylo_data[[1]]

### alphabetic sorting
phylo_best <- dist(cophenetic(phylo_best))
phylo_best <- as.matrix(phylo_best)
phylo_best <- as.dist(phylo_best[order(rownames(phylo_best)),order(colnames(phylo_best))])


################ Mantel test #####################################
simple.results.mantel <- mantel(xdis = vzdalenost, ydis = dist_morfo, method = "pearson", permutations = 999)
simple.results.mantel

plot(vzdalenost, dist_morfo)

dist_phylo_mantel <- mantel(xdis=vzdalenost, ydis=phylo_best, method="pearson", permutations=999)
dist_phylo_mantel

plot(vzdalenost, phylo_best)

################ Dendrograms #####################################

###### behavior-substrate dendro
dendro <- vzdalenost %>%
  hclust (method="ward.D2") %>%
  as.dendrogram()

###### behavior-fine substrate dendro
dendro_fine <- vzdalenost_fine %>%
  hclust (method="ward.D2") %>%
  as.dendrogram()

#### literature dendro behav
dendro_lit_behav<-dist_lit_behav %>%
  hclust(method="ward.D2")%>%
  as.dendrogram()


dendro_lit_substrate<-dist_lit_substrate %>%
  hclust(method="ward.D2")%>%
  as.dendrogram()

##### phylogeny dendro
detach(package:ape, unload=T)
dendro_phylo_best <- phylo_data[[1]] %>%
  as.dendrogram()%>%
  dendextend::rotate(c(1:6,7:18))

dendro_phylo_best<-rotate(dendro_phylo_best, c(1:6,7:18))
dendro_phylo_best

par(mar=c(3,1,1,12))
plot(dendro,horiz=T, main="dendro behav_substrate", )

plot(as.dendrogram(dendro_fine),horiz=T, main = "dendro behav_substrate_fine")

plot(as.dendrogram(dendro_phylo_best),horiz=T, main = "phylogeny", rotate(as.character(1:7)))

par(mar=c(2,1,1,12))
plot(as.dendrogram(dendro_lit_behav),horiz=T, main = "literature_behav")

plot(as.dendrogram(dendro_lit_substrate),horiz=T, main = "literature_substrate")



############################## tanglegram ######################################

tanglegram(dendro, dendro_fine, 
           common_subtrees_color_lines = TRUE, highlight_distinct_edges  = FALSE, highlight_branches_lwd=FALSE,
           margin_inner=14,
           lwd=2,
           sort=T,
           main_left="dendro",
           main_right="dendro_fine")



tangle_phylo <- tanglegram(dendro, dendro_phylo_best, 
                           common_subtrees_color_lines = TRUE,
                           highlight_distinct_edges  = FALSE,
                           sort=T,
                           highlight_branches_lwd=FALSE,
                           margin_inner=14,
                           lwd=2,
                           main_left="behavior",
                           main_right="phylogeny",
                           hang=F,)

rotate_tangle <- untangle_step_rotate_1side(
  dendro_phylo_best,
  dendro,
  L = 1.5,
  direction = c("forward", "backward"),
  k_seq = NULL,
  dend_heights_per_k,
  leaves_matching_method = c("labels", "order"))
plot(rotate_tangle)


tangle_morfo_behav <- tanglegram(dendro_morfo, dendro, 
                                 common_subtrees_color_lines = TRUE,
                                 highlight_distinct_edges  = FALSE,
                                 sort=T,
                                 highlight_branches_lwd=FALSE,
                                 margin_inner=14,
                                 lwd=2,
                                 main_left="morphology",
                                 main_right="behavior",
                                 hang=F,)

tangle_morfo_phylo <- tanglegram(dendro_morfo, dendro_phylo_best, 
                                 common_subtrees_color_lines = TRUE,
                                 highlight_distinct_edges  = FALSE,
                                 sort=T,
                                 highlight_branches_lwd=FALSE,
                                 margin_inner=14,
                                 lwd=2,
                                 main_left="morphology",
                                 main_right="phylogeny",
                                 hang=F,)


######################### aesthetics ####################################################################
op = par(bg = "white")
# vector of colors
labelColors_2 = c("#7A67EE", "#008B45", "#036564","#D46B37" )
labelLegend = c("probers", "leaf", "air", "ground")

# nařezat dendro do 4 gild
clusMember = cutree(dendro, 4)

# funkce pro obarvení gild
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors_2[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}

# dendrapply aplikace funkce pro guildy
dendro_pretty = dendrapply(dendro, colLab)
dendro_fine_pretty = dendrapply(dendro_fine, colLab)

### plot gildy barvy legenda
par(mar=c(5,1,1,12))
plot(dendro_pretty, main = "Foraging guilds", type = "rectangle", horiz = T)
legend("topleft", 
       legend = labelLegend, 
       col = labelColors_2, 
       pch = c(20,20,20,20), bty = "y",  pt.cex = 1.5, cex = 0.8 , 
       text.col = "black", horiz = FALSE,
       title="Guilds",
       inset = c(0, 0.1))

### kombo gildy_phylo estetika

tanglegram(dendro_pretty, dendro_phylo_best, 
           common_subtrees_color_lines = TRUE,
           highlight_distinct_edges  = FALSE,
           sort=T,
           highlight_branches_lwd=FALSE,
           margin_inner=14,
           margin_outer=7,
           lwd=3,
           main_left="behavior",
           main_right="phylogeny",
           hang=F,)


legend("topleft", 
       legend = labelLegend, 
       col = labelColors_2, 
       pch = c(20,20,20,20), bty = "n",  pt.cex = 1.5, cex = 0.8 , 
       text.col = "black", horiz = FALSE, inset = c(0, 0.1))






######################################################################################################################
################### Calculations #######################

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
cor(index_subset$BaS, index_subset$BaM, method="pearson")

model<-lm(index$BaS ~ index$BaM)
summary (model)

summary(lm(index_subset$BaS ~ index_subset$BaM))