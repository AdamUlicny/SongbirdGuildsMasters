#### Parulidae project, initiated 22.2.2016, V. Remes ####

source("/Users/Lada/Documents/stats, mat & program/R/R_projekty/prace se soubory/funkce_soubory.R")
setwd("~/Documents/clanky_a_projekty/Parulidae/PROJEKT_Parulidae")
require(ape)
require(geiger)

##### extract names from BLI shapefiles -----
Parul.jm.2011 <- SpeciesNamesCodes( GetFileNames("data/Parulidae_2011") )
Parul.jm.2014 <- SpeciesNamesCodes( GetFileNames("data/Parulidae_2014") )
write.csv(Parul.jm.2014, file="BLI_Parulidae.csv")

##### load phylogeny -----
# phylogeny for Parulidae from Barker et al. 2015 Auk, http://dx.doi.org/10.5061/dryad.pb787
load("data/phylogenies/emberizoid.supertrees.RData")  # all trees
load("data/phylogenies/Parulidae.RData")  # Parulidae all trees
trees <- read.nexus("data/phylogenies/parulids_all_trees.tre")  # Parulidae all trees in nexus format
mcct <- read.nexus("data/phylogenies/mcc_parulids.tre")  # Parulidae mcc tree
branching.times(mcct)[1]  # root height
mcct1 <- rescale(mcct, model="depth", 1)
branching.times(mcct1)[1]  # root height=1
plot(mcct$edge.length~mcct1$edge.length)  # OK
# mcc tree is virtually identical to the ML tree in Fig. 5 in MPE Lovette et al. 2010

##### prepare names of species -----
# t1 <- parulids[[1]]$tip.label  # nahodne vybrany strom
# t2 <- parulids[[100]]$tip.label  # nahodne vybrany strom
# identical(t1,t2)  # yes
require(readxl)
names <- read_excel(path="data/Families_Passerines.xlsx", sheet = "Parulidae")  # path: GACR2016-veda-Kuba prace
nam.lovette <- paste(names$Lovette2010_Genus[!is.na(names$Lovette2010_Genus)], names$Lovette2010_Species[!is.na(names$Lovette2010_Species)], sep="_")
names(nam.lovette) <- nam.lovette
nam.barker <- paste(names$Barker2015_Dryad_Genus[!is.na(names$Barker2015_Dryad_Genus)], names$Barker2015_Dryad_Species[!is.na(names$Barker2015_Dryad_Species)], sep="_")
names(nam.barker) <- nam.barker

# prune names to include only spp on phylogeny:
pruned.names <- treedata(phy = mcct, data = nam.barker)$data
name.check(phy=mcct, data=pruned.names)  # making sure that pruned.names are OK

# attach clade names:
names$Barker2015_Dryad <- paste(names$Barker2015_Dryad_Genus, names$Barker2015_Dryad_Species, sep="_")
pruned.names <- as.data.frame(pruned.names, stringsAsFactors=F)
names(pruned.names)[1] <- "Barker2015_Dryad"
pruned.names <- merge(x=pruned.names, y=names[, c(10,14)], by="Barker2015_Dryad", all.x=T, all.y=F)

# add clades data on parulid morphology with spp names
d <- merge(x=Parulidae.spp, y=pruned.names, by.x="phyloName", by.y="Barker2015_Dryad", all=T)
d <- na.omit(d)  # remove Basileuterus_leucophrys = I do not have data for this species
rownames(d) <- d$phyloName

# prune mcct1 tree = remove Basileuterus_leucophrys
mcct1pr <- treedata(phy=mcct1, data=d)$phy
name.check(phy=mcct1pr, data=d)  # overeni, ze data a phylogeny souhlasi = OK

##### plotting tree with traits -----
# adephylo
require(adephylo)
require(phylobase)
require(phylosignal)
X <- phylo4d(x=mcct1pr, tip.data=d)  # morphometric traits
X <- phylo4d(x=mcct1pr, tip.data=PCA_scores_log)  # pca or ppca axes
quartz(width = 7, height = 13)
par(mar=c(2,2,2,6))
table.phylo4d( X, treetype="phylogram", symbol="circles", ratio.tree=0.45, center=T, scale=T, legend=T, grid=T, box=F, cex.symbol=0.9, cex.label=0.6, cex.legend=0.8, var.label=c("", names(d)[2:7]),  )
table.phylo4d( X, treetype="phylogram", symbol="circles", ratio.tree=0.45, center=T, scale=T, legend=T, grid=T, box=F, cex.symbol=0.9, cex.label=0.6, cex.legend=0.8, var.label=colnames(scores.spp)  )
# phylosignal: wrappers for multiplot.phylo4d {phylosignal}, uses phylo4d object
quartz(width = 15, height = 15)
barplot(X, tree.ladderize=T, center=T, scale=T, tree.type="phylogram", trait=c("b_length","b_width","b_height","tarsus","tail_length","wing_length"), bar.col=cols)  # prepare cols below
dotplot(X, trait=c("b_length","b_width","b_height"))
gridplot(X, trait=c("b_length","b_width","b_height"))

# colors for phytools
tree <- ladderize(mcct1pr)
plotTree(tree, node.numbers=T)
tree <- paintSubTree(tree, node=203, state="2")  # BasalClade
tree <- paintSubTree(tree, node=113, state="3")  # Basileuterus
tree <- paintSubTree(tree, node=121, state="4")  # Cardellina
tree <- paintSubTree(tree, node=183, state="5")  # Geothlypis
tree <- paintSubTree(tree, node=125, state="6")  # Myioborus
tree <- paintSubTree(tree, node=136, state="7")  # Myiothlypis
tree <- paintSubTree(tree, node=196, state="8")  # Oreothlypis
tree <- paintSubTree(tree, node=149, state="9")  # Setophaga
cols <- c("black","deep skyblue3","deeppink2","dark magenta","seagreen","slate blue","sienna","gold","goldenrod"); names(cols) <- 1:9
plotSimmap(tree, color=cols, pts=F, lwd=3, node.numbers=F, fsize=0.7)

# phytools
require(phytools)
dt <- d[tree$tip.label, ]
v <- dt$tail_length
names(v) <- rownames(dt)
phylomorphospace(tree, X=dt[, c(5,7)], colors=cols, lwd=2.5, label="off", xlab="tarsus length", ylab="tail length")
phenogram(tree, x=v, fsize = 0.8, spread.labels=F, colors="darkgray")
phenogram(tree, x=v, fsize = 0.8, spread.labels=T, colors=cols, ylab="tail length")
fancyTree(tree, type="phenogram95", x=v)

##### plotting tree with clades -----
require(ape)
# colored tip labels:
cols <- c(rep("black",105))  # prepare colors for tip labels
cl <- d$Lovette2010_Clade  # based on clades in Lovette et al. 2010
names(cl) <- d$phyloName
cl <- cl[mcct1pr$tip.label]
cols[which(cl=="BasalClade")] <- "deep skyblue3"
cols[which(cl=="Basileuterus")] <- "deeppink2"
cols[which(cl=="Cardellina")] <- "dark magenta"
cols[which(cl=="Geothlypis")] <- "seagreen"
cols[which(cl=="Myioborus")] <- "slate blue"
cols[which(cl=="Myiothlypis")] <- "sienna"
cols[which(cl=="Oreothlypis")] <- "gold"
cols[which(cl=="Setophaga")] <- "goldenrod"

# vertical bars with clade names, can be colored:
quartz(width = 7, height = 13)
plot(ladderize(mcct1pr), edge.col="black", tip.col=cols, cex=0.7, x.lim=c(0, 1.4), no.margin=T, edge.width=2, label.offset=0.004)  # tip color by clade or uniformly black
bar <- 1.34
w <- 3
tx <- 1.37
c1 <- rep("black", 8)
c2 <- rep("black", 8)
c1 <- c2 <- c("slate blue","dark magenta","deeppink2","sienna","goldenrod","seagreen","gold","deep skyblue3")  # prepare colors for clade bars
segments(bar, 1, bar, 12, lwd=w, col=c1[1])
text(tx, 7, "Myioborus", srt=270, col=c2[1])
segments(bar, 13, bar, 17, lwd=w, col=c1[2])
text(tx+0.04, 15, "Cardellina", srt=270, col=c2[2])
segments(bar, 18, bar, 25, lwd=w, col=c1[3])
text(tx, 21.5, "Basileuterus", srt=270, col=c2[3])
segments(bar, 26, bar, 39, lwd=w, col=c1[4])
text(tx, 32.5, "Myiothlypis", srt=270, col=c2[4])
segments(bar, 40, bar, 73, lwd=w, col=c1[5])
text(tx, 57, "Setophaga", srt=270, col=c2[5])
segments(bar, 74, bar, 87, lwd=w, col=c1[6])
text(tx, 81, "Geothlypis", srt=270, col=c2[6])
segments(bar, 88, bar, 95, lwd=w, col=c1[7])
text(tx, 91.5, "Oreothlypis", srt=270, col=c2[7])
segments(bar, 96, bar, 103, lwd=w, col=c1[8])
text(tx, 99.5, "basal cl.", srt=270, col=c2[8])

# colored edges based on clades:
# get edges for clades
mcct1pr <- ladderize(mcct1pr)
cols <- rep("black", Nedge(mcct1pr))
cl <- d$Lovette2010_Clade  # based on clades in Lovette et al. 2010
names(cl) <- d$phyloName
cl <- cl[mcct1pr$tip.label]
Basal <- which.edge( mcct1pr, which(cl=="BasalClade") )
Basileuterus <- which.edge( mcct1pr, which(cl=="Basileuterus") )
Cardellina <- which.edge( mcct1pr, which(cl=="Cardellina") )
Geothlypis <- which.edge( mcct1pr, which(cl=="Geothlypis") )
Myioborus <- which.edge( mcct1pr, which(cl=="Myioborus") )
Myiothlypis <- which.edge( mcct1pr, which(cl=="Myiothlypis") )
Oreothlypis <- which.edge( mcct1pr, which(cl=="Oreothlypis") )
Setophaga <- which.edge( mcct1pr, which(cl=="Setophaga") )

cols[Basal] <- "deep skyblue3"
cols[Basileuterus] <- "deeppink2"
cols[Cardellina] <- "dark magenta"
cols[Geothlypis] <- "seagreen"
cols[Myioborus] <- "slate blue"
cols[Myiothlypis] <- "sienna"
cols[Oreothlypis] <- "gold"
cols[Setophaga] <- "goldenrod"

quartz(width = 7, height = 13)
plot(mcct1pr, edge.col=cols, tip.col="black", cex=0.7, x.lim=c(0, 1.4), no.margin=T, edge.width=3, label.offset=0.004)  # tip color by clade or uniformly black
# zatim nejde kombinovat s vertikalnimi carami, nesedi nektere klady asi kvuli ladderizaci

##### diversification analyses -----
require(ape)
require(laser)
require(phytools)

# ltt plots:
mcct <- read.nexus("/Users/Lada/Documents/clanky_a_projekty/Parulidae/PROJEKT_Parulidae/data/phylogenies/mcc_parulids.tre")  # Parulidae mcc tree
ltt.plot(mcct, log="y", ylab="Log number of lineages", xlab="Relative time")  # for one tree, in ape

btimes <- branching.times(mcct)  # ape
plotLtt(btimes)  # laser

trees <- read.nexus("/Users/Lada/Documents/clanky_a_projekty/Parulidae/PROJEKT_Parulidae/data/phylogenies/parulids_all_trees.tre") 
ltt95(trees[1:100])  # for a set of trees, in phytools
ltt95(trees[1:100], xaxis="flipped", log=TRUE)
pbtrees <- pbtree(n=100, scale=50, nsim=100)  # pure birth trees
ltt95(pbtrees, xaxis="flipped", log=TRUE)

# gamma statistic:
is.ultrametric(mcct)  # tree must be ultrametric
gammaStat(mcct)  # under constant divers, gamma~norm(0,1)
2*(1 - pnorm(abs(gammaStat(mcct))))  # 2-tailed test
1 - pnorm(abs(gammaStat(mcct)))  # 1-tailed test
require(laser)
btimes2 <- getBtimes(file="/Users/Lada/Documents/clanky_a_projekty/Parulidae/PROJEKT_Parulidae/data/phylogenies/mcc_parulids.txt", string=NULL)  # vector with branching times
gamStat(x=btimes2, return.list=TRUE)  # x=numeric vector of branching times
gamma.test <- mccrTest(CladeSize=106, NumberMissing=1, NumberOfReps=5000, ObservedGamma=NULL, fname=NULL)  # list with:
# null.gamma: The null distribution of gamma. You can plot a histogram or otherwise inspect these values...
# critical.value: The 0.05 percentile of the null distribution. This is the value corresponding to alpha = 0.05
# p.value: The actual p-value, only returned if ObservedGamma is supplied by user
hist(gamma.test$null.gamma)









