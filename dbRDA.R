# dbRDA for guilds based on morphology-behavior-phylogeny
library(vegan)
library(ape)

phylo_PC<-cophenetic.phylo(phylo_meta)
phylo_PC<-pcoa(phylo_PC)
explained <- phylo_PC$values$Relative_eig
cumulative<-cumsum(explained)
which(cumulative >= 0.80)[1]
bstick <- bstick(n = length(explained))
plot(explained, type = "b", ylim = c(0, max(explained)), main = "Broken Stick vs Observed Eigenvalues")
lines(bstick, col = "red", type = "b")
which(explained > bstick) 

phylo_PC <- phylo_PC$vectors[, 1:11] 
predictors <- cbind(morphology_Global_matrix, matrix_meta_full, phylo_PC)


dbrda_model<-capscale(dist_Bray_Global ~ ., data = predictors)
anova(dbrda_model)  # Tests overall significance
anova(dbrda_model, by = "terms")  # Tests individual variables


varpart(dist_Bray_Global, ~ morphology_Global_matrix, ~ phylo_PC, ~ matrix_meta_full)
