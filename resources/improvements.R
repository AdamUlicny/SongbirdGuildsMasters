### improvements

labelColors_2 = c("#7A67EE", "#008B45", "#036564","#D46B37" )
labelLegend = c("probers", "leaf", "air", "ground")

dist_phylo_mantel <- mantel(xdis=vzdalenost, ydis=phylo_best, method="pearson", permutations=999)
dist_phylo_mantel

plot(vzdalenost, phylo_best)