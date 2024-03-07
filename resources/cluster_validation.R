install.packages("fpc")
library(fpc)
clusters <- 5
clus.boot <- clusterboot(dist_mat_best, 
                         B=1000, 
                         clustermethod=hclustCBI, 
                         method="ward.D2", 
                         k=clusters, 
                         count=TRUE) 
### výstup
AvgJaccard <- clus.boot$bootmean
Instability <- clus.boot$bootbrd/1000
Clusters <- c(1:clusters)
Eval <- cbind(Clusters, AvgJaccard, Instability)
Eval


clus.boot_fine <- clusterboot(dist_mat_best_fine, 
                              B=1000, 
                              clustermethod=hclustCBI, 
                              method="ward.D2", 
                              k=clusters, 
                              count=TRUE) 

### výstup fine
AvgJaccard_fine <- clus.boot_fine$bootmean
Instability_fine <- clus.boot_fine$bootbrd/1000
Clusters_fine <- c(1:clusters)
Eval_fine <- cbind(Clusters_fine, AvgJaccard_fine, Instability_fine)
Eval_fine

