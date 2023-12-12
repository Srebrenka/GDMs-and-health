###################
#     PaLS 06b    #
###################

#rm(list = ls())

#################
#               #
#      Notes    #
#               #
#################

# Srebrenka Letina developed this script


#############
#  Purpose  #
#############

# choosing number of clusters for each network 
# creating and executing the function

#########################
#                       #
#    Load scripts       #
#                       #
#########################
source("02_pals_create_nets.R") 
# loaded via source:

#library(foreign)
library(dplyr)
library(statnet)
library(igraph)
library(intergraph)

#########################
#                       #
#    Load packages      #
#                       #
#########################

library(formattable)
library(openxlsx)
library(readxl)
library(sna)
library(plyr)
library(corrplot)
library(pracma)
library(fpc)
library(memisc)
#library(devtools)
#library(network)


#########################
#                       #
#     Main Script       #
#                       #
#########################



#-----------------------------------------------------------------------------
#          FUNCTION FOR DISCOVERING THE OPTIMAL NUMBER OF CLUSTERS
#------------------------------------------------------------------------------

opt_nof <- function(net, highest_nof) {
  # net - igraph object
  # highest_nof - the highest N of clusters to consider
  mat <- as.matrix(get.adjacency(net)) #arg: net
  profiles <- cbind(mat, t(mat)) 
  t_cor <- cor(t(profiles))
  cor_dist <- round((1-t_cor)/2 , 3)
  cor_dist <- as.dist(cor_dist)  
  
  clustobj <- hclust(cor_dist, method = "average")
  
  #Goodness of fit with the full range of cluster cuts
  nodes_mat <- dim(mat)[1]
  rsq_range <- matrix(0,nodes_mat-1,2)
  
  for (clusn in 2:highest_nof){
    tree_groups<-cutree(clustobj, clusn)
    blk <- blockmodel(mat, tree_groups, plabels = rownames(mat)) 
    dens_net <- blk$block.model
    preimage_net <- matrix(0,nodes_mat, nodes_mat)
    for (i in 1:nodes_mat){
      for (j in 1:nodes_mat){
        i_clus <- tree_groups[i]
        j_clus <- tree_groups[j]
        preimage_net[i,j] <- dens_net[i_clus, j_clus]
      }
    }
    preimage_net <- diag.remove(preimage_net)
    unpacked_net <- cbind(Reshape(mat, nodes_mat^2,1), Reshape(preimage_net, nodes_mat^2,1))
    unpacked_net <- unpacked_net[!is.na(unpacked_net[,2]),]
    rsq_range[clusn-1,1] <- clusn
    rsq_range[clusn-1,2] <- cor(unpacked_net)[1,2]^2
  }
  
  colnames(rsq_range) <- c("N of clusters", "R-squared") 
  rsq_range <-rsq_range[1:highest_nof ,] # get rid of na rows
  
  # Avg Jaccard and Instability
  clusters <- highest_nof
  clus.boot <- clusterboot(cor_dist, 
                           B=1000, # Number of bootstrap resamples
                           clustermethod=disthclustCBI, # for hierarchical clustering of matrix
                           method="average", 
                           k=clusters, 
                           count=FALSE) 
  AvgJaccard <- clus.boot$bootmean
  Instability <- clus.boot$bootbrd/1000
  Clusters <- c(1:clusters)
  Eval <- cbind(Clusters, AvgJaccard, Instability)
  
  # merging it with other indicies
  Eval <- as.data.frame(Eval)
  rsq_range <- as.data.frame(rsq_range)
  Eval$BlockR <- rsq_range$"R-squared"
  Eval$BlockR <- round(Eval$BlockR, 2)
  Eval$AvgJaccard <- round(Eval$AvgJaccard, 3)
  
  # the decision rule
  # take five highest in AJ, then select three with lower Instability 
  # and choose 2 with higher R and select one with smaller number of clusters 
  # if it has R value that is not worse by 3% 
  Eval2 <-Eval[order(Eval$AvgJaccard, decreasing = TRUE),]
  Eval2<- Eval2[1: 10,] # look at this depending on the N of run
  Eval2 <-Eval2[order(Eval2$BlockR, decreasing = TRUE),]
  Eval2<- Eval2[1:6,]
  Eval2 <-Eval2[order(Eval2$Instability, decreasing = FALSE),]
  Eval2<- Eval2[1:2,]
  
  # should I search for smaller clusters?
  # making the rule slightly biased to smaller N of clusters
  Eval2 <-Eval2[order(Eval2$Clusters, decreasing = FALSE),]
  #nofc =ifelse((Eval2$BlockR[2] - Eval2$BlockR[1]) >= 0.03, Eval2$Clusters[2], 
  #             Eval2$Clusters[1] )  # main result
  nofc =ifelse((Eval2$BlockR[2] > Eval2$BlockR[1]), Eval2$Clusters[2], 
               Eval2$Clusters[1] )
  resC = Eval2[Eval2$Clusters == nofc,]
  resC$Clusters # relevant info
  # save the info about values for that cluster
  resCc = c(resC$AvgJaccard, resC$Instability, resC$BlockR)
  the_best = c(max(Eval$AvgJaccard), min(Eval$Instability), max(Eval$BlockR))
  the_best # info of some interest
  endR = list(nofc, resCc, the_best)
  return(endR)
}

#------------------------------- checking
san <- opt_nof(palnet[[7]], highest_nof = 20)
san

#-------------------------------------------
# run it on all network and save results
#------------------------------------------
# i have calculated it with highest_nof = 10
cldat <- list()
for (i in 1:22){
  cmeb <- opt_nof(palnet[[i]], highest_nof = 20)
  cldat[[length(cldat) + 1]] = cmeb
  print(i)
}

length(cldat) # this among other contains the optimal N of clusters per net

for(i in 1:22){
  print(cldat[[i]][1])
}

# save all the info from the cldat list
res0 <- c()# the list of the opt. N of clusters per net

res1<- list() # to save
for (i in 1:22) {
  cc <- cldat[[i]][1]
  res0[i] <- cc
  cc <- c(i, cc)
  res1[[i]] = cc
}
# into dataset 
sq <- matrix(unlist(res1), ncol = 2, byrow = TRUE)
colnames(sq) <- c("net", "optNclusters")

int_in <- list()
for (i in 1:22) {
  cc <- cldat[[i]][2]
  int_in[[length(int_in) +1]] = cc
}

sqb <- matrix(unlist(int_in), ncol = 3, byrow = TRUE)
sqb2 <- as.data.frame(sqb)  
sqb2$ID <- 1:nrow(sqb2)
colnames(sqb2) <- c("AvgJaccard", "Instability", "R2", "net")

best_in <- list()
for (i in 1:22) {
  cc <- cldat[[i]][3]
  best_in[[length(best_in) +1]] = cc
}

sqc <- matrix(unlist(best_in), ncol = 3, byrow = TRUE)
sqc2 <- as.data.frame(sqc)  
sqc2$ID <- 1:nrow(sqc2)
colnames(sqc2) <- c("AvgJaccard_best", "Instability_bes", "R2_best", "net")

### merge these datasets
finr <- merge(sq, sqb2, by="net")
finr <- merge(finr, sqc2, by ="net")

write.table(finr, file = "optNclusters2.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)

write.xlsx(finr, "optNclusters2.xlsx")


