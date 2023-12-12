###################
#     PaLS 06a    #
###################
setwd("/Users/srebrenkaletina/Downloads")
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

#  choosing the clustering method 

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
#library(devtools)
#library(network)


#########################
#                       #
#     Main Script       #
#                       #
#########################

# size of nets
for (i in 1:22){
  print(gorder(palnet[[i]]))
}
#------------------------------------------------------------------------------
#                    INDIRECT APPROACH TO BLOCKMODELING
#------------------------------------------------------------------------------


###Profile similarity

# We are comparing the relational prfiles of nodes. 
# First we need to generate data that records out- and in-choices.

#get the adjacency matrix


#This is a directed graph, so we need to take the transpose also into account 
# (out and in-profiles)

# We need to assemble profiles into one matrix

profiles <- cbind(mat, t(mat))  #rows will show outging and then incoming ties - the profile
profiles


#We can compare profiles by calculating their vector distances (or similarities)

t_cor <- cor(t(profiles)) #cor is looking down the columns, so we transpose our profiles

round(diag.remove(t_cor), digits = 3)


#The correlation matrix with 'corrplot' See: [https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html]


###Clustering basics: agglomerative hierarchical clustering

#The first step in each iteration is to find the smallest distance, and merge those two items. 
#Methods differ about how the merged entity's new distances are calculated.
#The 'hclust' function computes all mergers and documents that on a dendrogram.

# Plot a 'hclust' object to see the dendrogram.


###Choice of clustering method via goodess of fit

#Cophenetic distance: the distance matrix computed based on dendrogram hierarchy traversal.

# We can compare the cophenetic distances to the original distances.


#We can calculate a goodness of fit (proportion of variance explained).

#How would a different method perform? Single link, for example.

#Fit for complete link


#----------------------------------------------------------------------------
#             choosing clustering method for all networks
#----------------------------------------------------------------------------

# go over each network profile with each clustering method

pernet=list()
for (k in 1:22){
  mat <- as.matrix(get.adjacency(palnet[[k]]))
  net_t <- network(mat, directed=TRUE)
  
  profiles <- cbind(mat, t(mat)) 
  t_cor <- cor(t(profiles))
  cor_dist <- (1-t_cor)/2     #first, we need distances out of correlations
  cor_dist <- as.dist(cor_dist)  #set object to be distances
  round(cor_dist, digits = 2) 
  
  melis <- list("complete", "single", "average", "median", "centroid",
                "ward.D", "ward.D2", "mcquitty")
  
  twores = list()
  for (i in 1:length(melis)){
    clust_example <- hclust(
      cor_dist,    #the distance matrix we cluster
      method = melis[[i]]   #the method of clustering
    )
    clust_example_coph <- as.matrix(cophenetic(clust_example))
    
    unpacked_dist <- cbind(Reshape(mat, dim(mat)[1]^2,1), Reshape(clust_example_coph, dim(clust_example_coph)[1]^2,1))
    example_rsq <- cor(unpacked_dist)[1,2]^2
    
    twores[[length(twores) + 1]] = example_rsq
  }
  pernet[[length(pernet) + 1]] = twores
}

sq <- matrix(unlist(pernet), ncol = 8, byrow = TRUE)
datcom <- as.data.frame(sq)

pernet[[1]]

colnames(datcom) = c("complete", "single", "average", "median", "centroid",
                     "ward.D", "ward.D2", "mcquitty")
datcom <- dplyr::mutate(datcom, ID = row_number())
datcom <- round(datcom, 3)
View(datcom)

write.xlsx(datcom, "hclust_method_reg.xlsx")

# we choose average



