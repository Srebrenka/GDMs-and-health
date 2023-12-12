###################
#     PaLS 09    #  
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

#  creating the function for extracting blocks info
# and creating datasets

#########################
#                       #
#    Load scripts       #
#                       #
#########################


# loaded via source:

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
#library(network)

#########################
#                       #
#     Main Script       #
#                       #
#########################


#-----------------------------------------------------------------------------
#         FUNCTION FOR INDIRECT BLOCKMODELING WITH KNOWN NUMBER OF CLUSTERS
#-----------------------------------------------------------------------------

net_blocks <- function(net, clusn, school) {
  # net - igraph object
  # clusn - N of clusters
  # school - "name" of school
  mat <- as.matrix(get.adjacency(net)) # ARG
  profiles <- cbind(mat, t(mat)) 
  t_cor <- cor(t(profiles))
  cor_dist <- (1-t_cor)/2     #first, we need distances out of correlations
  round(cor_dist, digits = 2) 
  cor_dist <- as.dist(cor_dist)  #set object to be distances
  
  clustobj <- hclust(cor_dist, method = "average")
  
  net_groups<-cutree(clustobj, clusn)
  # level1
  net_grdf <- as.data.frame(net_groups)
  vertex_attr(net, "groups") <- net_grdf$net_groups
  
  # this chunck is not necessary - only for vizualization
  blk <- blockmodel(mat, net_groups, plabels = rownames(mat))
  dens_net <- formattable(blk$block.model, digits=2, format = "f")
  #densities <- diag(dens_net) # extra info to be saved for this type of blocks
  
  # level3
  bas <- c(school, clusn)
  
  # level2
  comdat = list()
  for (k in unique(V(net)$groups)) {
    
    swg2 <- net
    ed1 <- gsize(swg2)
    swga<- igraph::delete.vertices(swg2, V(swg2)$groups != k)
    ed2 <- gsize(swga)
    # density of swga - community !!!!!!
    dencom <- round(edge_density(swga, loops = FALSE), 3)
    sizecom <- gorder(swga)
    swgb<- igraph::delete.vertices(swg2, V(swg2)$groups == k)
    ed3 <- gsize(swgb)
    
    edo <- ed1 - ed3 - ed2# careful, it can be zero
    if (edo==0) {
      rat <- 0
    } else {
      eda <- ed2 + edo
      rat <- round(edo/eda, 2)
    }
    eachcom <- c(k, sizecom, ed2, dencom, rat) 
    comdat[[length(comdat) + 1]] = eachcom
  }
  sq <- matrix(unlist(comdat), ncol = 5, byrow = TRUE)
  datcom <- as.data.frame(sq)
  colnames(datcom) = c("com_name","size", "edges","density","rat")
  datcom$schid <- school
  
  # saving the data
  va1<- vertex_attr(net)
  va1df <- data.frame(va1)
  va1df$schid <- school
  
  finalres <- list(basic = bas, data_com = datcom, dfres = va1df)
  return(finalres)
}

#------------------------------------------------------------------------------
#                function for creating MLM datasets 
#-------------------------------------------------------------------------------
# WITHOUT ADDING ATTRIBUTE DATA!!!!

level1datB<- function(resoff, abb) {
  # B for result of blockmodeling
  # resoff - the name of the object where results of net_cda are stored
  # abb - in "", the abbravation of the CDA
  # returns a dataframe with 11 variables (for blocks too?)
  
  # putting level 3 info about net size to level_1
  resoff$dfres$net_size_3 <- sum(resoff$data_com$size)
  # map other level 3 to level 1
  resoff$dfres$school_id_3 <- round(resoff$basic[1], 0)
  resoff$dfres$N_com_3 <- round(resoff$basic[2], 0)
  #resoff$dfres$Modularity_3 <- resoff$basic[3]
  # creating unique ID per com
  resoff$data_com$uniqueID_2 <- resoff$data_com$schid * 1000 + resoff$data_com$com_name
  # map it to level 1
  resoff$dfres$uniqueID_2 <- mapvalues(resoff$dfres$groups, 
                                       from = resoff$data_com$com_name,
                                       to = resoff$data_com$uniqueID_2)
  # map other level 2 to level 1
  resoff$dfres$com_size_2 <- mapvalues(resoff$dfres$uniqueID_2, 
                                       from = resoff$data_com$uniqueID_2,
                                       to = resoff$data_com$size)
  resoff$dfres$edges_2 <- mapvalues(resoff$dfres$uniqueID_2, 
                                    from = resoff$data_com$uniqueID_2,
                                    to = resoff$data_com$edges)
  resoff$dfres$density_2 <- mapvalues(resoff$dfres$uniqueID_2, 
                                      from = resoff$data_com$uniqueID_2,
                                      to = resoff$data_com$density)
  resoff$dfres$rat_2 <- mapvalues(resoff$dfres$uniqueID_2, 
                                  from = resoff$data_com$uniqueID_2,
                                  to = resoff$data_com$rat)
  # add com% of net size
  resoff$dfres$com_perN_2 <- round((resoff$dfres$com_size_2/resoff$dfres$net_size)*100, 1)
  # get rid of the duplicates and non-useful columns
  resoff$dfres$schid <- NULL
  resoff$dfres$groups <- NULL
  # 11 variables - put it in a special dataframe object, 10 for blocks
  cdadata <- resoff$dfres
  # add prefix but not for "name" variable bc it is ID
  # but name needs to be the same across datasets
  colnames(cdadata)[2:10] <- paste(abb, colnames(cdadata)[2:10], sep = "_") # ARGUMENT
  # change "name" to id_1
  names(cdadata)[names(cdadata) == 'name'] <- 'id_1'
  return(cdadata) 
}

xyx <- level1datB(bre, abb = "bia")
dim(xyx)
#------------------------------------------------------------------------------
#       Run the function with optimal number of clusters
#-------------------------------------------------------------------------------
# open dataset with info on opt. N of clusters

op <- read_excel("~/Downloads/optNclusters2.xlsx", col_names = T)
colnames(op)
op$net # in the right order
nc <- c(op$optNclusters) # we use this info in the function

b2all<-list() # results are here
for (i in 1:22){
  cda <- net_blocks(palnet[[i]], clusn=nc[i], school=i) 
  alldata <- level1datB(resoff = cda, abb = "bia")
  b2all[[length(b2all) + 1]] = alldata
}

# some sanity check

b2all[[1]]$bia_E_B1 <- NULL
names(b2all[[1]])[names(b2all[[1]]) == 'com_perN_2'] <- 'bia_com_perN_2'

findat <- rbindlist(b2all)
dim(findat) # 3649
colnames(findat)

write.table(findat, file = "BIA_res.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
write.xlsx(findat, "BIA_res.xlsx")


#------------------------- the end --------------------------------------
