###################
#     PaLS 05     #
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

# creating the code for finding communities 
# and creating datasets


#########################
#                       #
#    Load scripts       #
#                       #
#########################
source("02_pals_create_nets.R") 

# loaded via source:

library(foreign)
library(dplyr)
library(statnet)
library(igraph)
library(intergraph)

#########################
#                       #
#    Load packages      #
#                       #
#########################
library(plyr)
library(purrr)
library(openxlsx)
library(readxl)


#########################
#                       #
#     Main Script       #
#                       #
#########################



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#                              THE FUNCTION 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
net_cda <- function(net, school, for_dir = "no", seedN, FUN,... ) {
  if (for_dir == "no"){
    # # # # #  Symmetrisation # # # # #
    ex <- net
    
    adj <- as_adjacency_matrix(ex)
    mat <-as.matrix(adj)
    net <- network::network(mat, directed =TRUE)
    
    sw <- sna::symmetrize(net, rule="weak")
    
    dimnames(sw) <- dimnames(mat) # now it has names again
    
    swg <- graph_from_adjacency_matrix(sw, mode="undirected") # names still here!!!
  }
  else { swg <- net
  }
  set.seed(seedN)
  EB1 <-FUN(
    swg #,
    #modularity = TRUE, 
    #membership = TRUE
  )
  ncom <- length(table(EB1$membership)) # 9
  mod <- modularity(swg, EB1$membership) 
  
  vertex_attr(swg, "E_B1") <- EB1$membership
  
  bas <- c(school, ncom, round(mod, 3))
  
  
  comdat = list()
  for (k in unique(V(swg)$E_B1)) {
    swg2 <- swg
    ed1 <- gsize(swg2)
    swga<- igraph::delete.vertices(swg2, V(swg2)$E_B1 != k)
    ed2 <- gsize(swga)
    # density of swga - community !!!!!!
    dencom <- round(edge_density(swga, loops = FALSE), 3)
    sizecom <- gorder(swga)
    swgb<- igraph::delete.vertices(swg2, V(swg2)$E_B1 == k)
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
  va1<- vertex_attr(swg)
  va1df <- data.frame(va1)
  va1df$schid <- school
  
  finalres <- list(basic = bas, data_com = datcom, dfres = va1df)
  
  return(finalres)
}


# test

zet = net_cda(palnet[[1]], school =1, for_dir = "no",  seedN = 199876, 
              FUN=cluster_label_prop) 
zet$basic # net based
zet$data_com # community based
zet$dfres # node based

#------------------------------------------------------------------------------ 
#      function that organizes the output of "net_cda" to level 1 dataset 
#        (for each network)
#------------------------------------------------------------------------------- 


level1dat<- function(resoff, abb) {
  # resoff - the name of the object where results of net_cda are stored
  # abb - in "", the abbravation of the CDA
  # returns a dataframe with 11 variables (for blocks too?)
  
  # putting level 3 info about net size to level_1
  resoff$dfres$net_size_3 <- sum(resoff$data_com$size)
  # map other level 3 to level 1
  resoff$dfres$school_id_3 <- round(resoff$basic[1], 0)
  resoff$dfres$N_com_3 <- round(resoff$basic[2], 0)
  resoff$dfres$Modularity_3 <- resoff$basic[3]
  # creating unique ID per com
  resoff$data_com$uniqueID_2 <- resoff$data_com$schid * 1000 + resoff$data_com$com_name
  # map it to level 1
  resoff$dfres$uniqueID_2 <- mapvalues(resoff$dfres$E_B1, from = resoff$data_com$com_name,
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
  resoff$dfres$E_B1 <- NULL
  # 11 variables - put it in a special dataframe object
  cdadata <- resoff$dfres
  # add prefix but not for "name" variable bc it is ID
  # but name needs to be the same across datasets
  colnames(cdadata)[2:11] <- paste(abb, colnames(cdadata)[2:11], sep = "_") # ARGUMENT
  # change "name" to id_1
  names(cdadata)[names(cdadata) == 'name'] <- 'id_1'
  return(cdadata) 
}



##########################
# label propagation
#########################

lball<-list() # all 22 datasets
for (i in 1:22){
  cda <- net_cda(palnet[[i]], school =i, for_dir = "no", seedN = 199876, 
                 FUN=cluster_label_prop) 
  alldata <- level1dat(resoff = cda, abb = "lpu")
  lball[[length(lball) + 1]] = alldata
  
}



findat <- rbindlist(lball)


write.table(findat, file = "LPU_res.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
write.xlsx(findat, "LPU_res.xlsx")

##########################
# leading_eigen
##########################


leball<-list() # all 22 datasets
for (i in 1:22){
  cda <- net_cda(palnet[[i]], school =i, for_dir = "no", seedN = 199876, 
                 FUN=cluster_leading_eigen) 
  alldata <- level1dat(resoff = cda, abb = "leu")
  leball[[length(leball) + 1]] = alldata
  
}

# some sanity check
class(leball[[1]])
length(leball)
dim(leball[[15]])

findat <- rbindlist(leball)
dim(findat) # 3649
colnames(findat)

write.table(findat, file = "LPU_res.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
write.xlsx(findat, "LPU_res.xlsx")

##########################
# edge betweenness
##########################


eball<-list() # all 22 datasets
for (i in 1:22){
  cda <- net_cda(palnet[[i]], school =i, for_dir = "no", seedN = 199876, 
                 FUN=cluster_edge_betweenness,  modularity = TRUE, 
                 membership = TRUE) 
  alldata <- level1dat(resoff = cda, abb = "eu")
  eball[[length(eball) + 1]] = alldata
  
}

# some sanity check
class(eball[[1]])
length(eball)
dim(eball[[15]])

findat <- rbindlist(eball)
dim(findat) # 3649
colnames(findat)

write.table(findat, file = "EU_res.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
write.xlsx(findat, "EU_res.xlsx")

#--------------------------------------------------------------------------------

#####################################
# other algorithms
#####################################


##########################
# for directed - infomap
##########################

infall2<-list()
for (i in 1:22){
  cda <- net_cda(palnet[[i]], school =i, for_dir = "yes", seedN = 199876, 
                 FUN=cluster_infomap, e.weights = NULL, v.weights = NULL, nb.trials = 10,
                 modularity = TRUE)
  alldata <- level1dat(resoff = cda, abb = "imd")
  infall2[[length(infall2) + 1]] = alldata
  
}
findat <- rbindlist(infall2)

write.table(findat, file = "IMD_res.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
write.xlsx(findat, "IMD_res.xlsx")

##########################
# fast_greedy  - undirected only
##########################
cfg<-list()
for (i in 1:22){
  cda <- net_cda(palnet[[i]], school =i, for_dir = "no",  seedN = 199876, 
                 FUN=cluster_fast_greedy, weights = NULL) 
  alldata <- level1dat(resoff = cda, abb = "cfg")
  cfg[[length(cfg) + 1]] = alldata
  
}

findat <- rbindlist(cfg)
dim(findat)

write.table(findat, file = "CFG_res.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
write.xlsx(findat, "CFG_res.xlsx")

##########################
# louvian - undirected only
##########################

lou<-list()
for (i in 1:22){
  cda <- net_cda(palnet[[i]], school =i, for_dir = "no",  seedN = 199876, 
                 FUN=cluster_louvain, weights = NULL) 
  alldata <- level1dat(resoff = cda, abb = "lo")
  lou[[length(lou) + 1]] = alldata
  
}

findat <- rbindlist(lou)
dim(findat)

write.table(findat, file = "LO_res.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
write.xlsx(findat, "LO_res.xlsx")


##########################
# Walktrap for directed 
##########################

Walk2<-list() # all 22 datasets
for (i in 1:22){
  cda <- net_cda(palnet[[i]], school =i, for_dir = "yes",  seedN = 199876,  
                 FUN=cluster_walktrap, steps = 10) 
  alldata <- level1dat(resoff = cda, abb = "wd")
  Walk2[[length(Walk2) + 1]] = alldata
  
}

findat <- rbindlist(Walk2)
dim(findat)

write.table(findat, file = "WD_res.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
write.xlsx(findat, "WD_res.xlsx")













