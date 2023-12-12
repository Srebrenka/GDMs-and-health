###################
#     PaLS 08    #   
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

# creating the function for  SBM
# and creating datasets

#########################
#                       #
#    Load scripts       #
#                       #
#########################
source("/Users/srebrenkaletina/Downloads/Pals_02_create_nets.R") 

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
library(blockmodels) # sbm 
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


#------------------------------------------------------------------------------
#            FUNCTION FOR PARTITION WITH Stochastic Block Models (SBMs)
#------------------------------------------------------------------------------



net_blocksSBM <- function(net, school) {
  mat <- as.matrix(get.adjacency(net))
  set.seed(782333)
  sbm_out <- BM_bernoulli("SBM", 
                          mat, 
                          verbosity = 3, 
                          plotting = "",
                          exploration_factor = 5) # run a bernoulli block model on the matrix. 
  
  sbm_out$estimate()
  best_fit_assignments <- sbm_out$memberships[[which.min(sbm_out$ICL)]] 
  
  class_assignments <- apply(best_fit_assignments$Z, 1, which.max) 
  Ncom <- length(unique(class_assignments))
  mod <- modularity(net, class_assignments) 
  
  vertex_attr(net, "groups") <- class_assignments
  # level3
  bas <- c(round(school, 0), round(Ncom, 0), round(mod, 3))
  
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
  
  va1<- vertex_attr(net)
  va1df <- data.frame(va1)
  va1df$schid <- school
  
  finalres <- list(basic = bas, data_com = datcom, dfres = va1df)
  return(finalres)
}

#--------------------checking

sbms <- net_blocksSBM(net = palnet[[15]], school = 15)

sbms$basic
sbms$data_com
sbms$dfres
sum(sbms$data_com$size)
#-------------------------------------------------------------------------
#                Function for getting the dataset
#-------------------------------------------------------------------------
level1datSBM<- function(resoff, abb) {
  # SBM for result of SBM blockmodeling
  # resoff - the name of the object where results of net_cda are stored
  # abb - in "", the abbravation of the CDA
  # returns a dataframe with 11 variables (for blocks too?)
  
  # putting level 3 info about net size to level_1
  resoff$dfres$net_size_3 <- sum(resoff$data_com$size)
  # map other level 3 to level 1
  resoff$dfres$school_id_3 <- resoff$basic[1]
  resoff$dfres$N_com_3 <- resoff$basic[2]
  resoff$dfres$Modularity_3 <- resoff$basic[3]
  
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
  colnames(cdadata)[2:11] <- paste(abb, colnames(cdadata)[2:11], sep = "_") # ARGUMENT
  # change "name" to id_1
  names(cdadata)[names(cdadata) == 'name'] <- 'id_1'
  return(cdadata) 
}

#----------------------checking
xyx <- level1datSBM(sbms, abb = "sbm")
dim(xyx)
colnames(xyx)
#-----------------------running over all networks and saving the data
sball<-list() # results are here
for (i in 1:22){
  cda <- net_blocksSBM(palnet[[i]], school=i) 
  alldata <- level1datSBM(resoff = cda, abb = "sbm")
  sball[[length(sball) + 1]] = alldata
  print(i)
}

# some sanity check
length(sball)
dim(sball[[1]])
colnames(sball[[2]])

#names(b2all[[1]])[names(b2all[[1]]) == 'com_perN_2'] <- 'bia_com_perN_2'

findat <- rbindlist(sball)
dim(findat) # 3649
colnames(findat)

write.table(findat, file = "SBM_res.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
write.xlsx(findat, "SBM_res.xlsx")


#------------------------- the end --------------------------------------
