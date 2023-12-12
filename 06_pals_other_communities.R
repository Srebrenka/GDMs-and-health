###################
#     PaLS 06     #
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

# creating the code for finding communities for additional GDMs
# WITH SPECIFIC ARGUMENTS
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
#               MY FUNCTION 2.0 - leiden : 
# objective_function = "modularity", resolution_parameter = 2
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# GENDER INFO MUST BE READ WITH NETWORK DATA
net_cda_le2 <- function(net, school, for_dir = "no", seedN, FUN,... ) {
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
    swg, objective_function = "modularity", resolution_parameter = 2
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

#----------------------------------------------------------------------

ball<-list() # all 22 datasets
for (i in 1:22){
  cda <- net_cda_le1(palnet[[i]], school =i, for_dir = "no", seedN = 199876, 
                     FUN=cluster_leiden) 
  alldata <- level1dat(resoff = cda, abb = "le1")
  ball[[length(ball) + 1]] = alldata
  
}

findat <- rbindlist(ball)
write.table(findat, file = "LE1_res.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
write.xlsx(findat, "LE1_res.xlsx")





###############################################################################
#                           Clique Percolation
###############################################################################
# DIRECTED AND FOR CLIQUES k = 3 
library(CliquePercolation)
# source: https://cran.r-project.org/web/packages/CliquePercolation/vignettes/CliquePercolation.html
# https://cran.r-project.org/web/packages/CliquePercolation/CliquePercolation.pdf
# package not working in Mac :(

##################################################################################
#                               CPM FUNCTION
#################################################################################

# I modified it
source("cpAlgorithm2.R")
# arguments
# undirected network, school

net_pcm <- function(objgr, school){
  # # # # #  Symmetrisation # # # # #
  net <- objgr
  ex <- net
  adj <- as_adjacency_matrix(ex)
  mat <-as.matrix(adj)
  net <- network::network(mat, directed =TRUE)
  sw <- sna::symmetrize(net, rule="weak")
  dimnames(sw) <- dimnames(mat) # now it has names again
  swg <- graph_from_adjacency_matrix(sw, mode="undirected") 
  
  cp <- cpAlgorithm2(sw, 3) # we look at cliques of 3
  
  nccp <- length(cp$list.of.communities.labels)
  nodenames <- rownames(sw)
  memdf <- as.data.frame(nodenames)
  
  fc <- list()
  for(i in 1:nccp){
    fc1 <- cp$list.of.communities.labels[[i]]
    fc[[i]] <- fc1
  }
  for(i in 1:length(cp$isolated.nodes.labels)){
    fc2 <- cp$isolated.nodes.labels[i]
    fc[[nccp + i]] <- fc2
  }
  
  realfc<- length(fc)
  
  ec <- list()
  for(i in 1:realfc){
    kn <- ifelse(memdf$nodenames %in% fc[[i]], i, NA)
    ec[[i]] <- kn
  }
  
  see <- do.call(cbind, ec) # membership matrix
  
  rownames(see) <- memdf$nodenames
  # assign just one community membership
  c1 <- c()
  for(i in 1:nrow(see)){
    par <- see[i,]
    if (sum(!is.na(par)) == 1){
      c1[i] <- sum(par, na.rm = T)
    }
    if(sum(!is.na(par)) == 2){ # in two communities
      par<-c(par[!is.na(par)])
      
      com1 <- par[1]
      com2 <- par[2]
      
      see2 <- see[see[,com1] == com1, ]
      com1nodes <- c(rownames(see2))
      see3 <- see[see[,com2] == com2, ]
      com2nodes <- c(rownames(see3))
      
      cg1 <- igraph::induced.subgraph(objgr, which(V(objgr)$name %in% com1nodes))
      cg2 <- igraph::induced.subgraph(objgr, which(V(objgr)$name %in% com2nodes))
      cegc1 <- igraph::degree(cg1, v = rownames(see)[i], mode = "all")
      cegc2 <- igraph::degree(cg2, v = rownames(see)[i], mode = "all")
      
      if(cegc1 > cegc2) {c1[i]<- com1}
      if(cegc1 < cegc2) {c1[i]<- com2}
      
      if (cegc1 == cegc2) {
        ogc1 <- igraph::degree(cg1, v = rownames(see)[i], mode = "out")
        ogc2 <- igraph::degree(cg2, v = rownames(see)[i], mode = "out")
        if(ogc1 > ogc2) {c1[i] <- com1}
        if(ogc1 < ogc2) {c1[i] <- com2}
        if(ogc1 == ogc2) {
          c1s <- gorder(cg1)
          c2s <- gorder(cg2)
          if (c1s > c2s) {c1[i] <- com1}
          if (c1s < c2s) {c1[i] <- com2}
          if (c1s == c2s) { c1[i] <- sample(par, 1)}
        }
        
      }
      
    }
    
    if(sum(!is.na(par)) > 2){ # in more than two communities - taking just first two
      par<-c(par[!is.na(par)])
      
      com1 <- par[1]
      com2 <- par[2]
      
      see2 <- see[see[,com1] == com1, ]
      com1nodes <- c(rownames(see2))
      see3 <- see[see[,com2] == com2, ]
      com2nodes <- c(rownames(see3))
      
      cg1 <- igraph::induced.subgraph(objgr, which(V(objgr)$name %in% com1nodes))
      cg2 <- igraph::induced.subgraph(objgr, which(V(objgr)$name %in% com2nodes))
      
      cegc1 <- igraph::degree(cg1, v = rownames(see)[i], mode = "all")
      cegc2 <- igraph::degree(cg2, v = rownames(see)[i], mode = "all")
      
      if(cegc1 > cegc2) {c1[i]<- com1}
      if(cegc1 < cegc2) {c1[i]<- com2}
      
      if (cegc1 == cegc2) {
        ogc1 <- igraph::degree(cg1, v = rownames(see)[i], mode = "out")
        ogc2 <- igraph::degree(cg2, v = rownames(see)[i], mode = "out")
        if(ogc1 > ogc2) {c1[i] <- com1}
        if(ogc1 < ogc2) {c1[i] <- com2}
        if(ogc1 == ogc2) {
          c1s <- gorder(cg1)
          c2s <- gorder(cg2)
          if (c1s > c2s) {c1[i] <- com1}
          if (c1s < c2s) {c1[i] <- com2}
          if (c1s == c2s) { c1[i] <- sample(par, 1)}
        }
        
      }
      
    }
    
  }
  
  #c1# this is community membership
  
  # which/how many students belong to two communities?
  da <- ifelse(rownames(see)  %in% cp$shared.nodes.labels, 1, 0) # 1 means it is a broker
  nbrok <- sum(da) # how many brokers
  # which/how many communities have members with shared membership?
  
  sharedcom <- list()
  for(i in 1:nrow(see)){
    par <- see[i,]
    if(sum(!is.na(par)) >1){
      par<-c(par[!is.na(par)])
      sharedcom[[i]] <- par
    }
    
  }
  
  sharedcom <- unique(unlist(sharedcom)) # communities that share
  nsharecom <- length(sharedcom)
  
  nodenames <- rownames(sw)
  EB1 <- as.data.frame(nodenames)
  EB1$membership <- c1
  
  ncom <- length(table(EB1$membership)) # incl.isolates as communities
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
  datcom$schid <- school # 1 school
  
  # saving the data
  va1<- vertex_attr(swg)
  va1df <- data.frame(va1)
  va1df$schid <- school # 1 school
  
  finalres <- list(basic = bas, data_com = datcom, dfres = va1df,
                   # additional results for CPM
                   memb.mat = see, vec.brokers = da,
                   NofBrokers = nbrok, sharing.com = sharedcom,
                   NofShareCom = nsharecom)
  
  return(finalres)
  
}

test1 <- net_pcm(palnet[[1]], school=1)
head(test1)
nrow(test1$memb.mat)
dim(test1$memb.mat)
eball2<-list()
for (i in 1:22){
  cda <- net_pcm(palnet[[i]], school =i) 
  alldata <- level1dat(resoff = cda, abb = "pcm")
  eball2[[length(eball2) + 1]] = alldata
}

findat <- rbindlist(eball2)
dim(findat)
colnames(findat)

write.table(findat, file = "PCM_res.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
write.xlsx(findat, "PCM_res.xlsx")

#-------------------------------------------------------------------------------
#              Saving additional PCM data in long format etc.
#-------------------------------------------------------------------------------

dtpcm<- list() # additional pcm info in long format for each school
macroinfo <- list() # macro info about PCM
for(i in 1:22){
  # loop for each net/school
  test1 <- net_pcm(palnet[[i]], school =i) 
  NofBrokers <- test1$NofBrokers # N
  NofShcom <- length(test1$sharing.com) # coms, but not ids
  test1$dfres$brokers <- test1$vec.brokers # 1 if it is a broker
  test1$dfres$com_share_2 <-ifelse(test1$dfres$E_B1 %in% test1$sharing.com, 1, 0 ) 
  # 1 if com shares members
  colnames(test1$dfres)
  test1$data_com
  # creating unique ID per com
  test1$data_com$uniqueID_2 <- test1$data_com$schid * 1000 + test1$data_com$com_name
  # map it
  comname <- test1$data_com$com_name
  comid <- test1$data_com$uniqueID_2
  test1$dfres$uniqueID_2 <- mapvalues(test1$dfres$E_B1,
                                      from = comname,
                                      to = comid)
  test1$dfres$E_B1 <- NULL
  
  dtpcm[[i]] <-test1$dfres
  
  macroinfo[[i]] <- c(i, NofBrokers, NofShcom)
}


findat2 <- rbindlist(dtpcm)
write.xlsx(findat2, "PCM_add_data_level1.xlsx")


