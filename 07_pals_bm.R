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

#  creating a function for blockmodeling and creating datasets
# executing these functions

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

#library(formattable)
#library(pracma)
#library(devtools)
library(openxlsx)
library(readxl)
library(blockmodeling)
library(sna)
library(plyr)



#########################
#                       #
#     Main Script       #
#                       #
#########################

# read in attribute data
rad1 <- read_excel("~/Downloads/attribute_data_pals2FIN.xlsx", col_names = T)

rad1 <- mutate_all(rad1, function(x) as.numeric(as.character(x)))
rad1$SEX <- as.factor(rad1$SEX)
rad1$PALSQID <- as.factor(rad1$PALSQID)
rad1$FAS  <- as.factor(rad1$FAS)

# creating a variable for matching
rad1$code <- rad1$PALSQID

pals$code <- pals$palsqid

pals2 <- merge(rad1,pals, by = "code")



set.seed(100111)

###################################################################

#####           FUNCTION BLOCK STRUCTURAL

net_block_s <- function(net, school, seedN, repet) {
  mat <- as.matrix(get.adjacency(net))
  set.seed(seedN)
  EB1 <-optRandomParC(M=mat, k=2, rep=repet, # increase it
                      approaches = "hom", 
                      homFun = "ss", # change to "ad"
                      blocks= "com",
                      mingr = 10)
  ncom <- 2
  
  V(net)$opt.blocks <- EB1$best$best1$clu
  
  bas <- c(school, ncom, gorder(net))
  
  comdat=list()
  for (k in unique(V(net)$opt.blocks)) {
    swg2 <- net
    ed1 <- gsize(swg2)
    swga<- igraph::delete.vertices(net, V(net)$opt.blocks != k)
    ed2 <- gsize(swga)
    # density of swga - community !!!!!!
    dencom <- round(edge_density(swga, loops = FALSE), 3)
    sizecom <- gorder(swga)
    swgb<- igraph::delete.vertices(net, V(net)$opt.blocks == k)
    ed3 <- gsize(swgb)
    
    edo <- ed1 - ed3 - ed2# careful, it can be zero
    if (edo==0) {
      rat <- 0
    } else {
      eda <- ed2 + edo
      rat <- round(edo/eda, 2)
    }
    
    eachcom <- c(sizecom, ed2, dencom, rat) 
    comdat[[length(comdat) + 1]] = eachcom
  }
  
  sq <- matrix(unlist(comdat), ncol = 4, byrow = TRUE)
  datcom <- as.data.frame(sq)
  colnames(datcom) = c("size", "edges","density","rat")
  datcom$cID <- 1:nrow(datcom)
  datcom$schid <- school
  
  # saving the data
  va1<- vertex_attr(net)
  va1df <- data.frame(va1)
  va1df$REGE.class <- NULL
  va1df$schid <- school
  
  finalres <- list(basic = bas, data_com = datcom, dfres = va1df)
  
  return(finalres)
}

zet = net_block_s(net = palnet[[22]], school = 22, seedN = 197845, repet = 10)

zet$basic
zet$data_com
zet$dfres

#####           FUNCTION BLOCK REGULAR

net_block_r <- function(net, school, seedN, repet) {
  mat <- as.matrix(get.adjacency(net))
  set.seed(seedN)
  REGE.cat<-REGE.ownm.for(M=mat)$E 
  set.seed(seedN)
  EB1 <-optRandomParC(M=REGE.cat, 
                      k=2,        # look for 2 classes
                      rep=repet,
                      approaches = "hom",
                      homFun = "ss",
                      blocks="com",
                      mingr = 10)
  ncom <- 2
  
  V(net)$REGE.class <- EB1$best$best1$clu
  
  bas <- c(school, ncom, gorder(net))
  
  comdat=list()
  for (k in unique(V(net)$REGE.class)) {
    swg2 <- net
    ed1 <- gsize(swg2)
    swga<- igraph::delete.vertices(net, V(net)$REGE.class != k)
    ed2 <- gsize(swga)
    # density of swga - community !!!!!!
    dencom <- round(edge_density(swga, loops = FALSE), 3)
    sizecom <- gorder(swga)
    swgb<- igraph::delete.vertices(net, V(net)$REGE.class == k)
    ed3 <- gsize(swgb)
    
    edo <- ed1 - ed3 - ed2# careful, it can be zero
    if (edo==0) {
      rat <- 0
    } else {
      eda <- ed2 + edo
      rat <- round(edo/eda, 2)
    }
    
    eachcom <- c(sizecom, ed2, dencom, rat) 
    comdat[[length(comdat) + 1]] = eachcom
  }
  
  sq <- matrix(unlist(comdat), ncol = 4, byrow = TRUE)
  datcom <- as.data.frame(sq)
  colnames(datcom) = c("size", "edges","density","rat")
  datcom$cID <- 1:nrow(datcom)
  datcom$schid <- school
  
  # saving the data
  va1<- vertex_attr(net)
  va1df <- data.frame(va1)
  va1df$schid <- school
  
  finalres <- list(basic = bas, data_com = datcom, dfres = va1df)
  
  return(finalres)
}

#####           FUNCTION BLOCK REGULAR bin

net_block_r_bin <- function(net, school, seedN, repet) {
  mat <- as.matrix(get.adjacency(net))
  set.seed(seedN)
  REGE.cat<-REGE.ownm.for(M=mat)$E 
  set.seed(seedN)
  EB1 <-optRandomParC(M=REGE.cat, 
                      k=2,        # look for 2 classes
                      rep=repet,
                      approaches = "bin",
                      
                      blocks="reg",
                      mingr = 10)
  ncom <- 2
  
  V(net)$REGE.class <- EB1$best$best1$clu
  
  bas <- c(school, ncom, gorder(net))
  
  comdat=list()
  for (k in unique(V(net)$REGE.class)) {
    swg2 <- net
    ed1 <- gsize(swg2)
    swga<- igraph::delete.vertices(net, V(net)$REGE.class != k)
    ed2 <- gsize(swga)
    # density of swga - community !!!!!!
    dencom <- round(edge_density(swga, loops = FALSE), 3)
    sizecom <- gorder(swga)
    swgb<- igraph::delete.vertices(net, V(net)$REGE.class == k)
    ed3 <- gsize(swgb)
    
    edo <- ed1 - ed3 - ed2# careful, it can be zero
    if (edo==0) {
      rat <- 0
    } else {
      eda <- ed2 + edo
      rat <- round(edo/eda, 2)
    }
    
    eachcom <- c(sizecom, ed2, dencom, rat) 
    comdat[[length(comdat) + 1]] = eachcom
  }
  
  sq <- matrix(unlist(comdat), ncol = 4, byrow = TRUE)
  datcom <- as.data.frame(sq)
  colnames(datcom) = c("size", "edges","density","rat")
  datcom$cID <- 1:nrow(datcom)
  datcom$schid <- school
  
  # saving the data
  va1<- vertex_attr(net)
  va1df <- data.frame(va1)
  va1df$schid <- school
  
  finalres <- list(basic = bas, data_com = datcom, dfres = va1df)
  
  return(finalres)
}

zet = net_block_r(net = palnet[[22]], school = 22, seedN = 197845, repet = 10)

zet$basic
zet$data_com
zet$dfres

sball<-list()
for (i in 1:22){
  cmeb <- net_block_s(palnet[[i]], school =i, seedN = 197845, repet = 100) 
  sball[[length(sball) + 1]] = cmeb
  print(i)
}
sball[[1]]$dfres

#################################################
# var = opt.blocks in s, REGE.class in r
creat_dat <- function(the_list){
  
  eb.basic <- list()
  for (i in 1:22) {
    ebb <- the_list[[i]]$basic
    eb.basic[[length(eb.basic) +1]] = ebb
  }
  sq1 <- matrix(unlist(eb.basic), ncol = 3, byrow = TRUE)
  # level 3 - info per networks
  basic1 <- as.data.frame(sq1)  
  colnames(basic1) = c("school_3","Ncom_3", "netsize_3") # netsize 
  # level 2 - per communities
  eb.data_com <- list()
  for (i in 1:22) {
    ebb <- the_list[[i]]$data_com
    eb.data_com[[length(eb.data_com) +1]] = ebb
  }
  datacom1 <- ldply(eb.data_com, data.frame)
  
  colnames(datacom1) = c("com_size_2","edges_2", "com_density_2", 
                         "com_rat_2", "cID_2", "schid_2")
  datacom1$new_schid_2 <- datacom1$schid_2 * 1000 + datacom1$cID_2
  
  # level1
  eb.dfres <- list()
  for (i in 1:22) {
    ebb <- the_list[[i]]$dfres
    ebb$newEB <- i*1000 + ebb$REGE.class # to have unique ids in whole dataset
    eb.dfres[[length(eb.dfres) +1]] = ebb
  }
  edfs <- ldply(eb.dfres, data.frame)
  colnames(edfs) = c("name_1","COMid_1", "schid_1", "newCOMid_1")
  
  # merge with attribute data
  
  pals3 <- subset(pals2, select = c(code,
                                    SCHOOLID,
                                    PALSQID,
                                    SEX,
                                    DOBYR, 
                                    FAS,
                                    ethnicity,
                                    PBICONT,
                                    PBICARE,
                                    PC1,
                                    PC2,
                                    com1,
                                    com2))
  
  edfs$code <- edfs$name
  merdat <- merge(edfs, pals3, by="code", all.x = TRUE)
  # join all levels
  # from level 2 to level 1
  for (i in colnames(datacom1)){
    merdat[,i] <- mapvalues(merdat$newCOMid_1,
                            from = c(datacom1$new_schid_2),
                            to = c(datacom1[,i]))
  }
  # from level 3 to level 1
  
  for (i in colnames(basic1)){
    merdat[,i] <- mapvalues(merdat$schid_1,
                            from = c(basic1$school_3),
                            to = c(basic1[,i]))
  }
  merdat$SCHOOLID <- NULL # different school id
  merdat$schid_1 <- NULL # same as school_3
  merdat$schid_2 <- NULL # same as school_3
  merdat$un_com_id_2 <- merdat$new_schid_2 # renaming
  merdat$new_schid_2 <- NULL
  merdat$code <- NULL # same as name1, while PALSQID has missing 
  merdat$cID_2 <- NULL # same as COMid_1 - redundant, but leaving it in
  return(merdat)
}

# structural equivalence

sball<-list()
for (i in 1:22){
  cmeb <- net_block_s(palnet[[i]], school =i, seedN = 197845, repet = 100) 
  sball[[length(sball) + 1]] = cmeb
  
}
length(sball)
s_res <- creat_dat(sball)  

dim(s_res)

write.table(s_res, file = "blockS_res.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
write.xlsx(s_res, "blockS_res.xlsx") # *** do some checking

################

# regular equivalence 

rball<-list()
for (i in 1:22){
  cmeb <- net_block_r(palnet[[i]], school =i, seedN = 197845, repet = 100) 
  rball[[length(rball) + 1]] = cmeb
  print(i) 
  
}

r_res <- creat_dat(rball)  

dim(r_res)

write.table(r_res, file = "blockR_res.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
write.xlsx(r_res, "blockR_res.xlsx")

# regular equivalence binary

rball<-list()
for (i in 1:22){
  cmeb <- net_block_r_bin(palnet[[i]], school =i, seedN = 197845, repet = 100) 
  rball[[length(rball) + 1]] = cmeb
  print(i) 
  
}

r_res <- creat_dat(rball)  

dim(r_res)

write.table(r_res, file = "blockRbin_res.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
write.xlsx(r_res, "blockRbin_res.xlsx")
