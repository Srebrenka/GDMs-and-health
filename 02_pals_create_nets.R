###################
#     02 PaLS     #
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

# create networks


#########################
#                       #
#    Load packages      #
#                       #
#########################

library(foreign)
library(dplyr)
library(statnet)
library(igraph)
library(intergraph)

#########################
#                       #
#     Main Script       #
#                       #
#########################
setwd("/Users/srebrenkaletina/Downloads")

#Load data
pals <- read.spss("PalsData.SAV", to.data.frame = TRUE)
colnames(pals) <- tolower(colnames(pals))

dim(pals)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#                 NETWORKS

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

palnet = list() 

schools = levels(factor(as.character(pals$schoolid))) %>%
  data.frame %>%
  setNames(c("school"))

for (i in schools$school){
  #Pull out individual schools#
  school = pals %>% 
    filter(schoolid == i)
  
  sx1 <- subset(school, select = c(palsqid, fr1id, fr2id,fr3id, fr4id, fr5id, fr6id))
  edge1 <- na_if(sx1, NA)
  edge1 <- as.data.frame(apply(edge1, 2, function(x) as.numeric(x)))
  edgeY <- data.frame()
  temp <- data.frame()
  for (k in 2:ncol(edge1)) {
    temp <- cbind(edge1[, 1], edge1[, k])
    edgeY <- rbind(edgeY, temp)
  }
  
  colnames(edgeY) <- c("respondent_id", "alter")
  edgecleanY <- edgeY[which(!is.na(edgeY$alter)),] 
  # remove duplicate rows - multiple edges
  edgecleanY <- edgecleanY[!duplicated(edgecleanY), ]
  # remove rows with equal value in the columns - self-loops
  edgecleanY <- edgecleanY[(!(edgecleanY$respondent_id == edgecleanY$alter)),]
  
  palnet[[i]] <-graph_from_data_frame(edgecleanY, directed = TRUE, vertices = NULL) 
  
}

plot(palnet[[1]], vertex.size = 3, vertex.label = "")

# checking for loops
for (i in 1:22){
  loop <- sum(which_loop(palnet[[i]]))
  print(loop)
}


##################################
# dealing with pupils being in more than one classrom

# are there same pupils in more than one net?
puplist=list()
for(i in 1:22){
  pupils <- list(V(palnet[[i]])$name)
  puplist[[i]] <- pupils
}

fl2 <- purrr::flatten(puplist)
fl2 <- purrr::flatten(fl2) 
length(fl2)# 3674
length(unique(fl2)) # 3649
3674 - 3649 # 25 appears at least twice

more <- duplicated(fl2)
sum(more ==T) # 25 are appearing in two networks

flc <- unique(fl2[duplicated(fl2)])
flc2<-unlist(flc) # the list

## checking

af <- flc2[2] # "34055"

for (i in 1:22){
  print(af %in% V(palnet[[i]])$name)} # appears in two schools

# in which schools a node appears
for (i in 1:length(flc2)) {
  nod<- flc[i]
  for (k in 1:22){
    a=k
    if (nod %in% V(palnet[[k]])$name) {
      print(c(nod, a))}
  } 
}

############ step by step *** 
############ delete the nodes that appear twice based on the first two letters

# "24142" in 8 and 9
V(palnet[[9]])$name
gorder(palnet[[9]])# 219
palnet[[9]] = igraph::delete.vertices(palnet[[9]], "24142")
# "34295" in 14 and 16
V(palnet[[16]])$name
gorder(palnet[[16]])# 278
palnet[[16]] = igraph::delete.vertices(palnet[[16]], "34295")
# "36029" in 16 and 17
V(palnet[[17]])$name
gorder(palnet[[17]])# 104
palnet[[17]] = igraph::delete.vertices(palnet[[17]], "36029" )
# "40040" in 18 and 19
V(palnet[[19]])$name
gorder(palnet[[19]])# 98
palnet[[19]] = igraph::delete.vertices(palnet[[19]], "40040"  )
# "40112" in 18 and 19
V(palnet[[19]])$name
palnet[[19]] = igraph::delete.vertices(palnet[[19]], "40112" )
# "40017" in 18 and 19
V(palnet[[19]])$name
palnet[[19]] = igraph::delete.vertices(palnet[[19]], "40017" )
# "40116" in 18 and 19
V(palnet[[19]])$name
palnet[[19]] = igraph::delete.vertices(palnet[[19]], "40116" )
# "40014" in 18 and 19
V(palnet[[19]])$name
palnet[[19]] = igraph::delete.vertices(palnet[[19]], "40014" )
# "40009" in 18 and 19
V(palnet[[19]])$name
palnet[[19]] = igraph::delete.vertices(palnet[[19]], "40009" )
# "40002" in 18 and 19
V(palnet[[19]])$name
palnet[[19]] = igraph::delete.vertices(palnet[[19]], "40002" )
# "40108" in 18 and 19
V(palnet[[19]])$name
palnet[[19]] = igraph::delete.vertices(palnet[[19]], "40108" )
# "40008" in 18 and 19
V(palnet[[19]])$name
palnet[[19]] = igraph::delete.vertices(palnet[[19]], "40008" )
# "40013" in 18 and 19
V(palnet[[19]])$name
palnet[[19]] = igraph::delete.vertices(palnet[[19]], "40013" )
# "40004" in 18 and 19
V(palnet[[19]])$name
palnet[[19]] = igraph::delete.vertices(palnet[[19]], "40004" )
# "40010" in 18 and 19
V(palnet[[19]])$name
palnet[[19]] = igraph::delete.vertices(palnet[[19]], "40010" )
# "42088" in 11 and 20
V(palnet[[11]])$name
gorder(palnet[[11]]) # 160
palnet[[11]] = igraph::delete.vertices(palnet[[11]], "42088"  )
# "42114" in 11 and 20
V(palnet[[11]])$name
palnet[[11]] = igraph::delete.vertices(palnet[[11]], "42114")
# "42041" in 11 and 20
V(palnet[[11]])$name
palnet[[11]] = igraph::delete.vertices(palnet[[11]], "42041")
# "43027" in 20 and 21
V(palnet[[20]])$name
gorder(palnet[[20]]) # 144
palnet[[20]] = igraph::delete.vertices(palnet[[20]], "43027" )
# "42064" in 20 and 22
V(palnet[[22]])$name
gorder(palnet[[22]]) # 136
palnet[[22]] = igraph::delete.vertices(palnet[[22]], "42064" )
V(palnet[[1]])$name
gorder(palnet[[1]]) #
palnet[[1]] = igraph::delete.vertices(palnet[[1]], "34055" )
V(palnet[[16]])$name
gorder(palnet[[16]]) #
palnet[[16]] = igraph::delete.vertices(palnet[[16]], "34196" )
V(palnet[[19]])$name
gorder(palnet[[19]]) #
palnet[[19]] = igraph::delete.vertices(palnet[[19]], "40012" )
V(palnet[[19]])$name
gorder(palnet[[19]]) #
palnet[[19]] = igraph::delete.vertices(palnet[[19]], "40088" )
V(palnet[[21]])$name
gorder(palnet[[21]]) #
palnet[[21]] = igraph::delete.vertices(palnet[[21]], "45105" )

#############################
#From igraph to network
Nets= list()
for (i in 1:22){
  Nets[[i]] <- asNetwork(palnet[[i]]) 
} 



