###################
#     03 PaLS     #
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

# network properties and plots
# adding attribute gender


#########################
#                       #
#     Main Script       #
#                       #
#########################
setwd("/Users/srebrenkaletina/Downloads")

source("/Users/srebrenkaletina/Downloads/02_pals_create_nets.R") 

plot
#########################
#                       #
#    Loaded packages      #
#                       #
#########################
#library(foreign)
#library(dplyr)
#library(statnet)
#library(igraph)
#library(intergraph)

# load in this script:
library(readxl)
library(car)
# ?:
library(Hmisc)
library(furniture)
library(magrittr)
library(lubridate)
#######################################
#########  Basic network descriptives


for (i in 1:22) {
  print(igraph::is.connected(palnet[[i]]))
}
# 6 out of 22 are connected

# searching for isolates
for (i in 1:22) {
  iso_nodes = length(isolates(Nets[i][[1]]))
  print(iso_nodes)
} 
# No isolates - bc it is constructed from edge-list!!!!!!!

# calculating network properties

ListNP = list()

for (i in 1:22) {
  
  N = network.size(Nets[i][[1]])
  edges = igraph::ecount(palnet[i][[1]])
  iso_nodes = length(isolates(Nets[i][[1]]))
  density = round(gden(Nets[i][[1]])[1][[1]], 3)
  reciprocity = round(igraph :: reciprocity(palnet[i][[1]])[1][[1]], 3)
  transitivity = round(igraph:: transitivity(palnet[i][[1]])[1][[1]], 3)
  avg_out_deg = round(mean(igraph::degree(palnet[i][[1]], mode = 'out', loops = F, normalized = F)), 2)
  avg_in_deg = round(mean(igraph::degree(palnet[i][[1]], mode = 'in', loops = F, normalized = F)), 2)
  
  tranL <- mean(transitivity(palnet[[i]], type = "local"), na.rm = T)
  size<-gorder(palnet[[i]])
  degreeM <-mean(igraph::degree(palnet[[i]]), na.rm = T)
  central <-centr_degree(palnet[[i]])$centralization
  C <- igraph::is.connected(palnet[[i]])
  
  if (C == FALSE) {
    asp<-mean_distance(palnet[[i]], directed = F, unconnected = TRUE)
  } else {
    asp<-mean_distance(palnet[[i]], directed = F, unconnected = FALSE)
  }
  
  
  resL = list(N, size, edges, iso_nodes, density, reciprocity, transitivity, 
              tranL, avg_out_deg, avg_in_deg, degreeM, asp, central, C)
  
  ListNP[[length(ListNP) + 1]] = resL
}


s <- matrix(unlist(ListNP), ncol = 14, byrow = TRUE)
netdat <- as.data.frame(s)

colnames(netdat) = c("N", "size","edges", "isolates", "density", "reciprocity", 
                     "transitivityG","transitivityL", "avg_out_degree", 
                     "avg_in_degree", "all_degree", "ASP", "centralization" ,"connected")

write.table(netdat, file = "Network_properties_PaLSdec.csv", dec = ',' ,
            sep = ";",row.names = TRUE)
write.xlsx(netdat, "Network_properties_PaLSdec.xlsx")

# more nodes than data rows in pals.sav - bc students that 
# did not participate in the study have been nominated as friends

# NOTE:
# the true number of schoolmates is not know from the data
# some may be not in the attribute data and neither in network

# sanity check:
# check N of missing all friends (2 potential different reasons:
# not having friends or skipping that part (sna part) of the survey 


nams1=list()
nams2=list()
for (i in schools$school){
  school = pals %>% 
    filter(schoolid == i)
  net <- school %>% 
    dplyr::select(palsqid, fr1id, fr2id,fr3id, fr4id, fr5id, fr6id) 
  
  nam1 <- sum(net$fr1id == "missing - 20 missed all friends items")
  nams1[[i]] <- nam1
  nam2 <- sum(net$fr2id == "missing - 20 missed all friends items")
  nams2[[i]] <- nam2
}

nams1
nams2

# nams1 is the same as nams2, ok

# attribute data: gender
sum(is.na(pals$sex)) # 0
table(pals$sex)

# gender per school table
gen=list()
for (i in schools$school){
  #Pull out individual schools#
  school = pals %>% 
    filter(schoolid == i)
  females <- sum(school$sex=="female")
  males <- sum(school$sex=="male")
  schoolsize <- nrow(school)
  
  finres <- list(females, males, schoolsize)
  
  gen[[length(gen) + 1]] = finres
}

sq <- matrix(unlist(gen), ncol = 3, byrow = TRUE)
gendat <- as.data.frame(sq)

colnames(gendat) = c("females", "males","N")

head(gendat)
sum(gendat$females)

write.table(gendat, file = "genderSchools_PaLSdec.csv", dec = ',' ,
            sep = ";",row.names = TRUE)



# reading in attribute data

rad1 <- read_excel("~/Downloads/attribute_data_pals2_FIN.xlsx", col_names = T)

rad1 <- mutate_all(rad1, function(x) as.numeric(as.character(x)))
colnames(rad1)

# convert to factors
rad1$gender <- as.factor(rad1$gender)
rad1$PALSQID <- as.factor(rad1$PALSQID)
rad1$FAS  <- as.factor(rad1$FAS)

# creating a variable for matching
rad1$code <- rad1$PALSQID

pals$code <- pals$palsqid

pals2 <- merge(rad1,pals, by = "code")

table(pals2$gender)
# 0 is male


# add gender attribute
for (i in schools$school){
  #Pull out individual schools#
  school = pals2 %>% 
    filter(schoolid == i)
  school$name <- as.character(school$palsqid)
  tomatch <- school[school$name %in% V(palnet[[i]])$name,]
  tomatch <- tomatch[match(V(palnet[[i]])$name, tomatch$name), ]
  
  V(palnet[[i]])$gender <- tomatch$gender
}

V(palnet[[1]])$gender # added ok 

#>>>>>>>>>>>>>>>>>>>>>>>>  FIGURES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ploting one network (school 13)
V(palnet[[13]])$size <- 5
V(palnet[[13]])$color <- "orange"

plot(palnet[[13]], vertex.label='',
     edge.arrow.size = 0.2,
     arrow.mode = "-", 
     edge.color = "gray" )

# with gender info

V(palnet[[13]])$color <- ifelse(is.na(V(palnet[[13]])$gender)== TRUE, 'grey',
                                ifelse(V(palnet[[13]])$gender == 1, 'orange',
                                       ifelse(V(palnet[[13]])$gender == 0, 
                                              "skyblue",'black')))

# choosing layout
#lay <- layout.kamada.kawai(palnet[[13]]) # nice
#lay <- layout.fruchterman.reingold(palnet[[13]]) # nice
lay <- layout_nicely(palnet[[13]]) # I choose this one
plot(palnet[[13]], 
     vertex.size=3, 
     vertex.label="", 
     edge.arrow.size = 0.2,
     #     edge.arrow.width = 0.1,
     arrow.mode = "-",
     edge.color = "gray",
     vertex.color = V(palnet[[13]])$color,
     layout = lay, main = "Friendship network 13") 

# creating figures for all schools - with gender attribute
### pdf
pdf(file="all22decX.pdf",
    width = 10, height = 14)
par(mfrow=c(6,4), mar=c(0.9,0.9,0.9,0.9))
for (i in 1:length(palnet)) { 
  lay <- layout_nicely(palnet[[i]])
  V(palnet[[i]])$color <- ifelse(is.na(V(palnet[[i]])$gender)== TRUE, 'grey',
                                 ifelse(V(palnet[[i]])$gender == 1, 'orange',
                                        ifelse(V(palnet[[i]])$gender == 0, 
                                               "skyblue",'black')))
  plot(palnet[[i]], vertex.label='',
       vertex.size=6.5,
       vertex.label="", 
       edge.arrow.size = 0.2,
       arrow.mode = "-", 
       edge.color = "gray",
       vertex.color = V(palnet[[i]])$color,
       layout = lay,
       main = paste("School", i)
  ) }

dev.off()

##### png

png(filename="all22decX.png",
    width = 1220, height = 1220, units = "px")
par(mfrow=c(5,5))
for (i in 1:length(palnet)) { 
  lay <- layout_nicely(palnet[[i]])
  V(palnet[[i]])$color <- ifelse(is.na(V(palnet[[i]])$gender)== TRUE, 'grey',
                                 ifelse(V(palnet[[i]])$gender == 2, 'orange',
                                        ifelse(V(palnet[[i]])$gender == 1, 
                                               "skyblue",'black')))
  plot(palnet[[i]], vertex.label='',
       vertex.size=6,
       vertex.label="", 
       edge.arrow.size = 0.2,
       arrow.mode = "-", 
       edge.color = "gray",
       vertex.color = V(palnet[[i]])$color,
       layout = lay,
       main = paste("School", i)
  ) }
dev.off()

# png: all nets one by one 
png(filename="pals22nets%d.png")
for (i in 1:length(palnet)) { 
  lay <- layout_nicely(palnet[[i]])
  V(palnet[[i]])$color <- ifelse(is.na(V(palnet[[i]])$gender)== TRUE, 'grey',
                                 ifelse(V(palnet[[i]])$gender == 2, 'orange',
                                        ifelse(V(palnet[[i]])$gender == 1, 
                                               "skyblue",'black')))
  plot(palnet[[i]], vertex.label='',
       vertex.size=5,
       vertex.label="", 
       edge.arrow.size = 0.2,
       arrow.mode = "-", 
       edge.color = "gray",
       vertex.color = V(palnet[[i]])$color,
       layout = lay,
       main = paste("School", i)
  ) }
dev.off()






