###################
#     PaLS viz     #
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

# figure of one school with different GDMs


#########################
#                       #
#    Load packages      #
#                       #
#########################
library(readxl)
library(ggplot2)
# read in from another script

#########################
#                       #
#     Main Script       #
#                       #
#########################

source("/02_pals_create_nets.R")
# read in all GDMs partitions
mypath= "/PalsComSets"
multimerge = function(mypath){
  filenames=list.files(mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x) {read_excel(path = x, col_names = T)})
  Reduce(function(x,y) {merge(x,y, by="id_1")}, datalist)}# name is later id_1

tbl <- multimerge(mypath) 
dim(tbl)
head(tbl)
# Create CDA specific dataset
cdadat <- subset(tbl, select = c(id_1,  
                                 wd_school_id_3,
                                 bia_uniqueID_2,
                                 pcm_uniqueID_2,
                                 ed_uniqueID_2,
                                 cfg_uniqueID_2,
                                 imd_uniqueID_2,
                                 le1_uniqueID_2,
                                 lo_uniqueID_2,
                                 lpu_uniqueID_2,
                                 sbm_uniqueID_2,
                                 wd_uniqueID_2))

dim(cdadat)
WorkNames <- c("BIA", "PCM",  "ED", "CFG","IMD", "LE", "LO", "LPU", "SBM", "WD")
corNames <- c("BIA", "CP", "EB", "FG", "IM", "LE", "LO", "LP", "SBM", "WT")

# my way:
#choose school 3
cdadat3 <- cdadat[cdadat$wd_school_id_3 ==3, ]

nc <-c() # list of N of communities
allmemb_list <-list() # list of communities for each GDA
for(i in 1:10){
  k <- i+2
  cddat <- c(cdadat3[,k])
  ids <- c(cdadat3$id_1)
  ndf <-cbind(ids, cddat)
  ndf <- as.data.frame(ndf)
  colnames(ndf) <- c("id", "gda")
  #ndf$id <- as.numeric(ndf$id)
  #ndf$gda <- as.numeric(ndf$gda)
  nc[i] <-length(unique(ndf$gda)) 
  gr <- unique(ndf$gda)
  list_mem<- list()
  for(g in 1:length(gr)){
    coms2<- ndf[ndf$gda == gr[g], ]
    memb <- coms2$id
    list_mem[[g]] <- memb
  }
  allmemb_list[[i]]<- list_mem
}
nc
set.seed(100111)
plot(palnet[[3]], mark.groups=list_mem,
     vertex.size=5, 
     vertex.label="")

# make figure
pdf(file="school3_GDAs.pdf",
    width = 10, height = 5)
par(mfrow=c(2,5), mar=c(0.9,0.9,0.9,0.9))
for (i in 1:10) { 
  lay <- layout_nicely(palnet[[3]])
  set.seed(100111)
  plot(palnet[[3]], mark.groups=allmemb_list[[i]],
       vertex.label='',
       vertex.size=6.5,
       vertex.label="", 
       edge.arrow.size = 0.2,
       arrow.mode = "-", 
       #edge.color = "gray",
       layout = lay,
       main = paste(corNames[i], nc[i], "communities")
  ) }
dev.off()

#choose school 15 
cdadat15 <- cdadat[cdadat$wd_school_id_3 ==15, ]

nc <-c() # list of N of communities
allmemb_list <-list() # list of communities for each GDA
for(i in 1:10){
  k <- i+2
  cddat <- c(cdadat15[,k])
  ids <- c(cdadat15$id_1)
  ndf <-cbind(ids, cddat)
  ndf <- as.data.frame(ndf)
  colnames(ndf) <- c("id", "gda")
  #ndf$id <- as.numeric(ndf$id)
  #ndf$gda <- as.numeric(ndf$gda)
  nc[i] <-length(unique(ndf$gda)) 
  gr <- unique(ndf$gda)
  list_mem<- list()
  for(g in 1:length(gr)){
    coms2<- ndf[ndf$gda == gr[g], ]
    memb <- coms2$id
    list_mem[[g]] <- memb
  }
  allmemb_list[[i]]<- list_mem
}
nc


# make figure
pdf(file="school15_GDAs.pdf",
    width = 10, height = 5)
par(mfrow=c(2,5), mar=c(0.9,0.9,0.9,0.9))
for (i in 1:10) { 
  lay <- layout_nicely(palnet[[15]])
  set.seed(100111)
  plot(palnet[[15]], mark.groups=allmemb_list[[i]],
       vertex.label='',
       vertex.size=6.5,
       vertex.label="", 
       edge.arrow.size = 0.2,
       arrow.mode = "-", 
       #edge.color = "gray",
       layout = lay,
       main = paste(corNames[i], nc[i], "communities")
  ) }
dev.off()

#choose school 19
cdadat19 <- cdadat[cdadat$wd_school_id_3 ==19, ]

nc <-c() # list of N of communities
allmemb_list <-list() # list of communities for each GDA
for(i in 1:10){
  k <- i+2
  cddat <- c(cdadat19[,k])
  ids <- c(cdadat19$id_1)
  ndf <-cbind(ids, cddat)
  ndf <- as.data.frame(ndf)
  colnames(ndf) <- c("id", "gda")
  #ndf$id <- as.numeric(ndf$id)
  #ndf$gda <- as.numeric(ndf$gda)
  nc[i] <-length(unique(ndf$gda)) 
  gr <- unique(ndf$gda)
  list_mem<- list()
  for(g in 1:length(gr)){
    coms2<- ndf[ndf$gda == gr[g], ]
    memb <- coms2$id
    list_mem[[g]] <- memb
  }
  allmemb_list[[i]]<- list_mem
}
nc


# make figure
pdf(file="school19_GDAs.pdf",
    width = 10, height = 5)
par(mfrow=c(2,5), mar=c(0.9,0.9,0.9,0.9))
for (i in 1:10) { 
  lay <- layout_nicely(palnet[[19]])
  set.seed(100111)
  plot(palnet[[19]], mark.groups=allmemb_list[[i]],
       vertex.label='',
       vertex.size=6.5,
       vertex.label="", 
       edge.arrow.size = 0.2,
       arrow.mode = "-", 
       #edge.color = "gray",
       layout = lay,
       main = paste(corNames[i], nc[i], "communities")
  ) }
dev.off()

