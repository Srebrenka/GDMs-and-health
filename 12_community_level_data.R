###################
#     PaLS 12     #
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

# adding community level data


#########################
#                       #
#    Load packages      #
#                       #
#########################
library(plyr)
library(openxlsx)
library(readxl)

source("02_pals_create_nets.R") 


mypath= "/PalsComSets"
multimerge = function(mypath){
  filenames=list.files(mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x) {read_excel(path = x, col_names = T)})
  Reduce(function(x,y) {merge(x,y, by="id_1")}, datalist)}# name is later id_1

tbl <- multimerge(mypath) 
dim(tbl)


allcdas <- list(tbl$bia_uniqueID_2,
                tbl$cfg_uniqueID_2,
                tbl$ed_uniqueID_2,
                tbl$imd_uniqueID_2,
                tbl$le1_uniqueID_2,
                tbl$lo_uniqueID_2,
                tbl$lpu_uniqueID_2,
                tbl$pcm_uniqueID_2,
                tbl$sbm_uniqueID_2,
                tbl$wd_uniqueID_2) #,

namecdas <- list("bia","cfg","ed","imd","le1",
                 "lo","lpu","pcm","sbm","wd")

length(allcdas) 

start_time <- Sys.time()
datasets <- list()
cdainf <- list()
for(i in 1:length(allcdas)){
  gei <- get_ei(allcdas[[i]], var = "gender")
  gei$school_id <- NULL
  gei$cda <- namecdas[i]
  datasets[[i]] <-gei
  dc <- div_com(gei)
  cdainf[[i]] <- dc
}
end_time <- Sys.time()
end_time - start_time # Time difference of 1.027016 mins

length(datasets) # this is not a long format - all CDAs

findat2 <- rbindlist(datasets)
dim(findat2) # 6921   11
write.xlsx(findat2, "Com_level_data_allCDA_1.xlsx")
findat3 <- as.data.frame(do.call(base::rbind, cdainf))
dim(findat3) # 17 11
colnames(findat3) <- c("NCom", "Avg.com.size", "Max.com.size", "Min.com.size",
                       "NCsize_1", "NCsize_2", "NCsize_3",
                       "NCsize4_12","NCsize13_30", "NC31plus", "NcomWnon_res",
                       "Nofdiverse",  "perOfDivCom","CasesInDivCom", 
                       "perOfcasesInDiCom")
findat3$cda <- namecdas

findat3 <- findat3[, c(16, 1:15)]
write.xlsx(findat3, "CDA_info1.xlsx")

c1 = tbl[variabletbl == uc[i], ] # argument 1
nodeinc <-c1$id_1

sid <- c1$wd_school_id_3[4]
exg2 <- igraph::induced.subgraph(palnet[[sid]], which(V(palnet[[sid]])$name %in% nodeinc))
nnet <- gorder(exg2)

get_tau <- function(variabletbl){
  # variabletbl is column name for cda and its unique id, e.g., tbl$wd_uniqueID_2
  # var is variable on basis of which EI is calculated - needs to be added as attribute to palnets
  # tbl has to be open
  uc = unique(variabletbl)
  finres <- list()
  for (i in 1:length(uc)){
    c1 = tbl[variabletbl == uc[i], ] # argument 1
    nodeinc <-c1$id_1
    sid <- c1$wd_school_id_3[1]
    exg2 <- igraph::induced.subgraph(palnet[[sid]], which(V(palnet[[sid]])$name %in% nodeinc))
    nnet <- gorder(exg2)
    if((gsize(exg2) >2)) {
      
      # transitivity
      den2 = igraph::edge_density(exg2)
      netran <- igraph::transitivity(exg2, type = "global")
      trans_r <- c()
      
      for (l in 1:250) {
        set.seed(74590 + l)
        erg <- erdos.renyi.game(nnet, den2, type = "gnp", directed = TRUE)
        ergt <- igraph::transitivity(erg, type = "global")
        trans_r[l] <- ergt
        
      }
      
      trans_er <- mean(trans_r, na.rm = T)
      sd_trans_er <- round(sd(trans_r, na.rm = T), 3) 
      zTrans <-round((netran - trans_er)/sd_trans_er, 2)
      
      #  Producing a tau statistic
      tc_t <- triad_census(exg2)
      
      weighting_scheme <- c(0,0,0,1,1,-1,0,0,1,0,0,1,1,0,0,0) # for tau
      
      #  Producing a tau statistic
      tau_r <- c()
      outd <- igraph::degree(exg2,mode = "out")
      ind <- igraph::degree(exg2,mode = "in")
      for (k in 1:250) {
        # Generate random graphs with a given degree sequence
        confg <- sample_degseq(out.deg = outd, in.deg = ind, 
                               method = "simple.no.multiple" )
        tau_c <- sum(triad_census(confg) * weighting_scheme)
        tau_r[k] <- tau_c
        
      }
      tau_t <- sum(tc_t * weighting_scheme)
      # tau score (normalized with configuration models)
      tau_er <- mean(tau_r)
      sd_tau_er <- round(sd(tau_r), 3) 
      zTau <-round((tau_t - tau_er)/sd_tau_er, 2)
      
      # centralization (normalized?) 
      cdn <- centr_degree(
        graph = exg2,
        mode = "all",
        loops = FALSE,
        normalized = T)
      cnt <- cdn$centralization
      
    }
    else{
      tau_t <- NA
      zTrans <- NA
      zTau <- NA
      cnt <- NA
    }
    
    finres[[i]] <- c(sid, uc[i], round(netran, 3),round(zTrans, 3), round(tau_t, 3), round(zTau, 3), round(cnt, 3))
  }
  
  frdat <- as.data.frame(do.call(rbind, finres))
  colnames(frdat) <- c("school","uniqueID", "Transit", "zTrans", "Tau", "zTau", "CntN")
  return(frdat)
}

start_time <- Sys.time()
ts1 <- get_tau(tbl$wd_uniqueID_2)
end_time <- Sys.time()
end_time - start_time # Time difference of 1.55425 mins

head(ts1)
cor.test(ts1$Transit, ts1$zTrans)

sum(is.na(ts1$Transit))


start_time <- Sys.time()
datasetau <- list()
for(i in 1:length(allcdas)){
  gei <- get_tau(allcdas[[i]])
  gei$cda <- namecdas[i]
  datasetau[[i]] <-gei
}
end_time <- Sys.time()
end_time - start_time # Time difference of 30.43271 mins

findat4 <- rbindlist(datasetau)
dim(findat4) #  6921    7
head(findat4)
write.xlsx(findat4, "Com_level_data_allCDA_2_tauetc.xlsx")

histogram(ts1$zTrans)

histogram(gei1$ei)
histogram(gei1$orwg)
histogram(gei1$assort)
sum(is.na(gei1$orwg)) # 308
sum(is.na(gei1$ei)) # 31
sum(is.na(gei1$assort)) # 306!!!
max(gei1$ei, na.rm = T)
sum(gei1$net_prob) # 305

ge1 <- gei1[gei1$net_prob == 0,]




# get mixing parameter

get_mix_p <- function(variabletbl){
  uc <- unique(variabletbl)
  finres <- list()
  for(i in 1:length(uc)){
    c1 = tbl[variabletbl== uc[i], ]
    nodeinc <-c1$id_1
    sid <- c1$wd_school_id_3[1]
    ex <- palnet[[sid]]
    exg2 <- igraph::induced.subgraph(palnet[[sid]], which(V(palnet[[sid]])$name %in% nodeinc))
    nnet <- gorder(exg2)
    # lists:
    extdX <-list()
    totX <- list()
    mpnX <- list()
    for(k in 1:nnet){
      dt <- igraph::degree(ex, nodeinc[k], mode = "all")
      dt <- as.numeric(dt)
      dc <- igraph::degree(exg2, nodeinc[k], mode = "all")
      dc <- as.numeric(dc)
      extd <- dt- dc
      tot <- dt
      # node level mixing parameter (added)
      mpn <- round(extd/tot, 2)
      extdX[[k]] <- extd
      totX[[k]] <- tot
      mpnX[[k]] <- mpn
    }
    
    out <- as.data.frame(do.call(cbind, list(nodeinc, extdX, totX, mpnX)))
    colnames(out) <- c("id_1", "extd", "tot", "mpn")
    out$uniqueIDcom <- uc[i]
    out$school <- sid
    
    finres[[i]] <- out
  }
  finres <- as.data.frame(rbindlist(finres))
  return(finres)
}

ah <- get_mix_p(allcdas[[1]])
class(ah)
head(ah)


datasetmp <- list()
for(i in 1:length(allcdas)){
  gei <- get_mix_p(allcdas[[i]])
  gei$cda <- namecdas[i]
  datasetmp[[i]] <-gei
}

length(datasetmp) # 17 cdas
datasetmp2 <- rbindlist(datasetmp)
dim(datasetmp2)

write.xlsx(datasetmp2, "mix_par_all_cda.xlsx")

# get some figures
uc[55] # 5019
c1 = tbl[tbl$wd_uniqueID_2 == 5019, ]
nodeinc <-c1$id_1


uc <- unique(tbl$wd_uniqueID_2)
for(i in 300:330){
  c1 = tbl[tbl$wd_uniqueID_2 == uc[i], ]
  nodeinc <-c1$id_1
  print(i)
  print(length(nodeinc))
}
c1 = tbl[tbl$wd_uniqueID_2 == uc[320], ]
nodeinc <-c1$id_1
c1$wu_school_id_3[1]

#### EXAMPLE
bbb<- palnet[[19]]
exg2 <- igraph::induced.subgraph(bbb, which(V(bbb)$name %in% nodeinc))
plot(exg2, vertex.label="", vertex.color = "grey")
transitivity(exg2, type = "global")
uc[320]
txxx <- ts1[ts1$uniqueID == "19006",]
txxx$Transit

#Calculate transitivity again for WT
uc <- unique(tbl$wd_uniqueID_2)
com_ids<- c()
true.t<- c()
old.t<- c()
for(i in 1:length(uc)){
  c1 = tbl[tbl$wd_uniqueID_2 == uc[i], ]
  nodeinc <-c1$id_1
  sc <- c1$wu_school_id_3[1]
  bbb<- palnet[[sc]]
  exg2 <- igraph::induced.subgraph(bbb, which(V(bbb)$name %in% nodeinc))
  Transitiv <-transitivity(exg2, type = "global") 
  true.t[i] <- Transitiv
  id <-uc[i]
  txxx <- ts1[ts1$uniqueID == id,]
  txxx$Transit
  com_ids[i] <- id
  old.t[i] <- txxx$Transit
}

check <- as.data.frame(cbind(com_ids, true.t, old.t))
cor.test(check$true.t, check$old.t) # 0.90
head(check)
write.xlsx(check, "TrueTransitivityWT.xlsx")

nrow(check)


# calculate true transitivity for all 10 GDM
uc <- unique(tbl$wd_uniqueID_2)
com_ids<- c()
true.t<- c()
for(i in 1:length(uc)){
  c1 = tbl[tbl$wd_uniqueID_2 == uc[i], ]
  nodeinc <-c1$id_1
  sc <- c1$wu_school_id_3[1]
  bbb<- palnet[[sc]]
  exg2 <- igraph::induced.subgraph(bbb, which(V(bbb)$name %in% nodeinc))
  Transitiv <-transitivity(exg2, type = "global") 
  true.t[i] <- Transitiv
  id <-uc[i]
  com_ids[i] <- id
}
frdat <- cbind(com_ids, true.t)
frdat <- as.data.frame(frdat)
frdat$true.t <- round(frdat$true.t, 3)
head(frdat)

# the function
get_trans <- function(variabletbl){
  # variabletbl is column name for cda and its unique id, e.g., tbl$wd_uniqueID_2
  # var is variable on basis of which EI is calculated - needs to be added as attribute to palnets
  # tbl has to be open
  uc <- unique(variabletbl) # tbl$wd_uniqueID_2
  com_ids<- c()
  true.t<- c()
  for(i in 1:length(uc)){
    c1 = tbl[variabletbl == uc[i], ]# tbl$wd_uniqueID_2
    nodeinc <-c1$id_1
    sc <- c1$wu_school_id_3[1]
    bbb<- palnet[[sc]]
    exg2 <- igraph::induced.subgraph(bbb, which(V(bbb)$name %in% nodeinc))
    Transitiv <-transitivity(exg2, type = "global") 
    true.t[i] <- Transitiv
    id <-uc[i]
    com_ids[i] <- id
  }
  
  
  frdat <- cbind(com_ids, true.t)
  frdat <- as.data.frame(frdat)
  frdat$true.t <- round(frdat$true.t, 3)
  colnames(frdat) <- c("com_id", "Transit")
  return(frdat)
}

comdat2 <- as.data.frame(read_excel("Com_level_data_allCDA_2_tauetc.xlsx"))
comdat2 = comdat2[comdat2$cda == smallN[10], ]
head(comdat2)
nrow(comdat2) # 387
sum(is.na(comdat2$Tau)) # 48
sum(is.na(comdat2$zTau)) # 151

cdn <- comdat2[!is.na(comdat2$Tau) & is.na(comdat2$zTau), ]
nrow(cdn) # 103
head(cdn)

comstolook <- cdn$uniqueID

uc = unique(allcdas[[10]]) # wt
ch1 <- get_trans(allcdas[[2]])
naTra <- c()
for(i in 1:nrow(ch1)){
  if(is.na(ch1[i,2])){
    naTra[i] <- i} else{naTra[i]<- NA}
}

naTra <- naTra[!is.na(naTra)]
naTra


c1 = tbl[allcdas[[10]] == uc[1], ] #comstolook[2]
nodeinc <-c1$id_1
length(nodeinc)
sc <- c1$wu_school_id_3[1]
bbb<- palnet[[sc]]
exg2 <- igraph::induced.subgraph(bbb, which(V(bbb)$name %in% nodeinc))
Transitiv <-igraph::transitivity(exg2, type = "global") 
Transitiv
plot(exg2, vertex.size = 3, vertex.name = "")
naTra

cor.test(cdadat$Tau, cdadat$zTau_old)

sum(is.na(cdadat$zTau))
min(cdadat$zTau, na.rm = T)
sum(is.na(cdadat$zTau_old)) # explain in the article why needed to be done
# exclude NA value from the vector

#decomposing hierarchy

get_T <- function(variabletbl){
  # variabletbl is column name for cda and its unique id, e.g., tbl$wd_uniqueID_2
  # var is variable on basis of which EI is calculated - needs to be added as attribute to palnets
  # tbl has to be open
  uc = unique(variabletbl)
  finres <- list()
  for (i in 1:length(uc)){
    c1 = tbl[variabletbl == uc[i], ] # argument 1
    nodeinc <-c1$id_1
    sid <- c1$wd_school_id_3[1]
    exg2 <- igraph::induced.subgraph(palnet[[sid]], which(V(palnet[[sid]])$name %in% nodeinc))
    nnet <- gorder(exg2)
    if((gsize(exg2) >2)) {
      
      # transitivity
      den2 = igraph::edge_density(exg2)
      netran <- igraph::transitivity(exg2, type = "global")
      
      #  Producing a tau statistic
      tc_t <- triad_census(exg2)
      
      weighting_scheme <- c(0,0,0,1,1,-1,0,0,1,0,0,1,1,0,0,0) # for tau
      
      #  Producing a tau statistic
      tau_r <- c()
      outd <- igraph::degree(exg2,mode = "out")
      ind <- igraph::degree(exg2,mode = "in")
      for (k in 1:250) {
        # Generate random graphs with a given degree sequence
        confg <- sample_degseq(out.deg = outd, in.deg = ind, 
                               method = "simple.no.multiple" )
        tau_c <- sum(triad_census(confg) * weighting_scheme)
        tau_r[k] <- tau_c
        
      }
      tau_t <- sum(tc_t * weighting_scheme)
      # tau score (normalized with configuration models)
      tau_er <- mean(tau_r)
      sd_tau_er <- round(sd(tau_r), 3) 
      zTau <-round((tau_t - tau_er)/sd_tau_er, 2)
      
    }
    else{
      tau_t <- NA
      zTau <- NA
    }
    
    finres[[i]] <- c(sid, uc[i],  round(tau_t, 3), round(zTau, 3))
  }
  
  frdat <- as.data.frame(do.call(rbind, finres))
  colnames(frdat) <- c("school","uniqueID","Tau", "zTau")
  return(frdat)
}

# new simpler hierarchy measure
get_Tnew <- function(variabletbl){
  # variabletbl is column name for cda and its unique id, e.g., tbl$wd_uniqueID_2
  # var is variable on basis of which EI is calculated - needs to be added as attribute to palnets
  # tbl has to be open
  uc = unique(variabletbl)
  finres <- list()
  for (i in 1:length(uc)){
    c1 = tbl[variabletbl == uc[i], ] # argument 1
    nodeinc <-c1$id_1
    sid <- c1$wd_school_id_3[1]
    exg2 <- igraph::induced.subgraph(palnet[[sid]], which(V(palnet[[sid]])$name %in% nodeinc))
    nnet <- gorder(exg2)
    if((gsize(exg2) >2)) {
      
      #  Producing a tau statistic
      weighting_scheme <- c(0,0,0,1,1,-1,0,0,1,0,0,1,1,0,0,0) # for tau
      tc_t <- triad_census(exg2) 
      tau_t <- sum(tc_t * weighting_scheme)
      stc <- sum(triad_census(exg2)) # divide with number of triads
      Tau_avg <- tau_t/stc
    }
    else{
      Tau_avg <- NA
    }
    
    finres[[i]] <- c(sid, uc[i],  round(Tau_avg, 2))
  }
  
  frdat <- as.data.frame(do.call(rbind, finres))
  colnames(frdat) <- c("school","uniqueID","Tau_avg")
  return(frdat)
}

# calculate for all
dataseh <- list()
for(i in 1:length(allcdas)){
  gei <- get_Tnew(allcdas[[i]])
  gei$cda <- namecdas[[i]][1]
  dataseh[[i]] <-gei
}

finda <- as.data.frame(data.table::rbindlist(dataseh))
dim(finda) #  4474    3
head(finda)
write.xlsx(finda, "Hierarchy_new_com_level_data_10GDM.xlsx")

# new measure of Tau Z based on random networks
# new measure of centralization that looks at in-going ties only

get_Cnew <- function(variabletbl){
  # variabletbl is column name for cda and its unique id, e.g., tbl$wd_uniqueID_2
  # var is variable on basis of which EI is calculated - needs to be added as attribute to palnets
  # tbl has to be open
  uc = unique(variabletbl)
  finres <- list()
  for (i in 1:length(uc)){
    c1 = tbl[variabletbl == uc[i], ] # argument 1
    nodeinc <-c1$id_1
    sid <- c1$wd_school_id_3[1]
    exg2 <- igraph::induced.subgraph(palnet[[sid]], which(V(palnet[[sid]])$name %in% nodeinc))
    nnet <- gorder(exg2)
    if((gsize(exg2) >2)) {
      
      #  centralization
      cdn <- centr_degree(
        graph = exg2,
        mode = "in",
        loops = FALSE,
        normalized = T)
      cnt <- cdn$centralization
    }
    else{
      cnt <- NA
    }
    
    finres[[i]] <- c(sid, uc[i],  round(cnt, 3))
  }
  
  frdat <- as.data.frame(do.call(rbind, finres))
  colnames(frdat) <- c("school","uniqueID","cnt")
  return(frdat)
}

# calculate for all
datasec <- list()
for(i in 1:length(allcdas)){
  gei <- get_Cnew(allcdas[[i]])
  gei$cda <- namecdas[[i]][1]
  datasec[[i]] <-gei
}

finda <- as.data.frame(data.table::rbindlist(datasec))
dim(finda) #  4474    3
head(finda)
write.xlsx(finda, "Centralization_new_com_level_data_10GDM.xlsx")


newt <- get_Tnew(allcdas[[10]])
hist(newt$Tau_avg)
sum(is.na(newt$Tau_avg))
0/0


tt <-get_T(allcdas[[10]])

sum(is.na(tt$Tau))
sum(is.na(tt$zTau))

# look for one
c1 = tbl[allcdas[[10]] == comstolook[6], ] #comstolook[11]
nodeinc <-c1$id_1
length(nodeinc)
sc <- c1$wu_school_id_3[1]
bbb<- palnet[[sc]]
exg2 <- igraph::induced.subgraph(bbb, which(V(bbb)$name %in% nodeinc))

#  Producing a tau statistic
tc_t <- triad_census(exg2) 
weighting_scheme <- c(0,0,0,1,1,-1,0,0,1,0,0,1,1,0,0,0)
tau_t <- sum(tc_t * weighting_scheme)
stc <- sum(triad_census(exg2))
tau_t/stc

# for tau

#  Producing a tau statistic
tau_r <- c()
outd <- igraph::degree(exg2,mode = "out")
ind <- igraph::degree(exg2,mode = "in")

confg <- sample_degseq(out.deg = outd, in.deg = ind, 
                       method = "simple.no.multiple" )
config <-erdos.renyi.game(gorder(exg2), gsize(exg2),
                          type = "gnm", directed = TRUE, loops = FALSE)
plot(exg2)
tau_c <- sum(triad_census(confg) * weighting_scheme)
tau_c
tric <- list()
set.seed(119956)
for (k in 1:250) {
  # Generate random graphs with a given degree sequence
  #confg <- sample_degseq(out.deg = outd, in.deg = ind, 
  #                      method = "simple.no.multiple" )
  config <-erdos.renyi.game(gorder(exg2), edge_density(exg2, loops = FALSE),
                            type = "gnp", directed = TRUE, loops = FALSE)
  config <- igraph::simplify(config, remove.multiple = TRUE)
  tau_c <- sum(triad_census(confg) * weighting_scheme)
  tau_r[k] <- tau_c
  tric[[k]] <- triad_census(confg)
  
}
tau_r
gsize(exg2)

tric
tau_t <- sum(tc_t * weighting_scheme)
# tau score (normalized with configuration models)
tau_er <- mean(tau_r)
sd_tau_er <- round(sd(tau_r), 3) 
zTau <-round((tau_t - tau_er)/sd_tau_er, 2)


sid <- c1$wd_school_id_3[1]
vec <- V(exg2)[var]
g2 <- length(which(vec==2))
g1 <- length(which(vec==1))
gna <- length(which(is.na(vec)))
cs <- gorder(exg2)
ex <- igraph::delete.vertices(exg2, V(exg2)[is.na(gender)])
as_gen <- assortativity_nominal(ex, vertex_attr(ex)$gender, directed=T)

plot(ex, vertex.size = 3, vertex.name = "")



# now for all
start_time <- Sys.time()
datasetran <- list()
for(i in 1:length(allcdas)){
  gei <- get_trans(allcdas[[i]])
  gei$cda <- namecdas[[i]][1]
  datasetran[[i]] <-gei
}
end_time <- Sys.time()
end_time - start_time # Time difference of 30.43271 mins

finda <- as.data.frame(data.table::rbindlist(datasetran))
dim(finda) #  4474    3
head(finda)
write.xlsx(finda, "Transit_corr_com_level_data_10GDM.xlsx")

# how many communities per GDM?
table(finda$cda)
sum(is.na(finda$Transit)) # 1269 - while 1178 size 2 or 1
nrow(finda)

min(finda$Transit, na.rm = T)
max(finda$Transit, na.rm = T)
histogram(finda$Transit, na.rm = T)
# some checkings

comdat3 <- as.data.frame(read_excel("Transit_corr_com_level_data_10GDM.xlsx"))
comdat3 = comdat3[comdat2$cda == "wd", ]

comdat3$cda <- NULL

# merge the two

comdatf <- merge(comdat1, comdat2, by = "uniqueID_2")

# mapping
comid <- comdat3$uniqueID_2
Transit <- comdat3$Transit


tbl$Transit <- mapvalues(tbl$wd_uniqueID_2,
                         from = comid, to = Transit)



cdn <- centr_degree(
  graph = exg2,
  mode = "all",
  loops = FALSE,
  normalized = F) #,


cdn <- centr_degree(
  graph = exg2,
  mode = "all",
  loops = FALSE,
  normalized = T)


uc = unique(tbl$wd_uniqueID_2)

c1 = tbl[tbl$wd_uniqueID_2 == "wt", ] 



igraph::transitivity(exg2, type = "global") 
triad_census(exg2)
weighting_scheme <- c(0,0,0,1,1,-1,0,0,1,0,0,1,1,0,0,0) 
sum(triad_census(exg2) * weighting_scheme)

sid <- c1$wd_school_id_3[1]
vec <- V(exg2)[var]
g2 <- length(which(vec==2))
g1 <- length(which(vec==1))
gna <- length(which(is.na(vec)))
cs <- gorder(exg2)
ei_gen = ei(exg2, "gender") # argument 2
or_gen = orwg(ex, "gender")
ex <- igraph::delete.vertices(exg2, V(exg2)[is.na(gender)])
as_gen <- assortativity_nominal(ex, vertex_attr(ex)$gender, directed=T)

plot(ex, vertex.size = 3, vertex.name = "")

# the datasets joined and saved under "Com_level_data_allCDA_1.xlsx"
