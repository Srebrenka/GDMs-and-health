###################
#     PaLS 04     #
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

# assortativity and E-I measures

#########################
#                       #
#    Load packages      #
#                       #
#########################

library(foreign)
library(Hmisc)
library(foreign)
library(dplyr)
library(statnet)
library(furniture)
library(igraph)
library(magrittr)
library(intergraph)
library(lubridate)

#########################
#                       #
#     Main Script       #
#                       #
#########################


source("Pals_02_create_nets.R") 



################ make binary variables for calculating assortativity and EI


for (i in schools$school){
  #Pull out individual schools#
  school = pals2 %>% 
    filter(schoolid == i)
  school$name <- as.character(school$palsqid)
  tomatch <- school[school$name %in% V(palnet[[i]])$name,]
  tomatch <- tomatch[match(V(palnet[[i]])$name, tomatch$name), ]
  
  V(palnet[[i]])$eth <- tomatch$ethnicity
  V(palnet[[i]])$gender <- tomatch$SEX
  V(palnet[[i]])$cortisol <- tomatch$LOGCORT1
  V(palnet[[i]])$BMI <-tomatch$BMI # why is this not working????
  V(palnet[[i]])$physatt <- tomatch$PHYSATT
  V(palnet[[i]])$year  <- tomatch$DOBYR  
  V(palnet[[i]])$fasraw  <- tomatch$FASRAW 
  V(palnet[[i]])$fas  <- tomatch$FAS  
  V(palnet[[i]])$smoke  <- tomatch$SMOKENOW  
  V(palnet[[i]])$alc  <- tomatch$ALCFREQ  
  V(palnet[[i]])$drug.effects  <- tomatch$drug.effects  
  V(palnet[[i]])$drugs1  <- tomatch$used.drugs
  V(palnet[[i]])$drugs2  <- tomatch$used.drugs2  
  V(palnet[[i]])$selfest  <- tomatch$SELFEST  
  V(palnet[[i]])$worries  <- tomatch$wr.raw.t  
  V(palnet[[i]])$physmat  <- tomatch$PHYSMAT  
  V(palnet[[i]])$ukcat  <- tomatch$UK90CATS  
  V(palnet[[i]])$cortab  <- tomatch$CORTAB  
  V(palnet[[i]])$corttb  <- tomatch$CORTTB  
  V(palnet[[i]])$logcort2  <- tomatch$LOGCORT2  
  V(palnet[[i]])$tchgeton  <- tomatch$TCHGETON  
  V(palnet[[i]])$arguepar  <- tomatch$ARGUEPAR  
  V(palnet[[i]])$pbicont  <- tomatch$PBICONT  
  V(palnet[[i]])$pbicare  <- tomatch$PBICARE 
  #binary
  V(palnet[[i]])$useddrugsB<- tomatch$used.drugs
  V(palnet[[i]])$smokB<- tomatch$smokingB
  V(palnet[[i]])$seB<- tomatch$selfestB
  V(palnet[[i]])$worB<- tomatch$worriesB
  V(palnet[[i]])$GHQB<- tomatch$ghqB
  
}



ex <- palnet[[1]]
ex <- igraph::delete.vertices(ex, V(ex)[is.na(gender)])
as_genB <- assortativity(ex, vertex_attr(ex)$gender, directed=T) #  0.8743218
as_gen <- assortativity_nominal(ex, vertex_attr(ex)$gender, directed=T) # 0.8738871
ei(ex, "gender") # -0.8756477


##############################################
# some additional/optional
##############################################


############# CALCULATE ASSORTATIVITY FOR EVERY TRAIT ####################
# LIKE EI INDEX FOR CONTINUOUS TRAITS OR ASSORTATIVITY_NOMINAL
# https://igraph.org/r/html/1.2.6/assortativity.html

# calculate ei and assortativity_nominal
# assortativity degree for in and out separately
ascf <- list()
for (i in 1:22) {
  as_deg <- assortativity_degree(palnet[[i]], directed = TRUE)
  as_deg.out <- assortativity(palnet[[i]], 
                              (igraph::degree(palnet[i][[1]], mode = 'out', loops = F, normalized = F)),
                              directed = T)
  as_deg.in <- assortativity(palnet[[i]], 
                             (igraph::degree(palnet[i][[1]], mode = 'in', loops = F, normalized = F)),
                             directed = T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(BMI)])
  as_bmi <- assortativity(ex, vertex_attr(ex)$BMI, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(gender)])
  as_gen <- assortativity_nominal(ex, vertex_attr(ex)$gender, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(cortisol)])
  as_cor <- assortativity(ex, vertex_attr(ex)$cortisol, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(physatt)])
  as_pa <- assortativity(ex, vertex_attr(ex)$physatt, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(year)])
  as_y <- assortativity(ex, vertex_attr(ex)$year, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(fasraw)])
  as_fasr <- assortativity(ex, vertex_attr(ex)$fasraw, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(fas)])
  as_fas <- assortativity(ex, vertex_attr(ex)$fas, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(smoke)])
  as_smok.c <- assortativity(ex, vertex_attr(ex)$smoke, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(smoke)])
  as_smok.b <- assortativity_nominal(ex, as.factor(vertex_attr(ex)$smoke), directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(alc)])
  as_alc <- assortativity(ex, vertex_attr(ex)$alc, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(drug.effects)])
  as_dr.ef <- assortativity(ex, vertex_attr(ex)$drug.effects, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(drugs1)])
  as_dr1.c <- assortativity(ex, vertex_attr(ex)$drugs1, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(drugs1)])
  as_dr1.b <- assortativity_nominal(ex, as.factor(vertex_attr(ex)$drugs1), directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(drugs2)])
  as_dr2 <- assortativity(ex, vertex_attr(ex)$drugs2, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(selfest)])
  as_selfest <- assortativity(ex, vertex_attr(ex)$selfest, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(worries)])
  as_worr <- assortativity(ex, vertex_attr(ex)$worries, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(physmat)])
  as_phymat <- assortativity(ex, vertex_attr(ex)$physmat, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(ukcat)])
  as_ukcat <- assortativity(ex, vertex_attr(ex)$ukcat, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(cortab)])
  as_cortab <- assortativity(ex, vertex_attr(ex)$cortab, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(corttb)])
  as_corttb <- assortativity(ex, vertex_attr(ex)$corttb, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(logcort2)])
  as_logcort2 <- assortativity(ex, vertex_attr(ex)$logcort2, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(tchgeton)])
  as_tch <- assortativity(ex, vertex_attr(ex)$tchgeton, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(arguepar)])
  as_argue <- assortativity(ex, vertex_attr(ex)$arguepar, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(pbicont)])
  as_pcont <- assortativity(ex, vertex_attr(ex)$pbicont, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(pbicare)])
  as_pcare <- assortativity(ex, vertex_attr(ex)$pbicare, directed=T)
  
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(eth)])
  as_eth.b <- assortativity_nominal(ex, as.factor(vertex_attr(ex)$eth), directed=T)
  
  asres <-   list(as_deg, as_deg.out, as_deg.in, as_bmi, as_gen, as_cor, as_pa, as_y, as_fasr,
                  as_fas, as_smok.c, as_smok.b, as_alc, as_dr.ef, as_dr1.c, as_dr1.b, as_dr2,
                  as_selfest, as_worr, as_phymat, as_cortab, as_corttb, as_logcort2, as_tch,
                  as_argue, as_pcont, as_pcare, as_eth.b)
  
  ascf[[length(ascf) + 1]] = asres
}


sqa <- matrix(unlist(ascf), ncol = 28, byrow = TRUE)
asdat <- as.data.frame(sqa)

colnames(asdat) = c("as_deg", "as_deg.out", "as_deg.in", "as_bmi", "as_gen", "as_cor", "as_pa", "as_y", "as_fasr",
                    "as_fas", "as_smok.c", "as_smok.b", "as_alc", "as_dr.ef", "as_dr1.c", "as_dr1.b", "as_dr2",
                    "as_selfest", "as_worr", "as_phymat", "as_cortab", "as_corttb", "as_logcort2", "as_tch",
                    "as_argue", "as_pcont", "as_pcare", "as_eth.b")

write.table(asdat, file = "assortativity_PaLS1310.csv", dec = ',' ,
            sep = ";",row.names = TRUE)

# EI for binary variables
library(netseg)

ascfB <- list()
for (i in 1:22) {
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(gender)])
  ei_gen = ei(ex, "gender")
  ex <- palnet[[i]]
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(useddrugsB)])
  ei_drugs = ei(ex, "useddrugsB")
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(smokB)])
  ei_smok = ei(ex, "smokB")
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(seB)])
  ei_selfe = ei(ex, "seB")
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(worB)])
  ei_wor = ei(ex, "worB")
  ex <- igraph::delete.vertices(ex, V(ex)[is.na(GHQB)])
  ei_ghq = ei(ex, "GHQB")
  
  asresB<- list(ei_gen, ei_drugs, ei_smok, ei_selfe,
                ei_wor, ei_ghq)
  ascfB[[length(ascfB) + 1]] = asresB
}

sqa <- matrix(unlist(ascfB), ncol = 6, byrow = TRUE)
asdat <- as.data.frame(sqa)

colnames(asdat) = c("ei_gen", "ei_drugs", "ei_smok", "ei_selfe",
                    "ei_wor", "ei_ghq")

write.table(asdat, file = "EI_PaLS0511.csv", dec = ',' ,
            sep = ";",row.names = TRUE)

