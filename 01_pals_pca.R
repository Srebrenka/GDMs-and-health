###################
#     01 PaLS     #
###################

rm(list = ls())

#################
#               #
#      Notes    #
#               #
#################

# Srebrenka Letina developed this script


#############
#  Purpose  #
#############

# read and prepare Pals data
# create PCA scores
# create dataset with attributes for future analyses

#########################
#                       #
#    Load packages      #
#                       #
#########################

library(data.table)
library(readxl)
library(openxlsx)
library(dplyr)
library(haven)          # for reading spss
#library(foreign)
library(BBmisc)         # for normalization
library(psych)          # pca and describe
library(FactoMineR)     # pca
library(factoextra)     # pca
library(nortest)        # testing normality of distributions
#########################
#                       #
#     Main script       #
#                       #
#########################
#setwd()

# reading the dataset
dataset = read_sav("PalsData.SAV") # not avaialble without request

# Overview of the data
dim(dataset)
# 3194  438


# replacing NAs
dataset[dataset == -3] <- NA
dataset[dataset == 999] <- NA
# Worries variables
wrdat <- subset(dataset, select = c(WRSCHOOL,
                                    WRUNEMPL,
                                    WRWORK, 
                                    WRFAMGET, 
                                    WRROMANT,
                                    WRHEALTH,
                                    WROVERWT,
                                    WRTHIN,
                                    WRPREG,
                                    WRLOOKS))

psych::alpha(wrdat) # 0.78
# total score
dataset$wr.raw.t <- rowSums(wrdat, na.rm = FALSE)
hist(dataset$wr.raw.t)


# creating drug.effects scores
drugitems<- subset(dataset, select = c(DRGMORE,
                                       DRGFORGT,
                                       DRGTROUB))
psych::alpha(drugitems) # 0.72

dataset$drug.effects <- dataset$DRGMORE + dataset$DRGFORGT + 
  dataset$DRGTROUB

# creating using.drugs scores
drugset <- subset(dataset, select = c(CANNABIS,VALIUM, AMPHETAM, LSD, ECSTASY, 
                                      SOLVENTS, COCAINE, HEROIN, MAGICMUS))
drugset2 <- drugset

# 1 - ever used; 0 - never
drugset2 <- ifelse(is.na(drugset2), NA,
                   ifelse(drugset2 <= 4, 1, 0))
drugset2$drugs <- rowSums(drugset2, na.rm = TRUE)

# binary score: 1 - ever used any drug
drugset2$drugs <- ifelse(is.na(drugset2$drugs), NA,
                         ifelse(drugset2$drugs >=1, 1, 0))
dataset$used.drugs <- drugset2$drugs 

# three categories score 
drugset3 <- as.data.frame(drugset)
# 0 - never used
# 1 - less often or not past year
# 2 - every day or weekly

drugset3 <- ifelse(is.na(drugset3), NA,
                   ifelse(drugset3 <= 2, 2,
                          ifelse(drugset3 ==3 | drugset3 == 4, 1, 0)))
drugset3 <- as.data.frame(drugset3)
# any drug used 0, 1, or 2

drugset3$drugs<-pmax(drugset3$CANNABIS,drugset3$VALIUM, drugset3$AMPHETAM, 
                     drugset3$LSD, drugset3$ECSTASY, drugset3$SOLVENTS, 
                     drugset3$COCAINE, drugset3$HEROIN, drugset3$MAGICMUS)

# to control for NA
drugset3$drugs<-do.call(pmax, c(drugset3[1:9], list(na.rm=TRUE)))

dataset$used.drugs2 <- drugset3$drugs

# subseting to relevant variables
atdat <-  subset(dataset, select = c(SCHOOLID,
                                     PALSQID, # this is ID
                                     SEX,
                                     DOBYR, 
                                     FASRAW, 
                                     FAS, # socio
                                     SMOKENOW,
                                     ALCFREQ,
                                     drug.effects,
                                     used.drugs,
                                     used.drugs2,
                                     SELFEST,
                                     GHQLIKS,
                                     wr.raw.t,
                                     PBICONT,
                                     PBICARE))


atdat$gender <- ifelse(is.na(atdat$SEX), NA, # 1 is female
                       ifelse(atdat$SEX == 1, 0, 1))
table(atdat$gender)
atdat$SEX <- NULL
atdat$ethnicity <- ifelse(is.na(atdat$OWNETH), NA,
                          ifelse(atdat$OWNETH == 1, 1, 0))

table(atdat$ethnicity)


# age in months in dataset
dataset$AGEMONTH
dataset$age <- round(dataset$AGEMONTH/12, 0)
atdat$age <- dataset$age 
atdat$DOBYR <- NULL
# alpha for other scales

# GHQ-12
gsub <- subset(dataset, select = c(GHQ1, GHQ2, GHQ3, GHQ4, GHQ5, GHQ6, GHQ7, 
                                   GHQ8, GHQ9, GHQ10, GHQ11, GHQ12))
psych::alpha(gsub) # 0.85
# FASRAW
fsub <- subset(dataset, select = c(CARVAN, OWNROOM, NCOMPUTE, HOLIDAYS))
psych::alpha(fsub, check.keys = TRUE) # 0.49

# self-esteem scale 
ssub <- subset(dataset, select = c(SESURE, SENOTME, SELIKE, SELOWOP, SEFAIL, 
                                   SECHANGE, SEDOWELL, SESATISF, SELIKEME, SEGOODQL))
psych::alpha(ssub, check.keys = TRUE) # 0.86

# pbicare
csub <- subset(dataset, select = c(PBHELP, PBLOVE, PBUNDERS, PBBETTER))
psych::alpha(csub, check.keys = TRUE) # 0.72

# pbicont
ctsub <- subset(dataset, select = c(PBLET, PBOWNDEC, PBCONTROL, PBBABY))
psych::alpha(ctsub, check.keys = TRUE) # 0.59



#-----------------------------------------------------------------------------
#               PRINCIPAL COMPONENT ANALYSIS 
#-----------------------------------------------------------------------------
############################################################################
#             extracting PCs and creating component scores
#                    and saving it in a dataset
###########################################################################

# 7 variables selected for PCA

datv <-  subset(atdat, select = c(SMOKENOW,
                                  ALCFREQ,
                                  #used.drugs,
                                  used.drugs2,
                                  drug.effects,
                                  SELFEST,
                                  GHQLIKS,
                                  wr.raw.t))
#--------------------descriptives of variables

desvars <- psych:: describe(datv, skew = TRUE, ranges = TRUE)
desvars$variables <- rownames(desvars)
# save
write.xlsx(desvars, "descriptive of variables for PCA.xlsx")


datvn <- datv
colnames(datv)

# higher score means:
# worries - less worries (make it R)
# alcfreq - less frequently drinking (make it R)
# smokenow - smoking more
# used.drugs2 - more using drugs
# drug.effects - less effects (make it R)
# selfest - higher selfesteem (make it R)
# ghq - worse mental health


#-------------------------------------------------------
#                   standardize variables
#-------------------------------------------------------
datvnN <- BBmisc::normalize(datvn, method = "standardize")

cortest.bartlett(datvnN) # p = 0
KMO(datvnN) # Overall MSA =  0.71

# after standardization
colnames(datvnN) <- c("smoking", "drinking", "using_drugs", 
                      "drug_effects", "low_self_esteem",
                      "GHQ_scores_neg", "worries")
#Turn them all to negative
# self-esteem, worries, drug-effects and drinking - recode to neg. outcome
datvnN$drinking <- -1*datvnN$drinking
datvnN$worries <- -1*datvnN$worries
datvnN$low_self_esteem <- -1*datvnN$low_self_esteem 
datvnN$drug_effects <- -1*datvnN$drug_effects
#create_report(datvn)

# from psych package - extracting two pca scores

p1a_st <- principal(datvnN, nfactors=2, scores = TRUE, missing = TRUE,
                    impute = "median", rotate = "none")
p1a_st$values
p1a_st$weights

psych::alpha(datvnN, check.keys=TRUE) #  0.69  
# results are the same as for raw data, but the scores are not

library(hornpa)
hornpa(k = 7, size = nrow(datvn), reps = 500, seed = 1234) 

fs <- factor.scores(x = datvnN, f = p1a_st)

pcaweights <- as.data.frame(fs$weights)
write.xlsx(pcaweights, "PCA_weights.xlsx", rowNames = T)

fs1<- as.data.frame(fs$scores)
dim(fs1)
dim(atdat)
colnames(atdat)
head(fs1)
atdat2<- cbind(atdat, fs1)
colnames(atdat2)
atdat3 <- subset(atdat2, select = c(PALSQID,
                                    PC1,
                                    PC2))
final <- cbind(atdat3, datvnN) # all turned to neagtive!!!
colnames(final)
# if any value is NA sum is NA
final$rs1 <- final$smoking + final$drinking + final$using_drugs + final$drug_effects
final$rs2 <- final$low_self_esteem + final$low_self_esteem + final$worries

write.xlsx(final, "ids_and_pc_and_rs_turned_neg.xlsx")


datv2<- normalize(datv, method = "standardize")
datv2R <- datv2
# creating raw composite scores
# first make all variable in the same direction - higher score - worse outcome
datv2R$wr.raw.t <- -1 * datv2R$wr.raw.t 
datv2R$ALCFREQ <- -1 * datv2R$ALCFREQ
datv2R$drug.effects <- -1 * datv2R$drug.effects
datv2R$SELFEST <- -1 * datv2R$SELFEST

datv2R$com1 <- (datv2R$SMOKENOW + datv2R$ALCFREQ + datv2R$used.drugs2 + datv2R$drug.effects)/4
datv2R$com2 <- (datv2R$SELFEST + datv2R$GHQLIKS + datv2R$wr.raw.t)/3

# adding them to atdat2
atdat2$com1 <- datv2R$com1
atdat2$com2 <- datv2R$com2

datv2R$PC1 <- atdat2$PC1
datv2R$PC2 <- atdat2$PC2




atdat2 <- subset(atdat, select = c(SCHOOLID,
                                   PALSQID,
                                   gender,
                                   age,
                                   ethnicity,
                                   FASRAW,
                                   FAS,
                                   TCHGETON,
                                   ARGUEPAR,
                                   PBICONT,
                                   PBICARE,
                                   PC1,
                                   PC2,
                                   rcs1,
                                   rcs2))

dim(atdat2) # 3194   15
write.table(atdat, file = "attribute_data_pals2_FIN.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
# save attribute data
write.xlsx(atdat, "attribute_data_pals2_FIN.xlsx")

