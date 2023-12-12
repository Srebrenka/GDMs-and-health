###################
#     PaLS 14     #
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

# multi-level models
# some additional descriptives

# NOTE 1: names of variables are old names, later changed to be more suitable
# NOTE 2: This is for one GDM (WT: wlaktrap); and the script has to run for all GDMs

#########################
#                       #
#    Load packages      #
#                       #
#########################


library(ggplot2)  # Graphs
library(lme4)     # Estimating multilevel models 
library(lmerTest) # Printing p-values for the regression coefficients
library(lmtest)   # Likelihood ratio tests
library(optimx)  # required for some optimizer
library(openxlsx)
library(readxl)
library(car)
library(performance)
library(tidyr) # pipe operator
library(MuMIn)
library(stevemisc)
library(dplyr)
library(plyr)
library(texreg)
library(stringr)
library(tidyr)
library(purrr)
# optional:
library(merTools) # good plots
library(see)
library(corrplot)
library(MuMIn)
library(GGally)
library(ggeffects)
library(EMAtools)


###############################
# for each GDM 
##############################

# reading in the data for the GDM and level 1 data
cdadat <- as.data.frame(read_excel("GDMdata.xlsx", # use appropriate acronym instead GDM
                                   col_names = TRUE))

###################################
#  Descriptives for properties
####################################

# univariate descriptives
# delete all rows that have the same value in a column
df = cdadat[!duplicated(cdadat$cda_uniqueID_2),]

# get rid of L1 variables
df2<- subset(df, select= c(cda_uniqueID_2, PC1_tcs, PC2, 
                           cda_com_size_2, typegencom3, 
                           cda_rat_2, Transit, CntN, Tau))
# turn typegencom3 to numeric
df2$typegencom3 <- ifelse(df2$typegencom3 == "mixed", 0,
                          ifelse(df2$typegencom3 == "female", 1, -1))
# delete IDs
df2$cda_uniqueID_2 <- NULL
colnames(df2) <- c("SU", "MW", "Com.Size", "GCn",
                   "ROTC","Transit.", "Centr.", "Hierar.")
# corrplot
nrow(df2)
pdf(file = "Corrplot_DV_and_CPs_WT_level2.pdf",
    width =4.5, height = 4.5)   # The directory you want to save the file in
testRes = cor.mtest(df2,  use="pairwise.complete.obs", conf.level = 0.95)
M = cor(df2, use="pairwise.complete.obs", method = "pearson")
corrplot.mixed( M, p.mat = testRes$p,             
                lower = "number", 
                upper = "circle",
                tl.col = "black",
                number.cex = 0.75,
                insig =  "blank", 
                addCoef.col = 'black',
                tl.cex = 0.6,
                tl.pos = "lt")
dev.off()

# descriptives of com properties
# leave just numeric CPs (all except GCn)
df3 <- subset(df2, select = -c(SU, MW, GCn))
colnames(df3)

cda2des<- as.data.frame(psych::describe(df3, na.rm = TRUE))
cda2des <- round(cda2des, 2)
cda2des$vars <- rownames(cda2des)
write.xlsx(cda2des, "descriptive_level2_5comPro_WT.xlsx")

# GCn
raw <- table(df2$GCn)
prct <- round(prop.table(table(df2$GCn))*100, 1)

gcdes <- as.data.frame(cbind(raw, prct))
write.xlsx(gcdes, "descriptive_level2_gen_comp_WT.xlsx", rowNames=TRUE)

# distributions of numeric values

df3 %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()

# scatterplots

options(repr.plot.width = 20, repr.plot.height = 10)
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(size=0.8, alpha= 0.4) + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}
colnames(df2)
pdf(file = "Scatterplot_Level2_WT.pdf") #  height = 4, width = 5
# change plot size (optional)
options(repr.plot.width = 20, repr.plot.height = 10, digits =4)
g = ggpairs(df2,columns = 3:8, lower = list(continuous = my_fn))
g+ theme(axis.text = element_text(size = 5))
dev.off()

# level 1 descriptive - bivariate

subdat<- subset(cdadat, select = c(PC1_tcs,
                                   PC2, 
                                   cda_com_size_2,
                                   typegencom3, 
                                   cda_rat_2,
                                   Transit,
                                   CntN,
                                   Tau))
sum(is.na(subdat$typegencom3)) # 0
subdat$typegencom3 <- ifelse(subdat$typegencom3 == "mixed", 0,
                             ifelse(subdat$typegencom3 == "female", 1, -1))
table(subdat$typegencom3)
colnames(subdat) <- c("SU", "MW", "Com.Size", "GCn", 
                      "ROTC", "Transit.", "Centr.", "Hierar.")
# correlations between measures - level 1
# WT correlations


pdf(file = "Corrplot_DV_and_CPs_WT_level1.pdf",
    width =4.5, height = 4.5)   # The directory you want to save the file in
testRes = cor.mtest(subdat,  use="pairwise.complete.obs", conf.level = 0.95)
M = cor(subdat, use="pairwise.complete.obs", method = "pearson")
corrplot.mixed( M, p.mat = testRes$p,             
                lower = "number", 
                upper = "circle",
                tl.col = "black",
                number.cex = 0.75,
                insig =  "blank", 
                addCoef.col = 'black',
                tl.cex = 0.6,
                tl.pos = "lt")
dev.off()


####################################
#   M1, M1.1, M2,  M3, M4
###################################

#pc1 - SU

M1 <-lmer(PC1_tcs ~ 1+ (1|cda_uniqueID_2), # 0.11
          data   = cdadat,
          REML   = TRUE,
          control = lmerControl(optimizer ="Nelder_Mead"))  

M1_1 <-lmer(PC1_tcs ~ 1+  (1|cda_uniqueID_2) + (1|school_id_3) , 
            data   = cdadat,
            REML   = TRUE,
            control = lmerControl(optimizer ="Nelder_Mead"))

M2 <-lmer(PC1_tcs ~ gender + age_dich + ethnicity + Family_affluence + 
            Parent_control_cs + Parent_care_cs +
            (1|cda_uniqueID_2), 
          data   = cdadat,
          REML   = TRUE,
          control = lmerControl(optimizer ="Nelder_Mead"))

M3 <-lmer(PC1_tcs ~ gender + age_dich + ethnicity + Family_affluence + 
            Parent_control_cs + Parent_care_cs +
            cda_com_size_2cs + cda_rat_2 + typegencom3 + 
            Transit + CntN + Tau + 
            (1|cda_uniqueID_2), 
          data   = cdadat,
          REML   = TRUE,
          control = lmerControl(optimizer ="Nelder_Mead"))

M4 <- lmer(PC1_tcs ~ gender + age_dich + ethnicity + Family_affluence + 
             Parent_control_cs + Parent_care_cs +
             cda_com_size_2cs + cda_rat_2 + typegencom3 + 
             Transit + CntN + Tau +
             net_size_3cat + 
             cda_Modularity_3 + per_F_n1_3 +
             (1|cda_uniqueID_2), 
           data   = cdadat,
           REML   = TRUE,
           control = lmerControl(optimizer ="Nelder_Mead"))

wordreg(list(M1,
             M1_1,
             M2,
             M3,
             M4),
        custom.model.names = c("M1",
                               "M1_1",
                               "M2",
                               "M3",
                               "M4"),
        title            = "wd: Basic models pc1",
        caption.above      = TRUE,
        single.row         = TRUE, 
        file = "WT -Final table of basic models SU.doc")

# CIs - not for M2_1

test1 <- confint(M1, level = 0.95)
testy1 <- as.data.frame(test1)
testy1 <- round(testy1, 2)
colnames(testy1) <- paste("M1", colnames(testy1), sep = "_")
testy1$parameter <- row.names(testy1)

test1.1 <- confint(M1_1, level = 0.95)
testy1.1 <- as.data.frame(test1.1)
testy1.1 <- round(testy1.1, 2)
colnames(testy1.1) <- paste("M1_1", colnames(testy1.1), sep = "_")
testy1.1$parameter <- row.names(testy1.1)

test2 <- confint(M2, level = 0.95)
testy2 <- as.data.frame(test2)
testy2 <- round(testy2, 2)
colnames(testy2) <- paste("M2", colnames(testy2), sep = "_")
testy2$parameter <- row.names(testy2)

test3 <- confint(M3, level = 0.95)
testy3 <- as.data.frame(test3)
testy3 <- round(testy3, 2)
colnames(testy3) <- paste("M3", colnames(testy3), sep = "_")
testy3$parameter <- row.names(testy3)

test4 <- confint(M4, level = 0.95)
testy4 <- as.data.frame(test4)
testy4 <- round(testy4, 2)
colnames(testy4) <- paste("M4", colnames(testy4), sep = "_")
testy4$parameter <- row.names(testy4)

cis <- list(testy1, testy1.1, testy2, testy3, testy4)
cis[[3]]

cisf <- Reduce(function(...) merge(..., by="parameter", all=T), cis)

write.xlsx(cisf, file = "CI_estimates_models_pc1.xlsx")


# take out the estimates
Models_SU<- list(M1, M1_1, M2, M2_1, M3, M4)

estima <- list() # put estimates in df

for(i in 1:length(Models_SU)){
  Ms <- summary(Models_SU[[i]])
  dat <- as.data.frame(Ms$coefficients)
  dat<- subset(dat, select = c(Estimate))
  dat$Estimate <- round(dat$Estimate, 2)
  dat$parameter <- rownames(dat)
  estima[[i]] <- dat
  
}

# combine dfs

for(i in 1:6){
  colnames(estima[[i]]) <-paste(colnames(estima[[i]]), i, sep="_")
  colnames(estima[[i]])[2]<- "parameter"
}

dfes1 <- Reduce(function(...) merge(..., by= "parameter", all=T), estima)
dfes1 # not ordered though

write.xlsx(dfes1, file = "estimates_SU_all_models_WT_newTau.xlsx")

#pc2 - MW

M1a <-lmer(PC2 ~ 1+ (1|cda_uniqueID_2), 
           data   = cdadat,
           REML   = TRUE,
           control = lmerControl(optimizer ="Nelder_Mead"))  

M1_1a <-lmer(PC2 ~ 1+  (1|cda_uniqueID_2) + (1|school_id_3) , 
             data   = cdadat,
             REML   = TRUE,
             control = lmerControl(optimizer ="Nelder_Mead"))

M2a <-lmer(PC2 ~ gender + age_dich + ethnicity + Family_affluence + 
             Parent_control_cs + Parent_care_cs +
             (1|cda_uniqueID_2), 
           data   = cdadat,
           REML   = TRUE,
           control = lmerControl(optimizer ="Nelder_Mead"))

M2_1a <-lmer(PC2 ~ gender + age_dich + ethnicity + Family_affluence + 
               Parent_control_cs + Parent_care_cs +
               cda_com_size_2cs + cda_rat_2 + typegencom3 + 
               (1|cda_uniqueID_2), 
             data   = cdadat,
             REML   = TRUE,
             control = lmerControl(optimizer ="Nelder_Mead"))

M3a <-lmer(PC2 ~ gender + age_dich + ethnicity + Family_affluence + 
             Parent_control_cs + Parent_care_cs +
             cda_com_size_2cs + cda_rat_2 + typegencom3 + 
             Transit + CntN + Tau + 
             (1|cda_uniqueID_2), 
           data   = cdadat,
           REML   = TRUE,
           control = lmerControl(optimizer ="Nelder_Mead"))

M4a <- lmer(PC2 ~ gender + age_dich + ethnicity + Family_affluence + 
              Parent_control_cs + Parent_care_cs +
              cda_com_size_2cs + cda_rat_2 + typegencom3 + 
              Transit + CntN + Tau +
              net_size_3cat + 
              cda_Modularity_3 + per_F_n1_3 +
              (1|cda_uniqueID_2), 
            data   = cdadat,
            REML   = TRUE,
            control = lmerControl(optimizer ="Nelder_Mead"))

wordreg(list(M1a,
             M1_1a,
             M2a,
             M2_1a,
             M3a,
             M4a),
        custom.model.names = c("M1",
                               "M1_1",
                               "M2",
                               "M2_2",
                               "M3",
                               "M4"),
        title            = "wd: Basic models pc2",
        caption.above      = TRUE,
        single.row         = TRUE, 
        file = "WT -Final table of basic models MW.doc")

# CIs - not for M2_1
test1 <- confint(M1a, level = 0.95)
testy1 <- as.data.frame(test1)
testy1 <- round(testy1, 2)
colnames(testy1) <- paste("M1", colnames(testy1), sep = "_")
testy1$parameter <- row.names(testy1)

test1.1 <- confint(M1_1a, level = 0.95)
testy1.1 <- as.data.frame(test1.1)
testy1.1 <- round(testy1.1, 2)
colnames(testy1.1) <- paste("M1_1", colnames(testy1.1), sep = "_")
testy1.1$parameter <- row.names(testy1.1)

test2 <- confint(M2a, level = 0.95)
testy2 <- as.data.frame(test2)
testy2 <- round(testy2, 2)
colnames(testy2) <- paste("M2", colnames(testy2), sep = "_")
testy2$parameter <- row.names(testy2)

test3 <- confint(M3a, level = 0.95)
testy3 <- as.data.frame(test3)
testy3 <- round(testy3, 2)
colnames(testy3) <- paste("M3", colnames(testy3), sep = "_")
testy3$parameter <- row.names(testy3)

test4 <- confint(M4a, level = 0.95)
testy4 <- as.data.frame(test4)
testy4 <- round(testy4, 2)
colnames(testy4) <- paste("M4", colnames(testy4), sep = "_")
testy4$parameter <- row.names(testy4)

cis <- list(testy1, testy1.1, testy2, testy3, testy4)


cisf <- Reduce(function(...) merge(..., by="parameter", all=T), cis)

write.xlsx(cisf, file = "CI_estimates_models_pc2.xlsx")

# take out the estimates
Models_MW<- list(M1a, M1_1a, M2a, M2_1a, M3a, M4a)

estima <- list() # put estimates in df

for(i in 1:length(Models_MW)){
  Ms <- summary(Models_MW[[i]])
  dat <- as.data.frame(Ms$coefficients)
  dat<- subset(dat, select = c(Estimate))
  dat$Estimate <- round(dat$Estimate, 2)
  dat$parameter <- rownames(dat)
  estima[[i]] <- dat
  
}

# combine dfs

for(i in 1:6){
  colnames(estima[[i]]) <-paste(colnames(estima[[i]]), i, sep="_")
  colnames(estima[[i]])[2]<- "parameter"
}

dfes1 <- Reduce(function(...) merge(..., by= "parameter", all=T), estima)
dfes1 # not ordered though

write.xlsx(dfes1, file = "estimates_MW_all_models_WT_newTau.xlsx")


# make tables for all except M2_1 
# only in SM with estimates only (saved in word doc)

CI1 <- as.data.frame(read_excel( "CI_estimates_models_pc1.xlsx", col_names = T))
CI1 <- CI1[-(1:3), , drop = FALSE]
es1 <- as.data.frame(read_excel( "estimates_SU_all_models_WT_newTau.xlsx", col_names = T))
ordf <- as.data.frame(read_excel( "orderParameters.xlsx", col_names = T))
es1 <- merge(es1, ordf, by = "parameter")
sues <- merge(es1, CI1, by = "parameter")
sues <- sues[order(sues$order),]
sues$order <- NULL
sues
sues$m1CI <- str_c(sues$Estimate_1," ","[",sues$`M1_2.5 %`, ", ", sues$`M1_97.5 %`,"]")
sues$m1.1CI <- str_c(sues$Estimate_2," ","[",sues$`M1_1_2.5 %` , ", ", sues$`M1_1_97.5 %`,"]")
sues$m2CI <- str_c(sues$Estimate_3," ","[",sues$`M2_2.5 %`, ", ", sues$`M2_97.5 %`,"]")
# skip Estimate_4, that is M2_1, no CIs
sues$m3CI <- str_c(sues$Estimate_5," ","[",sues$`M3_2.5 %`, ", ", sues$`M3_97.5 %`,"]")
sues$m4CI <- str_c(sues$Estimate_6," ","[",sues$`M4_2.5 %`, ", ", sues$`M4_97.5 %`,"]")
colnames(sues)
sues <- sues[, c(1, 18, 19, 20, 21, 22)]
head(sues,1)
write.xlsx(sues, file = "for table SU.xlsx")

# pc2

CI1 <- as.data.frame(read_excel( "CI_estimates_models_pc2.xlsx", col_names = T))
CI1 <- CI1[-(1:3), , drop = FALSE]
es1 <- as.data.frame(read_excel( "estimates_MW_all_models_WT_newTau.xlsx", col_names = T))
ordf <- as.data.frame(read_excel( "orderParameters.xlsx", col_names = T))
es1 <- merge(es1, ordf, by = "parameter")
mwes <- merge(es1, CI1, by = "parameter")
mwes <- mwes[order(mwes$order),]
mwes$order <- NULL

mwes$m1CI <- str_c(mwes$Estimate_1," ","[",mwes$`M1_2.5 %`, ", ", mwes$`M1_97.5 %`,"]")
mwes$m1.1CI <- str_c(mwes$Estimate_2," ","[",mwes$`M1_1_2.5 %` , ", ", mwes$`M1_1_97.5 %`,"]")
mwes$m2CI <- str_c(mwes$Estimate_3," ","[",mwes$`M2_2.5 %`, ", ", mwes$`M2_97.5 %`,"]")
# skip Estimate_4, that is M2_1, no CIs
mwes$m3CI <- str_c(mwes$Estimate_5," ","[",mwes$`M3_2.5 %`, ", ", mwes$`M3_97.5 %`,"]")
mwes$m4CI <- str_c(mwes$Estimate_6," ","[",mwes$`M4_2.5 %`, ", ", mwes$`M4_97.5 %`,"]")
colnames(mwes)
mwes <- mwes[, c(1, 18, 19, 20, 21, 22)]
head(mwes,1)
write.xlsx(mwes, file = "for table MW.xlsx")


# Rs and both ICC for models (for tables)


SUs <- list()
for(i in 1:4){
  R2m <- round(r.squaredGLMM(Models_SU[[i]])[1],2)
  R2c <- round(r.squaredGLMM(Models_SU[[i]])[2],2)
  AdIcc <- round(as.numeric(icc(Models_SU[[i]])[1]),2)
  CoIcc <- round(as.numeric(icc(Models_SU[[i]])[2]),2)
  res<- c(ModelName[i], R2m, R2c, AdIcc, CoIcc)
  SUs[[i]] <- res
}

SUdf <- as.data.frame(do.call("rbind", SUs))
colnames(SUdf) <- c("model","R2m","R2c", "AdIcc", "CoIcc")

write.xlsx(SUdf, file = "table of Rs and ICCs per model - SU.xlsx")

# pc2
MWs <- list()
for(i in 1:6){
  R2m <- round(r.squaredGLMM(Models_MW[[i]])[1],2)
  R2c <- round(r.squaredGLMM(Models_MW[[i]])[2],2)
  AdIcc <- round(as.numeric(icc(Models_MW[[i]])[1]),2)
  CoIcc <- round(as.numeric(icc(Models_MW[[i]])[2]),2)
  res<- c(ModelName[i], R2m, R2c, AdIcc, CoIcc)
  MWs[[i]] <- res
}

MWdf <- as.data.frame(do.call("rbind", MWs))
colnames(MWdf) <- c("model","R2m","R2c", "AdIcc", "CoIcc")

write.xlsx(MWdf, file = "table of Rs and ICCs per model - MW.xlsx")

# CI for ICCs

#new lists without M1_1 and M2_1 - 4 basic models

ModelName2 <- c("M1", "M2", "M3", "M4")
Models_MW2 <- list(M1a,  M2a,  M3a, M4a)
Models_SU2 <- list(M1,  M2,  M3, M4)

#Calculate bootstrap distribution - 250 simulations 

calc.icc <- function(y) {
  sumy <- summary(y)
  if(is.null(sumy$varcor$cda_uniqueID_2[1])){ return(NA)}
  else{(sumy$varcor$cda_uniqueID_2[1]) / (sumy$varcor$cda_uniqueID_2[1] + sumy$sigma^2)}
} # this is not correct for models with slope, or two levels M1_2

# SU
cis <-list()
for(i in 1:length(Models_SU2)){
  
  boot.icc <- bootMer(Models_SU2[[i]], seed = 1088, calc.icc,
                      nsim=250) 
  CI1 <- quantile(boot.icc$t, 0.025, na.rm = T)
  CI1 <- as.numeric(CI1)
  CI2 <- quantile(boot.icc$t, 0.975, na.rm = T)
  CI2 <- as.numeric(CI2)
  CI <- c(CI1, CI2)
  cis[[i]] <- CI
  print(i)
}

cidf <- as.data.frame(do.call("rbind", cis))
colnames(cidf) = c("CI1", "CI2")
cidf <- round(cidf, 2) 
row.names(cidf) = ModelName2

write.xlsx(cidf, "CIs for SU.xlsx", colNames=TRUE, rowNames=TRUE)


# MW
cis <-list()
for(i in 1:length(Models_MW2)){
  
  boot.icc <- bootMer(Models_MW2[[i]], seed = 1088, calc.icc,
                      nsim=250) 
  CI1 <- quantile(boot.icc$t, 0.025, na.rm = T)
  CI1 <- as.numeric(CI1)
  CI2 <- quantile(boot.icc$t, 0.975, na.rm = T)
  CI2 <- as.numeric(CI2)
  CI <- c(CI1, CI2)
  cis[[i]] <- CI
  print(i)
}

cidf <- as.data.frame(do.call("rbind", cis))
colnames(cidf) = c("CI1", "CI2")
cidf <- round(cidf, 2) 
row.names(cidf) = ModelName2

write.xlsx(cidf, "CIs for MW.xlsx", colNames=TRUE, rowNames=TRUE)

# figure

fidat <- as.data.frame(read_excel("for icc fig1.xlsx", col_names = TRUE))
fidat$Outcome <- as.factor(fidat$Outcome)
ggplot(fidat, aes(x=Models, y=ICC, fill=Outcome)) +
  geom_bar(position=position_dodge(.9), stat="identity") + # legend.position="bottom"
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=lCI, ymax=uCI)) +
  scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
  theme(legend.position="bottom")


########################
# compare performance
#######################

# pairwaise preformance

# SU models
mc1 <- anova(M1, M1_1, refit = F)
mc1 <- round(mc1$`Pr(>Chisq)`[2],3)

mc2 <- anova(M1, M2, refit = F)
mc2 <- round(mc2$`Pr(>Chisq)`[2],3)


# to fit at the same dataset
fitdf <- subset(cdadat, select = c(id_1,
                                   school_id_3,
                                   cda_uniqueID_2,
                                   PC1_tcs,
                                   PC2,
                                   gender,
                                   ethnicity,
                                   age_dich,
                                   Family_affluence,
                                   Parent_control_cs,
                                   Parent_care_cs,
                                   cda_com_size_2cs,
                                   cda_rat_2,
                                   typegencom3,
                                   Transit,
                                   CntN,
                                   Tau))


# for SU
fitdf1 <- subset(fitdf, select=-c(PC2))
fitdf1 <- na.omit(fitdf1)

M2f <- lmer(PC1_tcs ~ gender + age_dich + ethnicity + Family_affluence +  
              Parent_control_cs + Parent_care_cs + (1 | cda_uniqueID_2),
            data   = fitdf1,
            REML   = TRUE,
            control = lmerControl(optimizer ="Nelder_Mead"))


M3f <- lmer(PC1_tcs ~ gender + age_dich + ethnicity + Family_affluence +  
              Parent_control_cs + Parent_care_cs  + 
              cda_com_size_2cs + cda_rat_2 +  
              typegencom3 + Transit + CntN + Tau +
              (1 | cda_uniqueID_2),
            data   = fitdf1,
            REML   = TRUE,
            control = lmerControl(optimizer ="Nelder_Mead"))


mc4 <- anova(M2f, M3f, refit = F)
mc4 <- round(mc4$`Pr(>Chisq)`[2],3)

mc5 <- anova(M3, M4, refit = F)
mc5 <- round(mc5$`Pr(>Chisq)`[2],3)

# save the comparisons
mc5
colist <- c(mc1, mc2,  mc4, mc5)
cmlist <- c("M1 & M1_1","M1 & M2", 
            "M2 & M3", "M3 & M4")
cm <- as.data.frame(cbind(cmlist, colist))
colnames(cm) <- c("Models", "Anova chi p")
write.xlsx(cm, "model pair comparisons SU WT.xlsx")

# MW models
kc1 <- anova(M1a, M1_1a, refit = F)
kc1 <- round(kc1$`Pr(>Chisq)`[2],3)

kc2 <- anova(M1a, M2a, refit = F)
kc2 <- round(kc2$`Pr(>Chisq)`[2],3)


# to fit at the same dataset
fitdf <- subset(cdadat, select = c(id_1,
                                   school_id_3,
                                   cda_uniqueID_2,
                                   PC1_tcs,
                                   PC2,
                                   gender,
                                   ethnicity,
                                   age_dich,
                                   Family_affluence,
                                   Parent_control_cs,
                                   Parent_care_cs,
                                   cda_com_size_2cs,
                                   cda_rat_2,
                                   typegencom3,
                                   Transit,
                                   CntN,
                                   Tau))


# for MW
fitdf2 <- subset(fitdf, select=-c(PC1_tcs))
fitdf2 <- na.omit(fitdf2)

M2af <- lmer(PC2 ~ gender + age_dich + ethnicity + Family_affluence +  
               Parent_control_cs + Parent_care_cs + (1 | cda_uniqueID_2),
             data   = fitdf2,
             REML   = TRUE,
             control = lmerControl(optimizer ="Nelder_Mead"))


M3af <- lmer(PC2 ~ gender + age_dich + ethnicity + Family_affluence +  
               Parent_control_cs + Parent_care_cs  + 
               cda_com_size_2cs + cda_rat_2 +  
               typegencom3 + Transit + CntN + Tau +
               (1 | cda_uniqueID_2),
             data   = fitdf2,
             REML   = TRUE,
             control = lmerControl(optimizer ="Nelder_Mead"))


kc4 <- anova(M2af, M3af, refit = F)
kc4 <- round(kc4$`Pr(>Chisq)`[2],3)

kc5 <- anova(M3a, M4a, refit = F)
kc5 <- round(kc5$`Pr(>Chisq)`[2],3)

# save the comparisons

colist <- c(kc1, kc2, kc4, kc5)
cmlist <- c("M1 & M1_1","M1 & M2", 
            "M2 & M3", "M3 & M4")
cm <- as.data.frame(cbind(cmlist, colist))
colnames(cm) <- c("Models", "Anova chi p")
write.xlsx(cm, "model pair comparisons MW WT.xlsx")
#- alternative: lrtest(M2, M3) 


# combined performance -  all models except for M1_1
# SU
cp1<- compare_performance(M1,  M2,  M3, M4, rank = T) 
cp1<- as.data.frame(cp1)
write.xlsx(cp1, "Performances_alll_basic_models_SU_WT.xlsx")

cp2<- compare_performance(M1a,  M2a, M3a, M4a, rank = T) 
cp2<- as.data.frame(cp2)
write.xlsx(cp2, "Performances_alll_basic_models_MW_WT.xlsx")


# diagnostics for M3

# make a subset and change labels first
modat <- subset(cdadat, select = c(id_1,
                                   school_id_3,
                                   cda_uniqueID_2,
                                   PC1_tcs,
                                   PC2,
                                   gender,
                                   ethnicity,
                                   age_dich,
                                   Family_affluence,
                                   Parent_control_cs,
                                   Parent_care_cs,
                                   cda_com_size_2, # add this
                                   cda_com_size_2cs,
                                   cda_rat_2,
                                   typegencom3,
                                   Transit,
                                   CntN,
                                   Tau))
head(modat)

colnames(modat) <- c("Id", "School", "Community", "SU",
                     "MW", "Gender", "Ethnicity", "Age", "Family.affluence",
                     "Parent.control", "Parent.care", "Community_size", "Com.size", "ROTC", 
                     "Gen.com.", "Transit.", "Centr.", "Hierar.")
# for SU
modat1 <- subset(modat, select=-c(MW))
modat1 <- na.omit(modat1)
# with clean labels
M3cl <-lmer(SU ~ Gender + Age + Ethnicity + Family.affluence + 
              Parent.control + Parent.care +
              Com.size + ROTC + Gen.com. + 
              Transit. + Centr. + Hierar. + 
              (1|Community), 
            data   = modat1,
            REML   = TRUE,
            control = lmerControl(optimizer ="Nelder_Mead"))

# Lavene's test

plot(resid(M3cl), modat1$SU)

# Assumption 2 Homogeneity of Variance - homoskedasticity assumption

modat1$M3cl.Res<- residuals(M3cl) #extracts the residuals and places them in a new column in our original data table
modat1$Abs.M3cl.Res <-abs(modat1$M3cl.Res) #creates a new column with the absolute value of the residuals
modat1$M3cl.Res2 <- modat1$Abs.M3cl.Res^2 #squares the absolute values of the residuals to provide the more robust estimate
Levene.M3cl <- lm(M3cl.Res2 ~ Community, data=modat1) #ANOVA of the squared residuals
Levene <- anova(Levene.M3cl)
Levene <- as.data.frame(Levene)
write.xlsx(Levene, "Lavene_test_SU_M3.xlsx")
# for MW
modat2 <- subset(modat, select=-c(SU))
modat2 <- na.omit(modat2)
# with clean labels
M3acl <-lmer(MW ~ Gender + Age + Ethnicity + Family.affluence + 
               Parent.control + Parent.care +
               Com.size + ROTC + Gen.com. + 
               Transit. + Centr. + Hierar. + 
               (1|Community), 
             data   = modat2,
             REML   = TRUE,
             control = lmerControl(optimizer ="Nelder_Mead"))

plot(resid(M3acl), modat2$MW)

# Assumption 2 Homogeneity of Variance - homoskedasticity assumption

modat2$M3acl.Res<- residuals(M3acl) #extracts the residuals and places them in a new column in our original data table
modat2$Abs.M3acl.Res <-abs(modat2$M3acl.Res) #creates a new column with the absolute value of the residuals
modat2$M3acl.Res2 <- modat2$Abs.M3acl.Res^2 #squares the absolute values of the residuals to provide the more robust estimate
Levene.M3acl <- lm(M3acl.Res2 ~ Community, data=modat2) #ANOVA of the squared residuals
Levene <- anova(Levene.M3acl) 
Levene <- as.data.frame(Levene)
write.xlsx(Levene, "Lavene_test_MW_M3.xlsx")


# other diagnostics
x1 <- check_collinearity(M3cl)
plot(x1)
x2 <- check_collinearity(M3acl)
plot(x2)

check1 <- check_normality(M3cl) # OK: residuals appear as normally distributed (p = 0.309).
plot(check1, type = "qq")
check2 <- check_normality(M3cl, effects = "random") #  OK: random effects appear as normally distributed (p = 0.073).
plot(check2)

check1a <- check_normality(M3acl) 
plot(check1a, type = "qq")
check2a <- check_normality(M3acl, effects = "random") 
plot(check2a)


# plots for fixed effects for M3

plot_model(M3cl, show.values = TRUE, digits = 2, 
           value.size = 3,
           dot.size = 1)

plot_model(M3acl, show.values = TRUE, digits = 2, 
           value.size = 3,
           dot.size = 1)

# cattepillar plots

pdf(file = "Caterpillar_WT_SU_M3.pdf",   
    width = 5, # The width of the plot in inches
    height = 3.5) 
plotREsim(REsim(M3cl))
dev.off()
pdf(file = "Caterpillar_WT_MW_M3.pdf",   
    width = 5, # The width of the plot in inches
    height = 3.5) 
plotREsim(REsim(M3acl))
dev.off()

plotREsim(REsim(M1))
plotREsim(REsim(M1a))

# plot predicted values by com. property

cps <- c("Com.size", "Gen.com.", "ROTC", "Transit.","Centr.", "Hierar.")
#create plot dataframe
# SU
plot_it <-function(model, com_pro){
  plot_data <- ggpredict(model, terms = c(com_pro))
  # plot_data
  #create plot
  plot_data %>%
    plot() + 
    #add ggplot2 as needed
    theme_blank() + ylim(c(-2.2,2.9)) + ggtitle("")
  
}

top1 <- plot_it(M3cl, cps[1])
top2 <- plot_it(M3cl, cps[2])
top3 <- plot_it(M3cl, cps[3])
top4 <- plot_it(M3cl, cps[4])
top5 <- plot_it(M3cl, cps[5])
top6 <- plot_it(M3cl, cps[6])

pdf(file = paste(com_pro, "_SU_M3_WT.pdf", sep = ""),   
    width = 15, # The width of the plot in inches
    height = 10.5)
(top1 + top2 + top3)/(top4 + top5 + top6)
dev.off()

top1 <- plot_it(M3acl, cps[1])
top2 <- plot_it(M3acl, cps[2])
top3 <- plot_it(M3acl, cps[3])
top4 <- plot_it(M3acl, cps[4])
top5 <- plot_it(M3acl, cps[5])
top6 <- plot_it(M3acl, cps[6])

pdf(file = paste(com_pro, "_MW_M3_WT.pdf", sep = ""),   
    width = 15, # The width of the plot in inches
    height = 10.5)
(top1 + top2 + top3)/(top4 + top5 + top6)
dev.off()


# some advanced descriptives
length(unique(cdadat$cda_uniqueID_2))
# aggregate SU and MW
td <- subset(cdadat, select = c(cda_uniqueID_2, PC1_tcs, PC2))
# means
ag1 <- td %>%
  group_by(cda_uniqueID_2) %>%
  summarize_all(mean, na.rm = TRUE)

ag1 <- as.data.frame(ag1)
nrow(ag1)# 387

# medians
ag2 <- td %>%
  group_by(cda_uniqueID_2) %>%
  summarize_all(median, na.rm = TRUE)

ag2 <- as.data.frame(ag2)

# combine it with community level data modat (labeled) - with Gen.com as numeric
df = modat[!duplicated(modat$Community),]
colnames(df)
df1<- subset(df, select= c(Community,
                           Community_size,  
                           Gen.com., ROTC, 
                           Transit., Centr., Hierar.))


# turn Gen.com. to numeric
df1$Gen.com. <- ifelse(df1$Gen.com. == "mixed", 0,
                       ifelse(df1$Gen.com. == "female", 1, -1))

forsc<- df1
forsc <- merge(df1, ag1, by="Community")

scatterplot(forsc$Community_size, forsc$SU, boxplots = F,
            xlab = "Community Size",
            ylab = "SU",
            ellipse = F,
            grid = T)

scatterplot(forsc$Community_size, forsc$MW, boxplots = F,
            xlab = "Community Size",
            ylab = "MW",
            ellipse = F,
            grid = T)

scatterplot(forsc$Gen.com., forsc$SU, boxplots = F,
            xlab = "Gender composition",
            ylab = "SU",
            ellipse = F,
            grid = T)

scatterplot(forsc$Gen.com., forsc$MW, boxplots = F,
            xlab = "Gender composition",
            ylab = "MW",
            ellipse = F,
            grid = T)


scatterplot(forsc$ROTC, forsc$SU, boxplots = F,
            xlab = "ROTC",
            ylab = "SU",
            ellipse = F,
            grid = T)

scatterplot(forsc$ROTC, forsc$MW, boxplots = F,
            xlab = "ROTC",
            ylab = "MW",
            ellipse = F,
            grid = T)

scatterplot(forsc$Transit., forsc$SU, boxplots = F,
            xlab = "Transitivity",
            ylab = "SU",
            ellipse = F,
            grid = T)

scatterplot(forsc$Transit., forsc$MW, boxplots = F,
            xlab = "Transitivity",
            ylab = "MW",
            ellipse = F,
            grid = T)

scatterplot(forsc$Centr., forsc$SU, boxplots = F,
            xlab = "Centralization",
            ylab = "SU",
            ellipse = F,
            grid = T)

scatterplot(forsc$Centr., forsc$MW, boxplots = F,
            xlab = "Centralization",
            ylab = "MW",
            ellipse = F,
            grid = T)

scatterplot(forsc$Hierar., forsc$SU, boxplots = F,
            xlab = "Hierarchy",
            ylab = "SU",
            ellipse = F,
            grid = T)

scatterplot(forsc$Hierar., forsc$MW, boxplots = F,
            xlab = "Hierarchy",
            ylab = "MW",
            ellipse = F,
            grid = T)

#################################################
# ICC for  models M3 and M4 
################################################

iccSU<- list()
iccMW <- list()

  #pc 1
  
  M5a <-lmer(PC1_tcs ~ gender + age_dich + ethnicity + Family_affluence + 
               Parent_control_cs + Parent_care_cs +
               cda_com_size_2cs + cda_rat_2 + typegencom3 + 
               Transit + CntN + Tau + 
               (1|cda_uniqueID_2), 
             data   = cdadat,
             REML   = TRUE,
             control = lmerControl(optimizer ="Nelder_Mead"))
  
  
  M6a <- lmer(PC1_tcs ~ gender + age_dich + ethnicity + Family_affluence + 
                Parent_control_cs + Parent_care_cs +
                cda_com_size_2cs + cda_rat_2 + typegencom3 + 
                Transit + CntN + Tau +
                net_size_3cat + 
                cda_Modularity_3 + per_F_n1_3 +
                (1|cda_uniqueID_2), 
              data   = cdadat,
              REML   = TRUE,
              control = lmerControl(optimizer ="Nelder_Mead"))
  
  
  #pc 2
  
  M5b <-lmer(PC2 ~ gender + age_dich + ethnicity + Family_affluence + 
               Parent_control_cs + Parent_care_cs +
               cda_com_size_2cs + cda_rat_2 + typegencom3 + 
               Transit + CntN + Tau + 
               (1|cda_uniqueID_2), 
             data   = cdadat,
             REML   = TRUE,
             control = lmerControl(optimizer ="Nelder_Mead"))
  
  
  M6b <- lmer(PC2 ~ gender + age_dich + ethnicity + Family_affluence + 
                Parent_control_cs + Parent_care_cs +
                cda_com_size_2cs + cda_rat_2 + typegencom3 + 
                Transit + CntN + Tau +
                net_size_3cat + 
                cda_Modularity_3 + per_F_n1_3 +
                (1|cda_uniqueID_2), 
              data   = cdadat,
              REML   = TRUE,
              control = lmerControl(optimizer ="Nelder_Mead"))
  
  AdIccM3a <- round(as.numeric(icc(M5a)[1]),2)
  AdIccM4a <- round(as.numeric(icc(M6a)[1]),2)
  AdIccM3b <- round(as.numeric(icc(M5b)[1]),2)
  AdIccM4b <- round(as.numeric(icc(M6b)[1]),2)
  
  iccSU[[k]] <- c(corNames[k], AdIccM3a, AdIccM4a)
  iccMW[[k]] <- c(corNames[k], AdIccM3b, AdIccM4b)






