###################
#     PaLS 13     #
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

# making imputations
# MAKING A MLM-READY DATASET - version 2.0, 26/01/2022
# needs to be done for each GDM (this is the example for Walktrap /wd/)

#########################
#                       #
#    Load scripts       #
#                       #
#########################


#########################
#                       #
#    Load packages      #
#                       #
#########################

library(ggplot2)  # Graphs
library(optimx)  # required for some optimizer
library(openxlsx)
library(readxl)
library(plyr)



#########################
#                       #
#     Main Script       #
#                       #
#########################

########################
#  imputation        
########################

# preparing variables for imputation
# SET GENDER AND ETHNICITY AS FACTORS

dat3 <- read.xlsx("attribute_data_pals2_FIN.xlsx")

dat3$gender <- factor(dat3$gender,
                      levels=c(0, 1),
                      labels=c("Male",
                               "Female"))

dat3$ethnicity <- factor(dat3$ethnicity,
                         levels=c(0, 1),
                         labels=c("non-white",
                                  "white"))
dat3$FAS<- as.factor(dat3$FAS)
dat3$school_n<- as.factor(dat3$school_n)

dat3$age_dich <- ifelse(is.na(dat3$age), NA,
                        ifelse(dat3$age == 15, 0, 1))
dat3$age <- NULL

dat3$PC1_t2 <- round(dat3$PC1_t2, 2)
xxx <- base::scale(dat3$PC1_t2, center = TRUE, scale = TRUE)
class(xxx)
xxx <- as.vector(xxx)
dat3$PC1_tcs <- round(xxx, 2)
dat3$PC1_t2<- NULL

nrow(dat2) # 3194
dat4<- merge(dat2, dat3, by="PALSQID", all.y=T) # bc we need just those with net data anyway
colnames(dat4)
dat4$school_n <- dat4$school_n.x
dat4$school_n.x <- NULL
dat4$school_n.y <- NULL
nrow(dat4)# 3148
# exclude only PCs
dat5<- subset(dat4, select=-c(PC1))
dat5$school_n <- as.factor(dat5$school_n)
dat5$PALSQID <- as.factor(dat5$PALSQID)
str(dat5)
table(dat5$used.drugs2)
#palsqid to factor?
for (i in 1:length(colnames(dat5))){
  print(colnames(dat5)[i])
  print(sum(is.na(dat5[[i]])))
}



# descriptive
# make a dataset with nice names for figures
dat5x <- dat5
# reorder
dat5x<- subset(dat5x, select =c(school_n,
                                age_dich,
                                gender,
                                ethnicity, 
                                FAS,
                                PBICONT,
                                PBICARE,
                                smoking,
                                ALCFREQ,
                                used.drugs2,
                                drug.effects,
                                GHQLIKS,
                                SELFEST,
                                wr.raw.t,
                                PC1_tcs, 
                                PC2))
colnames(dat5x) <- c("School", "Age", "gender", "ethnicity", "Familly_affluence",
                     "Parent_control", "Parent_care", "Smoking", "Drinking",
                     "Using_drugs", "Drug_effects", "GHQ", "Self-esteem", "Worries",
                     "PC1", "PC2")
library(visdat)
vis_dat(dat5x)


#------------------------------------------------------------------------------
#                    IMPUTATION PROCESS
#-----------------------------------------------------------------------------
imp <- mice(dat5, maxit=0, print=F)
imp$method # looks ok
pred <- quickpred(
  dat5,
  mincor = 0.01, 
  minpuc = 0.05,
  method = "pearson")


# just one dataset:
imp2 <- mice(dat5, m = 1, maxit = 40, 
             #predictorMatrix = predM, 
             #method = meth, 
             print =  FALSE, seed=500)
imp2$meth
imp2$loggedEvents
summary(imp2)

#check imputed values
imp2$imp$PC2

completeData <- complete(imp2)
dim(completeData)
colnames(completeData)
# save it for future use
write.xlsx(completeData, "One_imputed_dataset_maxit40.xlsx")

####################################################
# joining with other (network) data 
####################################################

# read net data - long format
mypath= "/PalsComSets"
multimerge = function(mypath){
  filenames=list.files(mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x) {read_excel(path = x, col_names = T)})
  Reduce(function(x,y) {merge(x,y, by="id_1")}, datalist)}# name is later id_1

tbl <- multimerge(mypath) 
dim(tbl)
head(tbl)



# read community level data1 and map it to long format in tbl

comdat1 <- as.data.frame(read_excel("Com_level_data_allCDA_1.xlsx"))

comdat1 = comdat1[comdat1$cda == "wd", ] # walktrap
colnames(comdat1)
# read community level data2 and map it to long format in tbl 

comdat2 <- as.data.frame(read_excel("Com_level_data_allCDA_2_tauetc.xlsx"))
comdat2 = comdat2[comdat2$cda == "wd", ]
comdat2$uniqueID_2 <-comdat2$uniqueID
comdat2$uniqueID <- NULL
comdat2$cda <- NULL
# merge the two

comdatf <- merge(comdat1, comdat2, by = "uniqueID_2")


# mapping
comid <- comdatf$uniqueID_2
propF <- comdatf$propF
Tau <- comdatf$Tau
CntN <- comdatf$CntN
Transit <- comdatf$Transit

tbl$propF <- mapvalues(tbl$wd_uniqueID_2,
                       from = comid, to = propF)
tbl$Tau <- mapvalues(tbl$wd_uniqueID_2,
                     from = comid, to = Tau)
tbl$CntN <- mapvalues(tbl$wd_uniqueID_2,
                      from = comid, to = CntN)
tbl$Transit <- mapvalues(tbl$wd_uniqueID_2,
                         from = comid, to = Transit)



# read attribute data (imputed)
impd <- as.data.frame(read_excel("~/Downloads/One_imputed_dataset_maxit40.xlsx", col_names = T))

impd$id_1 <- impd$PALSQID


# keep just necessary variables
impd<- subset(impd, select = c(id_1,
                               FAS,
                               PBICONT,
                               PBICARE,
                               gender, 
                               ethnicity,
                               age_dich,
                               PC1_tcs,
                               PC2))

# new dataset - connect att(imputed) data with netdata
nidt <- merge(tbl, impd, by="id_1", all.x = T)
nrow(nidt) # IMPUTED DATASET

# reading in data about proportion of girls in school
sch <- read_xlsx("schools_results.xlsx")
schid<-sch$school
propFsc <- as.numeric(sch$Girls_prop)
nidt$per_F_n1_3 <- mapvalues(nidt$wd_school_id_3, from = schid,
                             to = propFsc)


# creating modularity score for bia
source("02_pals_create_nets.R") 

mods=c() # schools$school
for (i in 1:22){
  #Pull out individual schools#
  school1 = nidt %>% 
    filter(bia_school_id_3 == i)
  school1$name <- as.character(school1$id_1)
  tomatch <- school1[school1$name %in% V(palnet[[i]])$name,]
  tomatch <- tomatch[match(V(palnet[[i]])$name, tomatch$name), ]
  
  V(palnet[[i]])$memb <- tomatch$bia_uniqueID_2
  mo <- modularity(palnet[[i]], V(palnet[[i]])$memb)
  mods[i] <- round(mo, 2)
}

netor <- seq(1, 22)

nidt$bia_Modularity_3 <- mapvalues(nidt$bia_school_id_3, 
                                   from = netor, 
                                   to = mods)



# Create CDA specific dataset
wddat <- base::subset(nidt, select = c(id_1,             # id level 1
                                       wd_uniqueID_2,    # id level 2
                                       wd_school_id_3,   # id level 3
                                       PC1_tcs,
                                       PC2,   
                                       PC1_tcs,
                                       # l 1
                                       gender,
                                       age_dich,
                                       ethnicity,
                                       FAS,
                                       PBICONT,
                                       PBICARE,  
                                       # l 2
                                       wd_com_size_2, 
                                       wd_rat_2,
                                       propF,
                                       Transit,
                                       Tau,
                                       CntN,
                                       wd_net_size_3,
                                       wd_N_com_3,     
                                       wd_Modularity_3,
                                       per_F_n1_3 ))

# Data preparations & pre-model diagnostics

# change names

names(wddat)[names(wddat) == 'PBICONT'] <- 'Parent_control'
names(wddat)[names(wddat) == 'PBICARE'] <- 'Parent_care'
names(wddat)[names(wddat) == 'FAS'] <- 'Family_affluence'
#names(wddat)[names(wddat) == 'FASRAW'] <- 'Family_affluence'

# cda specific

names(wddat)[names(wddat) == 'wd_uniqueID_2'] <- 'cda_uniqueID_2'
names(wddat)[names(wddat) == 'wd_school_id_3'] <- 'school_id_3'
names(wddat)[names(wddat) == 'wd_com_size_2'] <- 'cda_com_size_2'
names(wddat)[names(wddat) == 'wd_rat_2'] <- 'cda_rat_2'
names(wddat)[names(wddat) == 'wd_com_perN_2'] <- 'cda_com_perN_2'
names(wddat)[names(wddat) == 'wd_per_f_2'] <- 'cda_per_f_2'
# l 3
names(wddat)[names(wddat) == 'wd_net_size_3'] <- 'net_size_3'
names(wddat)[names(wddat) == 'wd_N_com_3'] <- 'cda_N_com_3'
names(wddat)[names(wddat) == 'wd_Modularity_3'] <- 'cda_Modularity_3'    


# Centering  &  (re)scaling 


# centering:scale(x, center = TRUE, scale = FALSE)
# centering and scaling: scale(x, center = TRUE, scale = TRUE)

# factor variables as numeric
wddat$Family_affN<- ifelse(wddat$Family_affluence == "1", 1, 
                           ifelse(wddat$Family_affluence == "2",2, 
                                  ifelse(wddat$Family_affluence == "3", 3, NA)))
wddat$genderN <- ifelse(wddat$gender == "Female", 1, 0)

# alternative variables
wddat$cda_com_size_2cat <- ifelse(wddat$cda_com_size_2 <= 5, 1,
                                  ifelse(wddat$cda_com_size_2 > 5 & wddat$cda_com_size_2 <= 10, 2,
                                         ifelse(wddat$cda_com_size_2 > 10 & wddat$cda_com_size_2 <= 20, 3,
                                                ifelse(wddat$cda_com_size_2 >20, 4, NA))))

wddat$net_size_3cat <- ifelse(wddat$net_size_3 <=140, 1,
                              ifelse(wddat$net_size_3 >140 & 
                                       wddat$net_size_3 <= 220, 2,
                                     ifelse(wddat$net_size_3 >220, 3, NA)))


# sanity check
fdf <- wddat[wddat$gender == "Female" , ]
min(fdf$propF, na.rm = T)
mdf<- wddat[wddat$gender == "Male" , ]
max(mdf$propF, na.rm = T)
# new variable
wddat$comgender2 <- ifelse(wddat$propF == 1 | wddat$propF == 0, "one_gender_com", 
                           "mixed_gen_com")

# as.numeric
wddat$comgender2N <- ifelse(wddat$propF == 1 | wddat$propF == 0, 0, 1)

wddat$typegencom3 <- ifelse(wddat$propF == 1, "female",
                            ifelse(wddat$propF == 0, "male", "mixed"))

# make less precise 
wddat$cda_Modularity_3 <-round(wddat$cda_Modularity_3, 2)

# on different scale
wddat$cda_com_size_2r <- wddat$cda_com_size_2/100
wddat$cda_com_size_2r <-round(wddat$cda_com_size_2r, 1)


# centering and scaling
# L1
wddat$Parent_care_cs <- scale(wddat$Parent_care, center = TRUE, scale = TRUE)
wddat$Parent_control_cs <- scale(wddat$Parent_control, center = TRUE, scale = TRUE)
wddat$Parent_care_cs <- as.vector(wddat$Parent_care_cs)
wddat$Parent_control_cs <- as.vector(wddat$Parent_control_cs)
wddat$Parent_care_cs <- round(wddat$Parent_care_cs, 1)
wddat$Parent_control_cs <- round(wddat$Parent_control_cs, 1)
# L2

wddat$cda_com_size_2cs <- as.vector(scale(wddat$cda_com_size_2, center = T, scale = T))
wddat$cda_com_size_2cs <- round(wddat$cda_com_size_2cs, 3)


wddat2 <- subset(wddat, select= c(PC1_tcs,
                                  PC2, 
                                  genderN,
                                  age_dich,
                                  Family_affN, 
                                  cda_com_size_2,
                                  cda_rat_2,
                                  Transit,
                                  CntN, 
                                  Tau,
                                  net_size_3,
                                  net_size_3cat,
                                  cda_Modularity_3,
                                  per_F_n1_3,
                                  propF))

colnames(wddat2)[1] <- "PC1"
colnames(wddat2)[6] <- "Com_size"
colnames(wddat2)[7] <- "ROCT"
colnames(wddat2)[13] <- "Modularity"
colnames(wddat2)[14] <- "perFschool"
colnames(wddat2)[15] <- "propFcom"


pdf(file = "Corrplot_all_vars_wd.pdf")   # The directory you want to save the file in

M = cor(wddat2, use="pairwise.complete.obs", method = "pearson")
corrplot.mixed( M,              
                lower = "number", 
                upper = "circle",
                tl.col = "black",
                number.cex = 0.6,
                tl.cex = 0.6,
                tl.pos = "lt")
dev.off()
# save the dataset
write.xlsx(wddat, "GDMdata.xlsx") # the "GDM" should be replaced with GDM acronym




# optional
######################################################################
#    some additional descriptives and info
######################################################################

# info on gender composition of communities

t1 <- as.data.frame(table(wddat$typegencom3)) # make it to table!!!
colnames(t1) <- c("type_gencom_wd", "Freq")
t1$Prct <- round((t1$Freq/nrow(wddat))*100, 2)
wddat$comgenN <- ifelse(wddat$propF == 1 | wddat$propF == 0, 0, 1)

wddat1<- wddat[!wddat$cda_com_size_2 ==1, ] 
wddat2<- wddat[!wddat$cda_com_size_2 <=2, ]
wddat3<- wddat[!wddat$cda_com_size_2 <=3, ]
wddat30p<- wddat[!wddat$cda_com_size_2 >30, ]



comdatf$typegencom3 <- ifelse(comdatf$propF == 1, "female",
                              ifelse(comdatf$propF == 0, "male", "mixed"))

t3 <- as.data.frame(table(comdatf$typegencom3))
colnames(t3) <- c("type_gencom_wd", "Freq_com")
t3$Prct_com <- round((t3$Freq_com/nrow(comdatf))*100, 2)


# exclude com size 1
comdatf1 <- comdatf[!comdatf$c_size == 1, ]
t4 <- as.data.frame(table(comdatf1$typegencom3)) # make it to table!!!
colnames(t4) <- c("type_gencom_wd", "Freq_com1")
t4$Prct_com1 <- round((t4$Freq_com1/nrow(comdatf1))*100, 2)
t3$Freq_com1 <- t4$Freq_com1
t3$Prct_com1 <- t4$Prct_com1

t3
t5<- merge(t1, t3, by = "type_gencom_wd")
write.xlsx(t1, "wd_outputs/type_gendercom_wd.xlsx")

t5

# make pie figures

p1 <- ggplot(t5, aes(x="", y=Prct, fill=type_gencom_wd)) +
  geom_bar(stat="identity", width=1, color="white") +
  geom_text(aes(label = Prct),
            position = position_stack(vjust = 0.5), color = "white", size=3) +
  coord_polar("y", start=0) +
  theme_void() + 
  scale_fill_brewer(palette="Set1")

p2 <- ggplot(t5, aes(x="", y=Prct_n1plus, fill=type_gencom_wd)) +
  geom_bar(stat="identity", width=1, color="white") +
  geom_text(aes(label = Prct_n1plus),
            position = position_stack(vjust = 0.5), color = "white", size=3) +
  coord_polar("y", start=0) +
  theme_void() + 
  scale_fill_brewer(palette="Set1")

p3<- ggplot(t5, aes(x="", y=Prct_com, fill=type_gencom_wd)) +
  geom_bar(stat="identity", width=1, color="white") +
  geom_text(aes(label = Prct_com),
            position = position_stack(vjust = 0.5), color = "white", size=3) +
  coord_polar("y", start=0) +
  theme_void() + 
  scale_fill_brewer(palette="Set1")

p4<- ggplot(t5, aes(x="", y=Prct_com1, fill=type_gencom_wd)) +
  geom_bar(stat="identity", width=1, color="white") +
  geom_text(aes(label = Prct_com1),
            position = position_stack(vjust = 0.5), color = "white", size=3) +
  coord_polar("y", start=0) +
  theme_void() + 
  scale_fill_brewer(palette="Set1")

pdf(file = "wd_outputs/gender_communities.pdf") 
p1 + p2 + p3 + p4 + plot_layout(ncol=2)
dev.off()


# other
# checking skewness
ivs_vec2 <- names(wddat)[c(12:28)] # 

# level 1
for (i in ivs_vec2) {
  print(i)
  print(shapiro.test(wddat[[i]]))
} # all sig.
for (i in ivs_vec2) {
  print(i)
  print(e1071::skewness(wddat[[i]], na.rm = TRUE))
} # moderate to strong for FASRAW, PBICONT, PBICARE






# some figures - not essential

pdf(file = "wd_outputs/Histogram of Modularity at level 1_wd.pdf",
    height = 4, width = 5) 
ggplot(data=wddat) +
  geom_histogram(aes(x=cda_Modularity_3, color=I("black"),fill=I("goldenrod")))
dev.off()

# Rat 2
pdf(file = "wd_outputs/Histogram of rat at level 1_wd.pdf",
    height = 4, width = 5) 
ggplot(data=wddat) +
  geom_histogram(aes(x=cda_rat_2, color=I("black"),fill=I("goldenrod")))
dev.off()


ggplot(data=wddat) +
  geom_histogram(aes(x=cda_rat_2, color=I("black"),fill=I("goldenrod")))
#### THIS IS CFA SPECIFIC TRANSORMATION!!! let's make it for all CDAs
##  absolute - better to change it to community size

# community size 2
pdf(file = "wd_outputs/Histogram of community size at level 1_wd.pdf",
    height = 4, width = 5) 
ggplot(data=wddat) +
  geom_histogram(aes(x=cda_com_size_2, color=I("black"),fill=I("goldenrod")))
dev.off()

#### THIS IS CFA SPECIFIC TRANSORMATION!!! let's make it for all CDAs
# similar pattern should be for other CDAs due to network level segregation
# absolute
pdf(file = "wd_outputs/Histogram of percentage of females in the community at level 1_wd.pdf",
    height = 4, width = 5) 
ggplot(data=wddat) +
  geom_histogram(aes(x=propF, color=I("black"),fill=I("goldenrod")))
dev.off()
#  THERE IS SEGRAGATION, GENDER HOMOPHILY - BUT LESS THAN FOR walktrap directed

# not sure about this:
# Transform cda_per_F_2 variable to another - same-gender percentage

# Rescaling - part 2: centering and scaling

# One person and dyads are not communities, and triads?


wddat1<- wddat[!wddat$cda_com_size_2 ==1, ]
wddat2<- wddat[!wddat$cda_com_size_2 <=2, ]
wddat3<- wddat[!wddat$cda_com_size_2 <=3, ]


# -----------------------------------------------------------------------------
#              Checking the relationship between l2, l3 predictors and DVs 
#                                   LEVEL 1
#------------------------------------------------------------------------------

# check the relationship between community size and perF!

# use Gally for this - on

# bigger community - less same gender percentage
cor.test(wddat1$propF,wddat1$cda_com_size_2, method="pearson") # -0.16


# VERY IMPORTANT:
wddat1 %>%  
  ggplot(aes(x=propF,y=cda_com_size_2, color= gender)) + geom_point(alpha=0.2) +
  geom_smooth(method = "loess") 

wddat1 %>%  
  ggplot(aes(x=propF,y=cda_com_size_2, color= gender)) + geom_point(alpha=0.2) +
  geom_smooth(method = "lm") # linear relationship shows different story !!!


# level 3
# Modularity
wddat %>%  
  ggplot(aes(x=cda_Modularity_3,y=PC1_tcs)) + geom_point(alpha=0.2) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "darkred") # not relation

wddat %>%  
  ggplot(aes(x=cda_Modularity_3,y=PC1_tcs, color=gender)) + geom_point(alpha=0.2) +
  geom_smooth(method = "loess") 

# N_com_3
wddat %>%  
  ggplot(aes(x=cda_N_com_3,y=PC1_tcs)) + geom_point(alpha=0.2) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "darkred") # no relation

wddat %>%  
  ggplot(aes(x=cda_N_com_3,y=PC1_tcs, color = gender)) + geom_point(alpha=0.2) +
  geom_smooth(method = "loess")

# per_F_n1_3
wddat %>%  
  ggplot(aes(x=per_F_n1_3,y=PC1_tcs)) + geom_point(alpha=0.2) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "darkred") # 

wddat %>%  
  ggplot(aes(x=per_F_n1_3,y=PC1_tcs, color = gender)) + geom_point(alpha=0.2) +
  geom_smooth(method = "loess")

# level 2

# cda_per_f_2 = THIS IS VERY IMPORTANT FIGURE THAT HIDES SOMETHING
wddat %>%  
  ggplot(aes(x=propF,y=PC1_tcs, color=gender)) + geom_point(alpha=0.2) +
  geom_smooth(method = "loess")# method=lm
wddat %>%  
  ggplot(aes(x=propF,y=PC1_tcs)) + geom_point(alpha=0.2) +
  geom_smooth(method = "loess")# overall - more women in community worse for PC1

wddat %>%  
  ggplot(aes(x=propF,y=PC2)) + geom_point(alpha=0.2) +
  geom_smooth(method = "loess")# overall - more women in community worse for PC2

wddat %>%  
  ggplot(aes(x=propF,y=PC2, color=gender)) + geom_point(alpha=0.2) +
  geom_smooth(method = "loess") # different for males and females !!!!!!!


wddat %>%  
  ggplot(aes(x=cda_com_size_2,y=propF)) + geom_point(alpha=0.2) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "darkred") # superinteresting!!!!
# bigger communities - less women, but not-linear!


# com size 2
wddat %>%  
  ggplot(aes(x=cda_com_size_2,y=PC1_tcs)) + geom_point(alpha=0.2) +
  geom_smooth(method = "loess") # non-linear
wddat %>%  
  ggplot(aes(x=cda_com_size_2,y=PC1_tcs, color = gender)) + geom_point(alpha=0.2) +
  geom_smooth(method = "loess") # non-linear and gender differences !!!!
wddat %>%  
  ggplot(aes(x=cda_com_size_2,y=PC2)) + geom_point(alpha=0.2) +
  geom_smooth(method = "loess") # non-linear

# rat 2 
wddat %>%  
  ggplot(aes(x=cda_rat_2,y=PC1_tcs)) + geom_point(alpha=0.2) +
  geom_smooth(method = "lm") 
wddat %>%  
  ggplot(aes(x=cda_rat_2,y=PC1_tcs, color= gender)) + geom_point(alpha=0.2) +
  geom_smooth(method = "lm") # gender

wddat %>%  
  ggplot(aes(x=cda_rat_2,y=PC2)) + geom_point(alpha=0.2) +
  geom_smooth(method = "lm") 
wddat %>%  
  ggplot(aes(x=cda_rat_2,y=PC2, color= gender)) + geom_point(alpha=0.2) +
  geom_smooth(method = "lm") # gender

colnames(wddat)

cor.test(geninfo$perF_att_data, geninfo$perF_net)


