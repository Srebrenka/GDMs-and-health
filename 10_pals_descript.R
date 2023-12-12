###################
#     PaLS 10     # 
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

# creating the full dataset for MLM
# creating some new variables on level 2 and 3
#  descriptive analysis of communities and health outcomes 



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

library(openxlsx)
library(readxl)
library(data.table)
library(dplyr)
library(plyr)
library(sjPlot)
library(psych)
library(ggplot2)

#########################
#                       #
#     Main Script       #
#                       #
#########################
# all datasets/outputs created in 05 -09 etc scripts are moved to folder named "PalsComSets" in Downloads

# read and merge all datasets

mypath= "/PalsComSets"
multimerge = function(mypath){
  filenames=list.files(mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x) {read_excel(path = x, col_names = T)})
  Reduce(function(x,y) {merge(x,y, by="id_1")}, datalist)}# name is later id_1

tbl <- multimerge(mypath) 
dim(tbl) # 3649  100 :)

# read in and merge attribute data

rad1 <- read_excel("~/Downloads/attribute_data_pals2_FIN.xlsx", col_names = T)
nrow(rad1) #  3194
rad1 <- mutate_all(rad1, function(x) as.numeric(as.character(x)))
rad1$PALSQID <- as.factor(rad1$PALSQID)
rad1$FAS  <- as.factor(rad1$FAS)
sum(is.na(rad1$PC1))

# creating a variable for matching
rad1$code <- rad1$PALSQID
tbl$code <- tbl$id_1
mydata <- merge(rad1,tbl, by = "code") # important


nrow(tbl) - nrow(rad1) # 455 has net data, but has no attribute data
nrow(rad1) - nrow(mydata) # 46 has attribute data, but no network data
colnames(mydata)



# make new variables-----------------------------------------------
# for percentage of F gender per school #*** directly from all data laetr
table(mydata$gender)
# read in from table
geninfo <- read_excel("~/Downloads/Pals_Main/newtable_pals.xlsx")
sum(is.na(mydata$per_F_n1_3)) # 46 
sum(is.na(mydata$bia_school_id_3)) # 46
mydata$per_F_n1_3 <- mapvalues(mydata$wu_school_id_3, 
                               from = geninfo$school, 
                               to = geninfo$perF_att_data)
mydata$per_F_n1_3 <- round(mydata$per_F_n1_3, 1)
max(mydata$per_F_n1_3, na.rm = T)

mydata$per_F_n2_3 <- mapvalues(mydata$wu_school_id_3, 
                               from = geninfo$school, 
                               to = geninfo$perF_net)
mydata$per_F_n2_3 <- round(mydata$per_F_n2_3, 1)
max(mydata$per_F_n2_3)
# level 2 descriptive - per community


# make new variables
# for percentage of F gender per community (level2) - PROBLEM !!!! ****

agfg1 <- aggregate(gender ~ cfg_uniqueID_2, mydata, FUN = sum, na.rm = TRUE)
agfg2 <- aggregate(cfg_net_size_3 ~ cfg_uniqueID_2, mydata, FUN = mean, na.rm = TRUE)
agfg3 <- aggregate(cfg_com_size_2 ~ cfg_uniqueID_2, mydata, FUN = mean, na.rm = TRUE)

agfg4 <- merge(agfg1, agfg2, by="cfg_uniqueID_2")
agfg4 <- merge(agfg4, agfg3, by="cfg_uniqueID_2")
dim(agfg4)
# checking
agfg4$check1 <- agfg4$cfg_com_size_2 - agfg4$gender

min(agfg4$check1)



# new variable creation starts here:
# to level 1 by level 2 variable
agfg1 <- aggregate(gender ~ cfg_uniqueID_2, mydata, FUN = sum, na.rm = TRUE)
mydata$cfg_sum_f_2 <- mapvalues(mydata$cfg_uniqueID_2, 
                                from = agfg1$cfg_uniqueID_2, 
                                to = agfg1$gender)
mydata$cfg_per_f_2 <- round((mydata$cfg_sum_f_2/mydata$cfg_com_size_2)*100, 2)
max(mydata$cfg_per_f_2) # max is 100

# repeat for other CDAs***
aglo1 <- aggregate(gender ~ lo_uniqueID_2, mydata, FUN = sum, na.rm = TRUE)
mydata$lo_sum_f_2 <- mapvalues(mydata$lo_uniqueID_2, 
                               from = aglo1$lo_uniqueID_2, 
                               to = aglo1$gender)
mydata$lo_per_f_2 <- round((mydata$lo_sum_f_2/mydata$lo_com_size_2)*100, 2)
max(mydata$lo_per_f_2)

ageu1 <- aggregate(gender ~ eu_uniqueID_2, mydata, FUN = sum, na.rm = TRUE)
mydata$eu_sum_f_2 <- mapvalues(mydata$eu_uniqueID_2, 
                               from = ageu1$eu_uniqueID_2, 
                               to = ageu1$gender)
mydata$eu_per_f_2 <- round((mydata$eu_sum_f_2/mydata$eu_com_size_2)*100, 2)
max(mydata$eu_per_f_2)

aged1 <- aggregate(gender ~ ed_uniqueID_2, mydata, FUN = sum, na.rm = TRUE)
mydata$ed_sum_f_2 <- mapvalues(mydata$ed_uniqueID_2, 
                               from = aged1$ed_uniqueID_2, 
                               to = aged1$gender)
mydata$ed_per_f_2 <- round((mydata$ed_sum_f_2/mydata$ed_com_size_2)*100, 2)
max(mydata$ed_per_f_2)

agwu1 <- aggregate(gender ~ wu_uniqueID_2, mydata, FUN = sum, na.rm = TRUE)
mydata$wu_sum_f_2 <- mapvalues(mydata$wu_uniqueID_2, 
                               from = agwu1$wu_uniqueID_2, 
                               to = agwu1$gender)
mydata$wu_per_f_2 <- round((mydata$wu_sum_f_2/mydata$wu_com_size_2)*100, 2)
max(mydata$wu_per_f_2)

agwd1 <- aggregate(gender ~ wd_uniqueID_2, mydata, FUN = sum, na.rm = TRUE)
mydata$wd_sum_f_2 <- mapvalues(mydata$wd_uniqueID_2, 
                               from = agwd1$wd_uniqueID_2, 
                               to = agwd1$gender)
mydata$wd_per_f_2 <- round((mydata$wd_sum_f_2/mydata$wd_com_size_2)*100, 2)
max(mydata$wd_per_f_2)

agimu1 <- aggregate(gender ~ imu_uniqueID_2, mydata, FUN = sum, na.rm = TRUE)
mydata$imu_sum_f_2 <- mapvalues(mydata$imu_uniqueID_2, 
                                from = agimu1$imu_uniqueID_2, 
                                to = agimu1$gender)
mydata$imu_per_f_2 <- round((mydata$imu_sum_f_2/mydata$imu_com_size_2)*100, 2)
max(mydata$imu_per_f_2)

agimd1 <- aggregate(gender ~ imd_uniqueID_2, mydata, FUN = sum, na.rm = TRUE)
mydata$imd_sum_f_2 <- mapvalues(mydata$imd_uniqueID_2, 
                                from = agimd1$imd_uniqueID_2, 
                                to = agimd1$gender)
mydata$imd_per_f_2 <- round((mydata$imd_sum_f_2/mydata$imd_com_size_2)*100, 2)
max(mydata$imd_per_f_2)

agbia1 <- aggregate(gender ~ bia_uniqueID_2, mydata, FUN = sum, na.rm = TRUE)
mydata$bia_sum_f_2 <- mapvalues(mydata$bia_uniqueID_2, 
                                from = agbia1$bia_uniqueID_2, 
                                to = agbia1$gender)
mydata$bia_per_f_2 <- round((mydata$bia_sum_f_2/mydata$bia_com_size_2)*100, 2)
max(mydata$bia_per_f_2)

agsbm1 <- aggregate(gender ~ sbm_uniqueID_2, mydata, FUN = sum, na.rm = TRUE)
mydata$sbm_sum_f_2 <- mapvalues(mydata$sbm_uniqueID_2, 
                                from = agsbm1$sbm_uniqueID_2, 
                                to = agsbm1$gender)
mydata$sbm_per_f_2 <- round((mydata$sbm_sum_f_2/mydata$sbm_com_size_2)*100, 2)
max(mydata$sbm_per_f_2)


#-----------------------------------------------------------------------------
# create Modularity for "bia" - later add it in the code
#----------------------------------------------------------------------------
mods=c() # schools$school
for (i in 1:22){
  #Pull out individual schools#
  school1 = mydata %>% 
    filter(bia_school_id_3 == i)
  school1$name <- as.character(school1$id_1)
  tomatch <- school1[school1$name %in% V(palnet[[i]])$name,]
  tomatch <- tomatch[match(V(palnet[[i]])$name, tomatch$name), ]
  
  V(palnet[[i]])$memb <- tomatch$bia_uniqueID_2
  mo <- modularity(palnet[[i]], V(palnet[[i]])$memb)
  mods[i] <- round(mo, 3)
}

netor <- seq(1, 22)

mydata$bia_Modularity_3 <- mapvalues(mydata$bia_school_id_3, 
                                     from = netor, 
                                     to = mods)
# it is there now
dim(mydata) # 3148  147
# SAVING THE DATA:
write.table(mydata, file = "all_ml_datDEC.csv", dec = ',' ,
            sep = ";",row.names = FALSE, col.names = TRUE)
write.xlsx(mydata, "all_ml_datDEC.xlsx")

# level 2 descriptive----------------------------------------------

# Ncom average, min and max per school
# transpose and make it average per CDA?
# add mod and Ncom to school descriptives

# N of com, Modularity and average com size per school (rows)
# and CDAs (columns)
ag1 <- aggregate(cbind(wd_net_size_3, # any
                       cfg_N_com_3,cfg_Modularity_3, cfg_com_size_2,
                       eu_N_com_3, eu_Modularity_3, eu_com_size_2,
                       ed_N_com_3, ed_Modularity_3, ed_com_size_2,
                       imu_N_com_3, imu_Modularity_3, imu_com_size_2,
                       imd_N_com_3, imd_Modularity_3, imd_com_size_2,
                       lo_N_com_3, lo_Modularity_3, lo_com_size_2,
                       wu_N_com_3, wu_Modularity_3, wu_com_size_2,
                       wd_N_com_3, wd_Modularity_3, wd_com_size_2, 
                       bia_N_com_3, bia_Modularity_3, bia_com_size_2,
                       sbm_N_com_3, sbm_Modularity_3, sbm_com_size_2) 
                 ~wd_school_id_3, # any com_school, mean is equal to true value
                 # for _N_com_3 & _Modularity_3, but new per _com_size_2
                 mydata, FUN = mean, na.rm = TRUE)
colnames(ag1)[3:32] <- paste("M", colnames(ag1)[3:32], sep = "_")

# minimum size of community per school and CDAs
ag2 <- aggregate(cbind(cfg_com_size_2,
                       eu_com_size_2,
                       ed_com_size_2,
                       imu_com_size_2,
                       imd_com_size_2,
                       lo_com_size_2,
                       wu_com_size_2,
                       wd_com_size_2,
                       bia_com_size_2,
                       sbm_com_size_2) ~ wd_school_id_3, mydata, FUN = min, na.rm = TRUE)
colnames(ag2)[2:11] <- paste("min", colnames(ag2)[2:11], sep = "_")
# maximum size of community per school and CDAs
ag3 <- aggregate(cbind(cfg_com_size_2,
                       eu_com_size_2,
                       ed_com_size_2,
                       imu_com_size_2,
                       imd_com_size_2,
                       lo_com_size_2,
                       wu_com_size_2,
                       wd_com_size_2,
                       bia_com_size_2,
                       sbm_com_size_2) ~ wd_school_id_3, mydata, FUN = max, na.rm = TRUE)
colnames(ag3)[2:11] <- paste("max", colnames(ag3)[2:11], sep = "_")


# merge all

agmer <- merge(ag1, ag2, by= "wd_school_id_3")
agmer <- merge(agmer, ag3, "wd_school_id_3")

write.xlsx(agmer, "level3_2_table.xlsx")




# level1 descriptive-------------------------------------------------
# correlation between level-1 covariates

# two versions, with additional control variables and without

level1a <- subset(mydata, select = c(PC1,
                                     PC2,
                                     rcs1,
                                     rcs2, 
                                     gender,
                                     ethnicity,
                                     age,
                                     FASRAW,
                                     PBICONT,
                                     PBICARE))


level1b <- subset(mydata, select = c(PC1,
                                     PC2,
                                     rcs1,
                                     rcs2, 
                                     gender,
                                     ethnicity,
                                     age,
                                     FASRAW,
                                     PBICONT,
                                     PBICARE))
library(corrplot)

M1 = cor(level1a, use="pairwise.complete.obs", method = "pearson")
corrplot.mixed( M1,              
                lower = "number", 
                upper = "circle",
                tl.col = "black",
                number.cex = 0.8,
                tl.cex = 0.7)

M1 = cor(level1b, use="pairwise.complete.obs", method = "pearson")
corrplot.mixed( M1,              
                lower = "number", 
                upper = "circle",
                tl.col = "black",
                number.cex = 0.8,
                tl.cex = 0.6,
                tl.pos = "lt")


# repeat this:
des2 <- describe(x =level1b, skew = T, ranges = T, check = T)
des2$variables <- rownames(des2)
write.xlsx(des2, "des_level_1_netsample.xlsx", row_names = TRUE)


#---------------------------------------------------------------------------
#       N of com sizes for all schools per GDM
#----------------------------------------------------------------------------
# aggregate per community ID for each CDA - make it a loop***
mydata <- read_excel("all_ml_datDEC.xlsx", col_names = TRUE)

ag_1 <- aggregate(wd_com_size_2~ wd_uniqueID_2, data=mydata, FUN = mean,
                  na.rm = T)
Mcs <- mean(ag_1$wd_com_size_2) # 9.42
comid <- unique(mydata$wd_uniqueID_2)
comN <- length(comid) # 387

p1 <- ggplot(data=ag_1, aes(wd_com_size_2)) + 
  geom_histogram(breaks=seq(1, 106, by=1), 
                 col="black", 
                 fill="grey", 
                 alpha = .2) + 
  labs(title="Sizes of 387 communities in 22 schools - Walktrap (directed)", 
       x="Size of community (M=9.42)", y="Frequency") +
  geom_vline(aes(xintercept = Mcs), # over all communities
             colour = "darkred", linetype ="longdash", size = .8)
# add median too?

tn22 <- length(comid)
for(k in 1:length(comid)){
  csid <- aggregate(my)
}

# descriptives level 3 per CDA

csdf <- read_excel("level3_2_table.xlsx", col_names = TRUE)
