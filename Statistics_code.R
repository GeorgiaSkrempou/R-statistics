####################################################################################
################################ VERY IMPORTANT ####################################
##### The names of the csv files must contain the word "firing" for the firing #####
##### frequency file and "CCIO" for the CCIO summary files. #######################
##### The code will not run correctly in a different case.  #######################
###################################################################################

# load packages #
library(psych) #describe
library(pastecs) #normality test
library(outliers) #Grubbs test
library(ggplot2) #graphs
library(extrafont) #fonts for graphs
library(grid) #don't know if I need this
library(gridExtra) #arrange graphs
library(tibble) #don't need this I think
library(dplyr)
library(naniar)
library(reshape2) #for melt function
library(stringr) #for detection of words in character strings
library(matrixStats)
library(tcltk2)
library(tcltk)
library(car) #levene test & ANOVA
library(data.table) #for transpose
library(afex) #repeated measures anova
library(rstatix) #shapiro test function
library(ez) #anova
library(lme4) #linear mixed models

library(DescTools)

library(Hmisc)
library(multcomp)
library(haven)


###########################################################################
########################### SECTION 1 #####################################
###################### File choosing ######################################
###########################################################################

# choose files and set working directory #
cat("\nChoose a pair of firing frequency and CCIO files for analysis. Do not chose more or less than 1 pair of files. ONLY csv FILES.\n")
fileLinks <- choose.files() # opens a window to choose the files
parent <- dirname(fileLinks) #sets working directory as the folder where you chose your files from
setwd(parent) #sets working directory as the folder where you chose your files from
getwd()

# performing various checks on the names of the data files to make sure we are reading the correct data #
if (length(fileLinks) == 2) { #check whether only 2 files are chosen
  Link1 <- fileLinks[1]
  Link2 <- fileLinks[2]
  #check if the files are firing frequency and CCIO files
  file.check.ccio.1 <- str_detect(Link1, "(CCIO)|(ccio)") #check if the file name contains the words CCIO or ccio
  file.check.firing.1 <- str_detect(Link1, "(Firing)|(firing)") #check if the file name contains the words Firing or firing
  file.check.ccio.2 <- str_detect(Link2,"(CCIO)|(ccio)")
  file.check.firing.2 <- str_detect(Link2, "(Firing)|(firing)")
  #check if the files are csv files
  file.check.csv.1 <- str_detect(Link1, "csv") 
  file.check.csv.2 <- str_detect(Link2, "csv")
  if (file.check.ccio.1 == TRUE & file.check.firing.2 == TRUE & file.check.csv.1 == TRUE & file.check.csv.2 == TRUE) {
    Link1<-fileLinks[2]
    Link2<-fileLinks[1]
    cat("You can proceed with reading the data")
  } else if (file.check.firing.1 == TRUE & file.check.ccio.2 ==TRUE & file.check.csv.1 == TRUE & file.check.csv.2 == TRUE) {
    cat("You can proceed with reading the data")
  } else {
    stop("There is something wrong with the data files. Check if the data you enteredare csv files and if the filenames are correct.")
  }
} else if (length(fileLinks) != 2) {
  stop("You must select exactly 2 files. Re-run section 1.")
}



##########################################################################################
######################### SECTION 2 ######################################################
############## Read FF data and select the good cells ####################################
##########################################################################################

# read data #
ff<-read.csv(Link1,header=T, check.names=FALSE) #ff for firing frequency
colnames(ff)[1] <- "sweep"

# Check data integrity # you can also do that by opening the file from the environment section -->
headTail(ff, 5,5)
dim(ff)
names(ff)

# clean data #
ff[ff == ""] <- NA #replace empty cells with NAs
all_na <- function(x) any(!is.na(x))
ff<-ff %>% select_if(all_na) #clean columns where all rows are NA
ff50<-ff[1:50,] #select the first 50 sweeps
ff2<-ff50 %>% select_if(~any(. >= 5)) #select only those cells that have >=5APs at least in one sweep

colFF<-colnames(ff)
colFF2<-colnames(ff2)
badCells.AP.names<-setdiff(colFF,colFF2) #these are the cells with <5AP (bad cells based on AP)

rowFF2<-nrow(ff2)
test<-sapply(ff2, function(x) sum(is.na(x))) #how many NAs in each column - to see which cells have too few sweeps and need to be excluded

#999 is good cell, 666 is bad cell
for (i in 1:ncol(ff2)) {
    check.NA<-is.na(ff2[rowFF2,i]) #check if the last value in each row is NA
    if (test[i] == 1 & check.NA == TRUE) { #checks whether there is only 1 NA in that column and whether this single NA is in the last row
      ff2[(rowFF2), i]<- ff2[(rowFF2-1), i] #replaces the last row (390pA) NA with the value of the previous row (380pA)
      ff2[(rowFF2+1), i]<- 999 #notes good cell
    }
    else if (test[i] >= 1 & check.NA == FALSE)    {
      stop("There is some problem with cell\t", colnames(ff2[i]), "\nNA values where they shouldn't be. Check if the data were entered correctly. \nAfter correcting, run the code from the beginning.")
    }
      else if (test[i]==0) {
      ff2[(rowFF2+1), i]<- 999 #notes good cell
    }
    else if (test[i] >1) {
      ff2[(rowFF2+1), i]<- 666 #notes bad cell)
      }  
}

#################################################################################
#################### SECTION 3 #################################################
########### select the good cells based on the firing frequency #################
#################################################################################

selectedCells<-which(ff2[nrow(ff2),]==999) #keep the columns of good cells in ff2
goodCells<-cbind(ff2[selectedCells]) #select only the good Cells
goodCells<-head(goodCells,-1) #remove the last unwanted row that we added before (the one with 666s and 999s)
goodCells.names <- c(colnames(goodCells)) # select the names of all the good cells ###maybe I have to remove sweep and injected
badCells.FF.names<-colnames(ff2)[which(ff2[nrow(ff2),]==666)] #keep the names of the bad cells (based on missing sweeps) for filtering the CCIO dataset


##################################################################################
################## SECTION 4 #####################################################
############ read the CCIO data ##################################################
####### exclude cells based on Erest, Rin, tau, cap z-scores #####################
##################################################################################


# read CCIO data for exclusion based on z-score#
data <- read.csv(Link2, header=T)

# Check data integrity #
headTail(data, 5,5)
dim(data)
names(data)



# clean data #
data[data == ""] <- NA #replace empty cells with NAs
all_na <- function(x) any(!is.na(x))
data <- data[!apply(is.na(data) | data == "", 1, all),] #clean rows where all columns are NA
row.names(data)<-c(1:nrow(data)) #reset the rownames to reflect the actual number of rownames
rownames(data) #check rownames
colnames(data) <- c("Cell_ID", "AP threshold (mV)", "AP amplitude (mV)", "AP half-width (ms)", "AHP amplitude (mV)", "Erest", "max voltage step", "voltage step", "V_sag", "% sag", "Rin", "current step tau (pA)", "tau", "rheobase (pA)", "capacitance", "firing mode","phenotype", "comments", "experimenter", "Georgia comments") # may need to change these names if i want to use a funtion
data$`% sag`<-replace(data$`% sag`, data$`% sag`<0,0) #replace all negative %sag with 0s

# check how many data have comment ""exclude", how many "interneuron", and which have missing values in Erest #
# Rin, tau and capacitance. All these cells will be excluded from the analysis #
excl<-which(data$comments == "exclude")
int<-which(data$comments == "interneuron")
misErest<-which(is.na(data$Erest))
misRin<-which(is.na(data$Rin))
mistau<-which(is.na(data$tau))
miscap<-which(is.na(data$capacitance))
length(excl)
length(int)
length(misErest)
length(misRin)
length(mistau)
length(miscap)
cat("Excluded cells:\n", length(excl), "based on the comments\n", length(int), "interneuron(s)\n", length(misErest), "based on missing Erest value\n", length(misRin), "based on missing Rin value\n", length(mistau), "based on missing tau value\n", length(miscap), "based on missing capacitance value\n")


data2<-subset(data, comments!="exclude" | is.na(comments)) #select all cells except comments "exclude"
data2<-subset(data2, comments !="interneuron" | is.na(comments)) #select all cells except interneurons


#function for exclusion of the cells that have NA values for Erest, Rin, tau and capacitance
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data2[completeVec, ])
}

data2<-completeFun(data2, c("Erest", "Rin", "tau", "capacitance")) #exclude the cells with NAs in Erest or Rin or tau or capacitance or firing mode#
if (all_na(data2$`firing mode`) ==TRUE) {
  data2<-completeFun(data2, "firing mode") #only if the column is not empty 9meaning there is a firing mode so P21 or adult), only then delete the rows that don't have a firing mode 
  
}

new_data <- data2[data2$Cell_ID %in% goodCells.names, ]#make new data variable based on the selected cells from the previous step (goodCells)
row.names(new_data)<-c(1:nrow(new_data)) #reset the rownames to reflect the actual number of rownames

#check if everything is okay#
check.cells.1 <-new_data$Cell_ID %in% goodCells.names #check whether all of the cells in new_data are contained in the selected cells MUST BE ALL TRUE
check.cells.2 <-new_data$Cell_ID %in% colnames(ff) #check whether all of the cells in new data are contained in the initial dataset MUST BE ALL TRUE
check.cells.3<-new_data$Cell_ID %in% badCells.FF.names #check whether the cells in the new_data contain any bad cells (FF based) MUST BE ALL FALSE
check.cells.4<-new_data$Cell_ID %in% badCells.AP.names #check whether the cells in the new_data contain any bad cells (AP based) MUST BE ALL FALSE

if (all(check.cells.1 == TRUE) & all(check.cells.2 == TRUE) & all(check.cells.3 == FALSE) & all(check.cells.4 == FALSE)) {
  cat("You can proceed to section 5")
} else {
    stop("There is some problem with the cell filtering. Check the cell IDs in the original files very carefully. \nMake sure the cells IDs in the Firing Frequency and CCIO files are exactly the same.\nRe-run section 4.")
}

###########################################################################
################### SECTION 5 #############################################
########### Finding outlying cells ########################################
########## Based on Erest, Rin, tau, capacitance ##########################
###########################################################################

# Shapiro test normality check #
filtering.variables<-c("Erest", "Rin", "tau", "capacitance")
for (i in filtering.variables)
{
  normTest<-stat.desc(new_data[i], basic=F, norm=T)
  if (normTest["normtest.p", i] < 0.05) {
    cat("\n", i, "is not normally distributed, SW = ", normTest["normtest.W", i], ", p-value = ", normTest["normtest.p", i], "\n")
  } else {
    cat("\n" ,i, "is normally distributed, SW = ", normTest["normtest.W", i], ", p-value = ", normTest["normtest.p", i], "\n")
  }
}  

# histograms #
dataLong <- reshape2::melt(new_data, id.vars = "Cell_ID", measure.vars = c("Erest", "Rin", "tau", "capacitance")) #create a long dataset

hist.all<-ggplot(dataLong, aes(x = value)) + 
  facet_wrap(~ variable, scales = "free", ncol = 2) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white",bins=5) + 
  theme_classic() +
  theme(text=element_text(size=15,family="Calibri Light",color="#2e2e2e")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)))

hist.all

# outliers with z-score use #
descriptives<-as.data.frame(psych::describe(new_data, na.rm = TRUE, omit = TRUE)) #calculate desccriptive statistics of the data

total.cells<-c()
total.rows<-c()

for (i in filtering.variables) {
  y<-paste("z.", i, sep = "")  
  dist.from.mean<-(abs(new_data[i] - descriptives[i,"mean"])) #calculate for every value the abs distance from the mean
  stdev<-descriptives[i, "sd"] #calculate the standard deviation
  new_data[y] <- dist.from.mean/stdev #calculate z-score
  out<-c(new_data[which(new_data[,y]>2),]) #values that should be excluded
  out.rows<-which(new_data[y]>2) #rows of values that should be excluded
  if (length(out.rows) != 0) {
    cat ("\n\nThe cells to be excluded based on the", i, "z-scores are:\n", out$Cell_ID, "in rows:", out.rows, "with", i, "value(s)", out[[i]], sep="  ")
  } else {
    cat("\n\nThere are no cells that should be excluded based on the", i, "z-scores")
  }
  total.cells<-c(total.cells, out$Cell_ID)
  total.rows<-c(total.rows, out.rows)
}

cat("You must exclude a total of", length(total.rows), "cells from this dataset. These are:\n", total.cells, "\nin rows:\n", total.rows, sep="\n") 

# boxplots - 1.5*IQR rule # this is just to look at the data outliers with the boxplot. should not be taken into account for outlier removal
overall.box<-ggplot(dataLong, aes(x = factor(variable), y = value, fill = factor(variable))) + 
  geom_boxplot() +
  stat_summary(aes(label = round(stat(y), 1)),geom = "text",fun.y = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },hjust = -1)

overall.box

# chosing the data without the outliers (outliers taken from z-score results) #
new_data$out<-new_data$Cell_ID %in% total.cells #note the outliers
final.CCIO.data <- new_data[!(new_data$out == TRUE),] #exclude the cells that are outliers in any of the Erest, Rin, tau, capacitance category
#### If the cells are from P21 or Adult layer 5/6, they have a firing mode. Wee need to separate ####
#### them based on their firing mode ####
final.CCIO.data<-final.CCIO.data[order(final.CCIO.data$`firing mode`),] #sort data based on firing mode
row.names(final.CCIO.data) <- c(1:nrow(final.CCIO.data))

final.ff.data<-goodCells %>% select(matches(final.CCIO.data$Cell_ID)) #select the non-outlier cells in the firing frequency dataset
final.ff.data<-cbind(goodCells["sweep"], goodCells["injected"], final.ff.data) #combine the sweep and injected current columns with the rest

#check whether the CCIO and FF files have the same N
N.cells.ff<-length(final.ff.data[,3:ncol(final.ff.data)])
N.cells.CCIO<-nrow(final.CCIO.data)
if (N.cells.ff == N.cells.CCIO) { # check whether CCIO and ff contain the same amount of cells
  N.cells<- N.cells.CCIO
  cat("Your N is:", N.cells, "cells.\n")
} else {
  stop("The number of cells in the CCIO file and firing frequency file is not the same. Something went wrong.\n")
}


#check if the final.CCIO.data contain the same cells as the final.ff.data#
cell.names.CCIO.dataset <- c(final.CCIO.data$Cell_ID)
cell.names.ff.dataset<-colnames(final.ff.data[3:ncol(final.ff.data)])
if (setequal(cell.names.CCIO.dataset, cell.names.ff.dataset) == FALSE) {
  stop("The cells contained in the CCIO dataset are not the same as in the Firing frequency dataset.\nSomething went wrong.\n")
} else if (setequal(cell.names.CCIO.dataset, cell.names.ff.dataset) == TRUE) {
  cat("You can proceed to section 6\n")
}


#########################################################################
##################### SECTION 6 #########################################
######### Normality check after the exclusion of outliers ###############
#########Save the new clean data ########################################
#########################################################################

#normality check without the outliers - check if now they are normal
# Shapiro test normality check #
filtering.variables<-c("Erest", "Rin", "tau", "capacitance")
for (i in filtering.variables)
{
  normTest<-stat.desc(final.CCIO.data[i], basic=F, norm=T)
  if (normTest["normtest.p", i] < 0.05) {
    cat("\n\n", i, "is not normally distributed, SW = ", normTest["normtest.W", i], ", p-value = ", normTest["normtest.p", i])
  } else {
    cat("\n\n" ,i, "is normally distributed, SW = ", normTest["normtest.W", i], ", p-value = ", normTest["normtest.p", i])
  }
}  

## histograms ##
dataLong.final <- reshape2::melt(final.CCIO.data, id.vars = "Cell_ID", measure.vars = c("Erest", "Rin", "tau", "capacitance")) #create a long dataset

hist.all<-ggplot(dataLong.final, aes(x = value)) + 
  facet_wrap(~ variable, scales = "free", ncol = 2) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white",bins=5) + 
  theme_classic() +
  theme(text=element_text(size=15,family="Calibri Light",color="#2e2e2e")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)))

hist.all
#ggsave("hist.all.png", hist.all) 
 
# write clean data # set names as precious + the word FINAL in the end
cat("\nSave the Firing frequency clean data\n")
write.csv(final.ff.data, file=file.choose(), row.names=FALSE) #set the filename - in create file press yes
cat("\nSave the CCIO clean data\n")
write.csv(final.CCIO.data, file=file.choose(), row.names=FALSE)


##########################################################################
######### at this point the output files are ready #######################
######## for graph making in graphpad ###################################
############################################################################

#############################################################################
############## NEXT STEP STATISTICS ######################################
##########################################################################



#########################################################
######### STATISTICS ACTIVE PASSIVE PROPERTIES ##########
############### FACTORIAL ANOVA #########################
#########################################################

#########################################################
############### DESCRIPTIVE STATISTICS ##################
#########################################################


fileLinks.ages<-data.frame()
mouse.ages<-c("P7", "P14", "P21", "Adult")
si.conditions<-c("control", "KD")

for (age in mouse.ages) {
  
  for (condition in si.conditions){
  
  cat("\nChoose the" , age, condition,  "FINAL CCIO file")
 
  fileLinks.ages[age,condition] <- file.choose() # opens a window to choose the files
  
  }
}

setwd(tk_choose.dir())
getwd()

# read data #
P7.Control<-read.csv(fileLinks.ages[1,1], header = TRUE, check.names = FALSE)
P7.KD<-read.csv(fileLinks.ages[1,2], header = TRUE, check.names = FALSE)
P14.Control<-read.csv(fileLinks.ages[2,1], header = TRUE, check.names = FALSE)
P14.KD<-read.csv(fileLinks.ages[2,2], header = TRUE, check.names = FALSE)
P21.Control<-read.csv(fileLinks.ages[3,1], header = TRUE, check.names = FALSE)
P21.KD<-read.csv(fileLinks.ages[3,2], header = TRUE, check.names = FALSE)
adult.Control<-read.csv(fileLinks.ages[4,1], header = TRUE, check.names = FALSE)
adult.KD<-read.csv(fileLinks.ages[4,2], header = TRUE, check.names = FALSE)

# ading the age variable ### this should be added previously or when the data were composed
P7.Control$age<-"P7"
P7.KD$age<-"P7"
P14.Control$age<-"P14"
P14.KD$age<-"P14"
P21.Control$age<-"P21"
P21.KD$age<-"P21"
adult.Control$age<-"adult"
adult.KD$age<-"adult"

# combine everything into one dataset #
complete.dataset<-rbind.data.frame(P7.Control, P7.KD, P14.Control, P14.KD, P21.Control, P21.KD, adult.Control, adult.KD)

# Make age and phenotype factors and adjust the orders of their levels
complete.dataset$age <- factor(complete.dataset$age, levels = c("P7","P14","P21","adult"))
complete.dataset$phenotype <- factor(complete.dataset$phenotype, levels = c("siControl","siCB1"))

# Check new order
print(levels(complete.dataset$age))
print(levels(complete.dataset$phenotype))

# separate the data into 8 different groups (4 ages x 2 siConditions)
complete.dataset$group<-ifelse(complete.dataset$age == "P7" & complete.dataset$phenotype == "siControl", "group1", 
                               ifelse(complete.dataset$age == "P7" & complete.dataset$phenotype == "siCB1", "group2",
                               ifelse(complete.dataset$age == "P14" & complete.dataset$phenotype == "siControl", "group3",
                               ifelse(complete.dataset$age == "P14" & complete.dataset$phenotype == "siCB1", "group4",
                               ifelse(complete.dataset$age == "P21" & complete.dataset$phenotype == "siControl", "group5",
                               ifelse(complete.dataset$age == "P21" & complete.dataset$phenotype == "siCB1", "group6",
                               ifelse(complete.dataset$age == "adult" & complete.dataset$phenotype == "siControl", "group7",
                               ifelse(complete.dataset$age == "adult" & complete.dataset$phenotype == "siCB1", "group8", NA))))))))



#keep only the columns we are going to test#
statistics.dataset<-cbind.data.frame(Cell_ID = complete.dataset$Cell_ID, 
                                     group = complete.dataset$group, 
                                     age = complete.dataset$age, 
                                     phenotype = complete.dataset$phenotype, 
                                     APthreshold = complete.dataset$`AP threshold (mV)`, 
                                     APamplitude = complete.dataset$`AP amplitude (mV)`, 
                                     APhalf_width = complete.dataset$`AP half-width (ms)`, 
                                     AHPamplitude = complete.dataset$`AHP amplitude (mV)`,
                                     sag = complete.dataset$`% sag`,
                                     Erest = complete.dataset$Erest,
                                     Rin = complete.dataset$Rin, 
                                     tau = complete.dataset$tau, 
                                     capacitance = complete.dataset$capacitance, 
                                     stringsAsFactors = TRUE)

# check if the character columns are factors
is.factor(statistics.dataset$age)
is.factor(statistics.dataset$phenotype)
is.factor(statistics.dataset$group)
is.factor(statistics.dataset$APthreshold)

#make group a factor #
statistics.dataset$group<-as.factor(statistics.dataset$group)

#check the levels of the factors in the AP dataset
levels(statistics.dataset$age)
levels(statistics.dataset$phenotype)
levels(statistics.dataset$group)

# check assumptions #

# Create separate goup data objects
list2env(split(statistics.dataset, statistics.dataset$group), envir = .GlobalEnv)

# TEMPORARY solution because group 1 has N=2, we have to exclude it
statistics.dataset<-subset(statistics.dataset, group!="group1" & group!="group2") #select all groups except 1 and 2 (P7)


#function for testing homogeneity of variance #
homogeneity.variances<-function(x, y){
  homo<-leveneTest(x, as.factor(statistics.dataset$group), center = mean)
  print(homo)
  if (homo$`Pr(>F)`[1] > 0.05) {
    cat("The assumption of homogeneity of variances for", y, "is satisfied, levene test p value =", homo$`Pr(>F)`[1], "\n\n")
  } else {
    cat("The assumption of homogeneity of variance for", y, "is violated, levene test p value =", homo$`Pr(>F)`[1], "\n\n")
  }
}


homogeneity.variances(statistics.dataset$APthreshold, "APthreshold")
homogeneity.variances(statistics.dataset$APamplitude, "APamplitude")
homogeneity.variances(statistics.dataset$APhalf_width, "APhalf_width")
homogeneity.variances(statistics.dataset$AHPamplitude, "AHPamplitude")
homogeneity.variances(statistics.dataset$sag, "sag")
homogeneity.variances(statistics.dataset$Erest, "Erest")
homogeneity.variances(statistics.dataset$Rin, "Rin")
homogeneity.variances(statistics.dataset$tau, "tau")
homogeneity.variances(statistics.dataset$capacitance, "capacitance")

#  homogeneity of variances test and normality tests for each group separately #
sink("active passive properties descriptive statistics.txt")

cat("#########################################################################\n############## Levene test for homogeneity of variances  ################\n#########################################################################\n\n")
homogeneity.variances(statistics.dataset$APthreshold, "APthreshold")
homogeneity.variances(statistics.dataset$APamplitude, "APamplitude")
homogeneity.variances(statistics.dataset$APhalf_width, "APhalf_width")
homogeneity.variances(statistics.dataset$AHPamplitude, "AHPamplitude")
homogeneity.variances(statistics.dataset$sag, "sag")
homogeneity.variances(statistics.dataset$Erest, "Erest")
homogeneity.variances(statistics.dataset$Rin, "Rin")
homogeneity.variances(statistics.dataset$tau, "tau")
homogeneity.variances(statistics.dataset$capacitance, "capacitance")


cat("#########################################################################\n############# Descriptive statistics and normality check ################\n#########################################################################\n\n")
descriptives.by.group<-by(statistics.dataset, statistics.dataset$group, stat.desc, basic = TRUE, norm=TRUE)
descriptives.by.group

sink()

#create long datasets and histograms for each group#

hist.function<-function(x,y){
  group.dataLong<-reshape2::melt(x, id.vars = "Cell_ID", measure.vars = c("APthreshold", "APamplitude", "APhalf_width", "AHPamplitude", "sag", "Erest", "Rin", "tau", "capacitance"))
  hist.group<-ggplot(group.dataLong, aes(x = value)) + 
    facet_wrap(~ variable, scales = "free", ncol = 3) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white",bins=5) + 
    theme_classic() +
    ggtitle(y) 
  print(hist.group)
}

hist.function(group2, "group1") #don't use it in statistics for now
hist.function(group2, "group2")  #don't use it in statistics for now
hist.function(group3, "group3")
hist.function(group4, "group4")
hist.function(group5, "group5")
hist.function(group6, "group6")
hist.function(group7, "group7")
hist.function(group8, "group8")





#########################################################
############### FACTORIAL ANOVA #########################
#########################################################

my.anova.model.function<-function(x,y, z){
  cat(z)
  anova.model<-lm(x~age + phenotype + age:phenotype, data = y)
  car::Anova(anova.model, type=3)
}

sink("active passive properties factorial anova.txt")
print(my.anova.model.function(statistics.dataset$APthreshold, statistics.dataset, "\n\nAPthreshold\n"))
print(my.anova.model.function(statistics.dataset$APamplitude, statistics.dataset, "\n\nAPamplitude\n"))
print(my.anova.model.function(statistics.dataset$APhalf_width, statistics.dataset, "\n\nAPhalf-width\n"))
print(my.anova.model.function(statistics.dataset$AHPamplitude, statistics.dataset, "\n\nAHPamplitude\n"))
print(my.anova.model.function(statistics.dataset$sag, statistics.dataset, "\n\nsag\n"))
print(my.anova.model.function(statistics.dataset$Erest, statistics.dataset, "\n\nErest\n"))
print(my.anova.model.function(statistics.dataset$Rin, statistics.dataset, "\n\nRin\n"))
print(my.anova.model.function(statistics.dataset$tau, statistics.dataset, "\n\ntau\n"))
print(my.anova.model.function(statistics.dataset$capacitance, statistics.dataset, "\n\ncapacitance\n"))
sink()

#Calculate the adjusted a
# new a = a / number of tests conducted #
0.05/9

statistics.dataset %>% 
  ggplot() +
  aes(x = age, y = APthreshold, color = phenotype) +
  geom_line(aes(group = phenotype)) 




# Interaction plots #

#function#
plots.interaction.function<-function(dataset, x.axon.variable, y.axon.variable, linecolor.variable, group.variable, x.axon.title, y.axon.title) {
  plots.interaction <- ggplot(dataset) +
    aes(x = x.axon.variable, y = y.axon.variable, color = linecolor.variable, group = group.variable) +
    stat_summary(fun.y=mean,geom="point", size = 2) +
    stat_summary(fun.y=mean, geom="line", size = 0.7) +
    stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=.2) +
    labs(x=x.axon.title, y=y.axon.title)
  
  plots.interaction
}

plots.interaction.function(statistics.dataset, 
                           statistics.dataset$phenotype, 
                           statistics.dataset$APthreshold, 
                           statistics.dataset$age, 
                           statistics.dataset$age, 
                           "", 
                           "AP threshold")


plots.interaction.function(statistics.dataset, 
                           statistics.dataset$phenotype, 
                           statistics.dataset$APamplitude, 
                           statistics.dataset$age, 
                           statistics.dataset$age, 
                           "", 
                           "AP amplitude")


plots.interaction.function(statistics.dataset, 
                           statistics.dataset$phenotype, 
                           statistics.dataset$APhalf_width, 
                           statistics.dataset$age, 
                           statistics.dataset$age, 
                           "", 
                           "AP half-width")


plots.interaction.function(statistics.dataset, 
                           statistics.dataset$phenotype, 
                           statistics.dataset$AHPamplitude, 
                           statistics.dataset$age, 
                           statistics.dataset$age, 
                           "", 
                           "AHP amplitude")



plots.interaction.function(statistics.dataset, 
                           statistics.dataset$phenotype, 
                           statistics.dataset$sag, 
                           statistics.dataset$age, 
                           statistics.dataset$age, 
                           "", 
                           "% sag")



plots.interaction.function(statistics.dataset, 
                           statistics.dataset$phenotype, 
                           statistics.dataset$Erest, 
                           statistics.dataset$age, 
                           statistics.dataset$age, 
                           "", 
                           "Erest")

plots.interaction.function(statistics.dataset, 
                           statistics.dataset$phenotype, 
                           statistics.dataset$Rin, 
                           statistics.dataset$age, 
                           statistics.dataset$age, 
                           "", 
                           "Rin")


plots.interaction.function(statistics.dataset, 
                           statistics.dataset$phenotype, 
                           statistics.dataset$tau, 
                           statistics.dataset$age, 
                           statistics.dataset$age, 
                           "", 
                           "tau")


plots.interaction.function(statistics.dataset, 
                           statistics.dataset$phenotype, 
                           statistics.dataset$capacitance, 
                           statistics.dataset$age, 
                           statistics.dataset$age, 
                           "", 
                           "capacitance")



#########################################################
############## STATISTICS FIRING FREQUENCY ##############
######## 2-WAY REPEATED MEASURES ANOVA ##################
#########################################################

#########################################################
############### DESCRIPTIVE STATISTICS ##################
#########################################################

# choose directory
setwd(tk_choose.dir())
getwd()

read.ff.data.function<-function(age.and.phenotype) {

  injected<-seq(from = -100, to = 390, by=10)
  injected<-paste(rep("inj", 50), injected, rep("pA", 50))
  column.names.wide<-c("Cell_ID", "age", "phenotype", injected)

  if (age.and.phenotype=="P7.Control") {
    
    P7.Control<-read.csv(fileLinks.ages["P7","siControl"], header = TRUE, check.names = FALSE)
    P7.Control<-P7.Control[,-c(1:2)] # throw out the "sweep" and injected column
    N.cells.P7.control<-ncol(P7.Control) #register the number of cells 
    Cell_ID<-colnames(P7.Control) #select the cell ID columns
    wide.data<-transpose(P7.Control) #transpose the data 
    wide.data<-cbind.data.frame(Cell_ID, age,phenotype,wide.data) #cbind the Cell id, age and phenotype columns
    names(wide.data)<-column.names.wide #add column names

  } else if (age.and.phenotype=="P7.KD") {
    
    P7.KD<-read.csv(fileLinks.ages["P7","siCB1"], header = TRUE, check.names = FALSE)
    P7.KD<-P7.KD[,-c(1:2)]
    N.cells.P7.KD<-ncol(P7.KD) 
    Cell_ID<-colnames(P7.KD) 
    wide.data<-transpose(P7.KD)
    wide.data<-cbind.data.frame(Cell_ID, age,phenotype,wide.data)
    names(wide.data)<-column.names.wide
    
  } else if (age.and.phenotype=="P14.Control") {
    
    P14.Control<-read.csv(fileLinks.ages["P14","siControl"], header = TRUE, check.names = FALSE)
    P14.Control<-P14.Control[,-c(1:2)] 
    N.cells.P14.Control<-ncol(P14.Control)
    Cell_ID<-colnames(P14.Control)
    wide.data<-transpose(P14.Control)
    wide.data<-cbind.data.frame(Cell_ID, age,phenotype,wide.data)
    names(wide.data)<-column.names.wide
    
  } else if (age.and.phenotype=="P14.KD") {
    
    P14.KD<-read.csv(fileLinks.ages["P14","siCB1"], header = TRUE, check.names = FALSE)
    P14.KD<-P14.KD[,-c(1:2)]
    N.cells.P14.KD<-ncol(P14.KD) 
    Cell_ID<-colnames(P14.KD)
    wide.data<-transpose(P14.KD)
    wide.data<-cbind.data.frame(Cell_ID, age,phenotype,wide.data)
    names(wide.data)<-column.names.wide
    
  } else if (age.and.phenotype=="P21.Control") {
    
    P21.Control<-read.csv(fileLinks.ages["P21","siControl"], header = TRUE, check.names = FALSE)
    P21.Control<-P21.Control[,-c(1:2)] 
    N.cells.P21.Control<-ncol(P21.Control)
    Cell_ID<-colnames(P21.Control)
    wide.data<-transpose(P21.Control)
    wide.data<-cbind.data.frame(Cell_ID, age,phenotype,wide.data)
    names(wide.data)<-column.names.wide
    
  } else if (age.and.phenotype=="P21.KD") {
    
    P21.KD<-read.csv(fileLinks.ages["P21","siCB1"], header = TRUE, check.names = FALSE)
    P21.KD<-P21.KD[,-c(1:2)] 
    N.cells.P21.KD<-ncol(P21.KD)
    Cell_ID<-colnames(P21.KD) 
    wide.data<-transpose(P21.KD)
    wide.data<-cbind.data.frame(Cell_ID, age,phenotype,wide.data)
    names(wide.data)<-column.names.wide
    
  } else if (age.and.phenotype=="adult.Control") {
    
    adult.Control<-read.csv(fileLinks.ages["adult","siControl"], header = TRUE, check.names = FALSE)
    adult.Control<-adult.Control[,-c(1:2)] 
    N.cells.adult.Control<-ncol(adult.Control)
    Cell_ID<-colnames(adult.Control)
    wide.data<-transpose(adult.Control)
    wide.data<-cbind.data.frame(Cell_ID, age,phenotype,wide.data)
    names(wide.data)<-column.names.wide
    
  }  else if (age.and.phenotype=="adult.KD") {
    
    adult.KD<-read.csv(fileLinks.ages["adult","siCB1"], header = TRUE, check.names = FALSE)
    adult.KD<-adult.KD[,-c(1:2)] 
    N.cells.adult.KD<-ncol(adult.KD)
    Cell_ID<-colnames(adult.KD) 
    wide.data<-transpose(adult.KD)
    wide.data<-cbind.data.frame(Cell_ID, age,phenotype,wide.data)
    names(wide.data)<-column.names.wide
    
  }
  return(wide.data)
}


#read data
fileLinks.ages<-data.frame()
mouse.ages<-c("P7", "P14", "P21", "adult")
phenotypes<-c("siControl", "siCB1")
injected<-seq(from = -100, to = 390, by=10)
injected<-paste(rep("inj", 50), injected, rep("pA", 50))
column.names.wide<-c("Cell_ID", "age", "phenotype", injected)



for (age in mouse.ages) {
  
  for (phenotype in phenotypes){
    
    cat("\nChoose the" , age, phenotype,  "FINAL Firing Frequency file")
    
    fileLinks.ages[age,phenotype] <- file.choose() # opens a window to choose the files
    if (phenotype == "siControl" & age == "P7") {
      wide.P7.Control<-read.ff.data.function("P7.Control")
    } else if (phenotype == "siCB1" & age == "P7"){
      wide.P7.KD<-    read.ff.data.function("P7.KD")
    } else if (phenotype == "siControl" & age == "P14"){
      wide.P14.Control<-read.ff.data.function("P14.Control")
    } else if (phenotype == "siCB1" & age == "P14") {
      wide.P14.KD<-read.ff.data.function("P14.KD")
    } else if (phenotype == "siControl" & age == "P21"){
      wide.P21.Control<-read.ff.data.function("P21.Control")
    } else if (phenotype == "siCB1" & age == "P21") {
      wide.P21.KD<-read.ff.data.function("P21.KD")
    } else if (phenotype == "siControl" & age == "adult") {
      wide.adult.Control<-read.ff.data.function("adult.Control")
    } else if (phenotype == "siCB1" & age == "adult") {
      wide.adult.KD<-read.ff.data.function("adult.KD")
    } 
  }
}


# read original data as well - good check to compare them #
#P7.Control<-read.csv(fileLinks.ages[1,1], header = TRUE, check.names = FALSE)
#P7.KD<-read.csv(fileLinks.ages[1,2], header = TRUE, check.names = FALSE)
#P14.Control<-read.csv(fileLinks.ages[2,1], header = TRUE, check.names = FALSE)
#P14.KD<-read.csv(fileLinks.ages[2,2], header = TRUE, check.names = FALSE)
#P21.Control<-read.csv(fileLinks.ages[3,1], header = TRUE, check.names = FALSE)
#P21.KD<-read.csv(fileLinks.ages[3,2], header = TRUE, check.names = FALSE)
#adult.Control<-read.csv(fileLinks.ages[4,1], header = TRUE, check.names = FALSE)
#adult.KD<-read.csv(fileLinks.ages[4,2], header = TRUE, check.names = FALSE)

#combine all data together and make phenotype and age as factors & make groups based on age and phenotype
wide.ff.data<-rbind.data.frame(wide.P7.Control, wide.P7.KD, wide.P14.Control, wide.P14.KD, wide.P21.Control, wide.P21.KD, wide.adult.Control, wide.adult.KD)

wide.ff.data$group<-ifelse(wide.ff.data$age == "P7" & wide.ff.data$phenotype == "siControl", "P7 siControl", 
                    ifelse(wide.ff.data$age == "P7" & wide.ff.data$phenotype == "siCB1", "P7 siCB1",
                    ifelse(wide.ff.data$age == "P14" & wide.ff.data$phenotype == "siControl", "P14 siControl",
                    ifelse(wide.ff.data$age == "P14" & wide.ff.data$phenotype == "siCB1", "P14 siCB1",
                    ifelse(wide.ff.data$age == "P21" & wide.ff.data$phenotype == "siControl", "P21 siControl",
                    ifelse(wide.ff.data$age == "P21" & wide.ff.data$phenotype == "siCB1", "P21 siCB1",
                    ifelse(wide.ff.data$age == "adult" & wide.ff.data$phenotype == "siControl", "adult siControl",
                    ifelse(wide.ff.data$age == "adult" & wide.ff.data$phenotype == "siCB1", "adult siCB1", NA))))))))

#set the order of the levels of the factors
wide.ff.data$age <- factor(wide.ff.data$age, levels = c("P7","P14","P21","adult"))
wide.ff.data$phenotype <- factor(wide.ff.data$phenotype, levels = c("siControl","siCB1"))
wide.ff.data$group <- factor(wide.ff.data$group, levels = c("P7 siControl","P7 siCB1", "P14 siControl", "P14 siCB1", "P21 siControl", "P21 siCB1", "adult siControl", "adult siCB1"))
wide.ff.data <- wide.ff.data %>%  select(Cell_ID, age, phenotype, group, everything()) #re order the columns so group column is at the beginning
#wide.ff.data[,5:54] <- lapply(wide.ff.data[,5:54], function(x) as.numeric(as.character(x))) #make the columns of the injected current values numeric
class(wide.ff.data$`inj -90 pA`)

# Check new order of levels of factors
print(levels(wide.ff.data$age))
print(levels(wide.ff.data$phenotype))
print(levels(wide.ff.data$group))


#print a table to see how many cells we have in each category
table(wide.ff.data$age, wide.ff.data$phenotype)


#make a long dataset

long.ff.data <- reshape2::melt(wide.ff.data, id.vars = c("Cell_ID", "age", "phenotype", "group"), measure.vars = injected) #create a long dataset
colnames(long.ff.data)<-c("Cell_ID", "age", "phenotype", "group","injected.current", "firing.frequency")
print(levels(long.ff.data$age))
print(levels(long.ff.data$phenotype))
print(levels(long.ff.data$group))
class(long.ff.data$`firing.frequency`)
class(long.ff.data$injected.current)
levels(long.ff.data$injected.current)


# make line plot
lineff <- ggplot(long.ff.data,aes(injected.current,firing.frequency, colour=group)) +
  stat_summary(fun=mean,geom="point", size=1) +
  stat_summary(fun=mean, geom="line", aes(group=group), linetype="solid", size=1) +
  stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=.2) +
  theme_classic()+
  labs(x="Injected current (pA)", y="Mean firing frequency (Hz)")

grid.draw(lineff) # interactive device 


#make datasets for normality tests - wide format and filter out any columns where all values are zero - we can't perform normlity tests there
norm.P7.Control<-Filter(function(x) !(all(x==0)), wide.P7.Control)
norm.P7.KD<-Filter(function(x) !(all(x==0)), wide.P7.KD)
norm.P14.Control<-Filter(function(x) !(all(x==0)), wide.P14.Control)
norm.P14.KD<-Filter(function(x) !(all(x==0)), wide.P14.KD)
norm.P21.Control<-Filter(function(x) !(all(x==0)), wide.P21.Control)
norm.P21.KD<-Filter(function(x) !(all(x==0)), wide.P21.KD)
norm.adult.Control<-Filter(function(x) !(all(x==0)), wide.adult.Control)
norm.adult.KD<-Filter(function(x) !(all(x==0)), wide.adult.KD)


# subset the data in a long format based on ages to do separate mixed repeated measures ANOVAs (two-way RMA)
groupP7<-subset(long.ff.data, age=="P7")
groupP14<-subset(long.ff.data, age=="P14")
groupP21<-subset(long.ff.data, age=="P21")
groupAdult<-subset(long.ff.data, age=="adult")


#function for levene test within each injected current 
levene.function<-function(dataset.Control, dataset.KD, age){
  
  data<-rbind.data.frame(dataset.Control, dataset.KD)
  
  cat("###########", age, "###########")
  for (i in 4:53) {
    cat("\n\n---------------------------------------------------------------------")
    cat("\n", names(data[i]))
    cat("\n---------------------------------------------------------------------\n\n")
    print(leveneTest(data[,i], as.factor(data$phenotype), center = mean))
  }
}

#Box's test for Homogeneity of covariances matrix function
box.function<-function(dataset, age) {
  box.results<-box_m(dataset[, "firing.frequency", drop = FALSE], dataset$phenotype)
  cat("\n\n---------------------------------------------------------------------\n")
  cat(age)
  cat("\n---------------------------------------------------------------------\n\n")
  print(box.results)
}

# write file with descriptive statistics
sink("firing frequency descriptive statistics.txt")
  cat("#########################################################################\n############# Descriptive statistics and normality check ################\n#########################################################################\n\n")
  cat("\n\nP7 siControl\n\n")
  stat.desc(norm.P7.Control, basic=TRUE, norm = TRUE)
  cat("\n\nP7 siCB1\n\n")
  stat.desc(norm.P7.KD, basic=TRUE, norm = TRUE)
  cat("\n\nP14 siControl\n\n")
  stat.desc(norm.P14.Control, basic=TRUE, norm = TRUE)
  cat("\n\nP14 siCB1\n\n")
  stat.desc(norm.P14.KD, basic=TRUE, norm = TRUE)
  cat("\n\nP21 siControl\n\n")
  stat.desc(norm.P21.Control, basic=TRUE, norm = TRUE)
  cat("\n\nP21 siCB1\n\n")
  stat.desc(norm.P21.KD, basic=TRUE, norm = TRUE)
  cat("\n\nadult siControl\n\n")
  stat.desc(norm.adult.Control, basic=TRUE, norm = TRUE)
  cat("\n\nadult siCB1\n\n")
  stat.desc(norm.adult.KD, basic=TRUE, norm = TRUE)
  cat("\n\n#########################################################################\n############## Levene test for homogeneity of variances  ################\n#########################################################################\n\n")
  levene.function(wide.P7.Control, wide.P7.KD, "P7")
  levene.function(wide.P14.Control, wide.P14.KD, "P14")
  levene.function(wide.P21.Control, wide.P21.KD, "P21")
  levene.function(wide.adult.Control, wide.adult.KD, "adult")
  cat("\n\n#########################################################################\n########### Box's test for Homogeneity of covariances matrix  ############\n#########################################################################\n\n")
  box.function(groupP7, "P7")
  box.function(groupP14, "P14")
  box.function(groupP21, "P21")
  box.function(groupAdult, "adult")
sink()


#########################################################
######## 2-WAY REPEATED MEASURES ANOVA ##################
#########################################################

#function for conducting repeated measures ANOVA
### sphericity cannot be calculated because we have so many within levels compared to Cells
my.RMA.function<-function(dataset, age) {
  ffmod <-aov_ez("Cell_ID", "firing.frequency", dataset, 
                    between = c("phenotype"), within = c("injected.current"), 
                    anova_table=list(correction = "GG", es = "pes"), print.formula=T, 
                    type=c("3"), check_contrasts=T, return = c("afex_aov"))
  cat("\n\n---------------------------------------------------------------------\n")
  cat(age)
  cat("\n---------------------------------------------------------------------\n\n")
  print(get_anova_table(ffmod))
  cat("\n---------------------------------------------------------------------\n")
  summary(ffmod)
}

sink("firing frequency two-way repeated measures ANOVA.txt")
  my.RMA.function(groupP7, "P7")
  my.RMA.function(groupP14, "P14")
  my.RMA.function(groupP21, "P21")
  my.RMA.function(groupAdult, "adult")
sink()


# LINEAR MIXED MODELS
library(lme4)
lmeModel = lmer(firing.frequency ~ phenotype*injected.current + (1|Cell_ID), data=groupP21)
anova(lmeModel)

# make line plot
lineff <- ggplot(groupAdult,aes(injected.current,firing.frequency, colour=phenotype)) +
  stat_summary(fun=mean,geom="point", size=1) +
  stat_summary(fun=mean, geom="line", aes(group=phenotype), linetype="solid", size=1) +
  stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=.2) +
  theme_classic()+
  labs(x="Injected current (pA)", y="Mean firing frequency (Hz)")+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("P14")
  

grid.draw(lineff) # interactive device 

