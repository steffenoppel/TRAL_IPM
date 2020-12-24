############################################################################
######## DATA PREPARATION FOR INTEGRATED POPULATION MODEL     ##############
############################################################################

### written by Steffen Oppel in December 2018
### steffen.oppel@rspb.org.uk
### uses existing CMR and breeding database to extract data

## updated on 31 March 2020

library(tidyverse)
library(lubridate)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select


#############################################################################
##   1. SPECIFY THE SPECIES AND START YEAR FOR WHICH YOU WANT A SUMMARY ####
#############################################################################
SP<-"TRAL"
start<-2000




###################################################################################
##   2. READ IN DATA FROM DATABASES AND FILTER DATA FOR SPECIES OF INTEREST ####
###################################################################################

## run the RODBC import of nest and count data in a 32-bit version of R
system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database\\RODBC_count_import.r")), wait = TRUE, invisible = FALSE, intern = T)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database"), silent=T)
load("GOUGH_seabird_data.RData")

## run the RODBC import of CMR data in a 32-bit version of R
system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival\\RODBC_CMR_import.R")), wait = TRUE, invisible = FALSE, intern = T)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival"), silent=T)
load("GOUGH_seabird_CMR_data.RData")

## filter data for the selected species
contacts<-contacts %>% filter(SpeciesCode==SP)
ages<-ages %>% filter(SpeciesCode==SP)
nests<-nests %>% filter(Species==SP)
counts<-counts %>% filter(Species==SP)


## look at data
head(nests)  ## nest monitoring data
head(counts)  ## seabird count data



### CHECK DATA
unique(nests$Species)
unique(nests$Colony)
unique(nests$StageFound)
unique(nests$LastStage)

unique(counts$Species)
unique(counts$Colony)
unique(counts$Breed_Stage)
unique(counts$Cohort)



#############################################################################
##   3. PREPARE THE BREEDING SUCCESS DATA FROM NEST RECORDS #################
#############################################################################
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM"), silent=T)


### remove only partially monitored nests
exclude <- nests %>%
  filter(LastStage=="INCU") %>%
  filter(SUCCESS==1)


### summary of breeding success per year from nests
FECUND<-nests %>% filter(Year<2020) %>% filter(Species==SP) %>% mutate(count=1) %>%
  filter(!NestID %in% exclude$NestID) %>%
  group_by(Year,Colony) %>%
  summarise(n_nests=sum(count),BREED_SUCC=mean(SUCCESS, na.rm=T))

### inspect the dubious summaries
FECUND %>% filter(BREED_SUCC==0)		## Tafelkop 2011 appears to have no monitoring data, data from 2019 incomplete
FECUND %>% filter(BREED_SUCC==1)		## Green Hill North 2007 - only one check in August - worthless

FECUND<-FECUND %>% filter(Colony=="Gonydale")


### PLOT TO SPOT ANY OUTLIERS OF BREEDING SUCCESS
ggplot(FECUND, aes(x=Year,y=BREED_SUCC)) +geom_point(size=2, color='darkred')+geom_smooth(method='lm') 
fwrite(FECUND,"TRAL_breed_success_2006_2019.csv")



#############################################################################
##   4. PREPARE THE POPULATION COUNT DATA FROM COUNT RECORDS ################
#############################################################################
#try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdBreedingSuccess"), silent=T)

head(counts)

### summary of population counts of breeding pairs per year and colony
POPSIZE<-counts %>% filter(Species==SP) %>%
  mutate(Colony= as.character(Colony)) %>%
  mutate(Colony= if_else(Colony=="Green Hill South","Green Hill",Colony)) %>%
  mutate(Colony= if_else(Colony=="Green Hill North","Green Hill",Colony)) %>%
  mutate(Year=year(Date)) %>%
  filter(Year>start) %>%
  filter(Breed_Stage=="INCU") %>%
  filter(Cohort %in% c("INCU","TERR","AON")) %>%
  filter(!(Colony=='Albatross Plain'& Year==2009 & Method=='COLO')) %>%
  group_by(Year,Colony) %>%
  summarise(N=sum(Number, na.rm=T)) %>%
  spread(key=Colony, value=N)
POPSIZE



### summary of population counts of fledglings per year and colony
CHICKCOUNT<-counts %>% filter(Species==SP) %>%
  mutate(Colony= as.character(Colony)) %>%
  filter(!(Colony=="Green Hill North" & year(Date)==2017)) %>%
  filter(!(Colony=="Gonydale" & year(Date)==2017 & month(Date)==8)) %>%     ### remove the double count in 2017
  mutate(Colony= if_else(Colony=="Green Hill South","Green Hill",Colony)) %>%
  mutate(Colony= if_else(Colony=="Green Hill North","Green Hill",Colony)) %>%
  mutate(Year=year(Date)) %>%
  filter(Year>start) %>%
  filter(Breed_Stage %in% c("CHIC","FLED")) %>%
  filter(Cohort %in% c("CHIC","FLED")) %>%
  group_by(Year,Colony) %>%
  summarise(N=sum(Number, na.rm=T)) %>%
  spread(key=Colony, value=N)
CHICKCOUNT



### MAKE SURE BOTH MATRICES HAVE IDENTICAL DIMENSIONS
## no adult counts in 2002 and 2003
## no chick counts in 2013 and 2020
dim(POPSIZE)
dim(CHICKCOUNT)

POPSIZE <- POPSIZE %>% full_join(CHICKCOUNT[,1], by="Year") %>%
		arrange(Year)

CHICKCOUNT <- CHICKCOUNT %>% full_join(POPSIZE[,1], by="Year") %>%
		arrange(Year)

BS<-CHICKCOUNT/POPSIZE

counts %>% filter(Colony=='Gonydale' & year(Date)==2017)

dim(POPSIZE)
dim(CHICKCOUNT)
fwrite(CHICKCOUNT,"TRAL_CHIC_counts_2001_2020.csv")
fwrite(POPSIZE,"TRAL_INCU_counts_2001_2020.csv")




#############################################################################
##   5. PREPARE THE MARK-RECAPTURE DATA FOR SURVIVAL ANALYSIS ###############
#############################################################################

### checked on 25 Dec 2018: CMR database does not have satisfactory level of age or breeding status assignment for vast majority of contacts
## including 'age' and 'breeding status' reduces number of contacts from ~49000 to ~6000
## including the side on which a federal band was applied reduces the AYNA contacts from 25201 to 7786
## the only AYNA ringed as 'Chick' are from the 2015-16 season
## given the gross inadequacies of past records it is safer to simply use 0/1 contacts and assume all birds are breeding
## created new query called 'metalside' which extracts the recorded side of the body for the metal ring - all birds ringed as chicks should have "L"


head(contacts)  ## CMR data
dim(contacts)

### REMOVE RECORDS FROM BEFORE THE SET START YEAR
contacts<-contacts %>%
  filter(year(Date_Time)>start) 
dim(contacts)


### FIND MISSING DATA FOR SEASON AND REPLACE BASED ON DATE
## changed on 13 July 2019 - use a single calendar year as encounter occasion, not a season (season goes Jan - Nov anyway)

contacts<-contacts %>%
  # mutate(Contact_Season=as.character(Contact_Season)) %>%
  # mutate(MO=month(Date_Time)) %>%
  # mutate(Season1=paste(Contact_Year,(as.numeric(substr(Contact_Year,3,4))+1),sep="-")) %>%
  # mutate(Season2=paste(Contact_Year-1,substr(Contact_Year,3,4),sep="-")) %>%
  # mutate(Contact_Season=if_else(is.na(Contact_Season), if_else(MO>6,Season1,Season2),Contact_Season)) %>%
  mutate(Contact_Year=if_else(is.na(Contact_Year),as.integer(year(Date_Time)),Contact_Year)) %>%
  select(BirdID,Location,Contact_Year)
dim(contacts)
head(contacts)


### CREATE SIMPLE ENCOUNTER HISTORY (0/1)
TRAL_EH<-contacts %>% select(BirdID,Contact_Year) %>%
  mutate(count=1) %>%
  group_by(BirdID,Contact_Year) %>%
  summarise(STATE=max(count)) %>%
  spread(key=Contact_Year, value=STATE, fill=0)
TRAL_EH



#############################################################################
##   6. AGE ASSIGNMENT OF BIRDS FOR SURVIVAL ANALYSIS ###############
#############################################################################

### EXTRACT AGE AT DEPLOYMENT FROM DATABASE
### ACCOUNT FOR THE TIME LAG BETWEEN DEPLOYMENT AND FIRST ENCOUNTER IN ENCOUNTER HISTORY
## BIRDS MARKED AS 'CHICK' PRIOR TO 1994 WILL BE ADULT IN ENCOUNTER HISTORY!

head(ages)
unique(ages$Age)
ages$AGE<-ifelse(ages$Age %in% c("Chick","Fledgling"),0,1)


minage<-ages %>% arrange(BirdID, Date_Time) %>%
  group_by(BirdID) %>%
  summarise(ndepl=length(Age), AGE=min(AGE), Age=first(Age), Sex=first(Sex), Date_Time=first(Date_Time)) %>%
  mutate(AGE = ifelse(year(Date_Time)<1995,1,AGE))      ## BIRDS MARKED AS 'CHICK' PRIOR TO 1995 WILL BE ADULT IN ENCOUNTER HISTORY!
head(minage)

### INSERT AGE INTO ENCOUNTER HISTORY
TRAL_EH$AGE<-minage$AGE[match(TRAL_EH$BirdID,minage$BirdID)]



### FOR BIRDS WITH NO AGE ASSIGNMENT IN DATABASE
### TRY TO ASSIGN AGE AT FIRST MARK FROM SIDE OF BODY OR DATE OF MARKING

head(metalside)
birdage<-metalside %>%   dplyr::filter(SpeciesCode==SP) %>%
  mutate(minage=if_else(Side_Of_Body=="L" & Side_Of_Body_Status=="Known"
                                             ,0,1)) %>%
  mutate(minage=if_else(Contact_Month %in% c(8,9,10,11) & Mark_Status=="Deployed"
                        ,0,1)) %>%
  group_by(SpeciesCode,BirdID) %>%
  summarise(MinAge=min(minage), Date_Time=min(Date_Time)) %>%
  mutate(MinAge = ifelse(year(Date_Time)<1995,1,MinAge))      ## BIRDS MARKED AS 'CHICK' PRIOR TO 1995 WILL BE ADULT IN ENCOUNTER HISTORY!
dim(birdage)


### INSERT AGE INTO ENCOUNTER HISTORY
TRAL_EH$AGE[is.na(TRAL_EH$AGE)]<-birdage$MinAge[match(TRAL_EH$BirdID[is.na(TRAL_EH$AGE)],birdage$BirdID)]
TRAL_EH$AGE[is.na(TRAL_EH$AGE)]<-1 ## fill in the remaining NAs as 'adult'










### CALCULATE PROPORTION OF JUVENILE VS ADULT
table(TRAL_EH$AGE)

TRAL_EH %>% gather(key='Year',value='Contact',-BirdID,-AGE) %>%
  group_by(BirdID,AGE) %>%
  summarise(n_capt=sum(Contact)) %>%
  #group_by(AGE) %>%
  #summarise(n_capt=mean(n_capt))
  
  ggplot() + geom_histogram(aes(x=n_capt)) + facet_wrap(~AGE)





### CALCULATE TIME BETWEEN MARKING AND FIRST RECAPTURE
## clearly shows that there is a problem in assigning birds ringed as JUVENILE!!


head(contacts)

TRAL_EH %>% gather(key='Year',value='Contact',-BirdID,-AGE) %>%
  filter(Contact==1) %>%
  mutate(Year=as.numeric(Year)) %>%
  arrange(BirdID,AGE,Year) %>%
  mutate(gap=dplyr::lead(Year, k=1)-Year) %>%
  filter(!is.na(gap)) %>%
  group_by(BirdID,AGE) %>%
  summarise(PostDeplLag=first(gap)) %>%
  
  ggplot() + geom_histogram(aes(x=PostDeplLag)) + facet_wrap(~AGE)



### IDENTIFY THE PROBLEM BIRDS AND ENCOUNTER HISTORIES
suspects<-TRAL_EH %>% gather(key='Year',value='Contact',-BirdID,-AGE) %>%
  filter(Contact==1) %>%
  mutate(Year=as.numeric(Year)) %>%
  arrange(BirdID,AGE,Year) %>%
  mutate(gap=dplyr::lead(Year, k=1)-Year) %>%
  filter(!is.na(gap)) %>%
  group_by(BirdID,AGE) %>%
  summarise(PostDeplLag=first(gap)) %>%
  filter(AGE==0, PostDeplLag<5) %>%
  left_join(metalside[metalside$Mark_Status=="Deployed",], by="BirdID") %>%
  left_join(minage[,c(1,2,4,5)], by="BirdID") %>%
  select(BirdID,Contact_Year,Date_Time,AGE,Age,Sex,Side_Of_Body,Side_Of_Body_Status,PostDeplLag,ndepl)
  

#fwrite(suspects,"TRAL_suspicious_aged_indiv.csv")






##### FIX AGE ERRORS AND EXPORT ENCOUNTER HISTORY #############################
wrongage<- suspects %>% mutate(month=month(Date_Time)) %>%
  filter(month<8)
TRAL_EH$AGE[TRAL_EH$BirdID %in% wrongage$BirdID]<-1


### EXPORT ENCOUNTER HISTORY
dim(TRAL_EH)
names(TRAL_EH)
#try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival"), silent=T)
fwrite(TRAL_EH[,c(1,22,2:21)],"TRAL_simple_encounter_history_2000_2020.csv")
