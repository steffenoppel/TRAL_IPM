############################################################################
######## DATA PREPARATION FOR INTEGRATED POPULATION MODEL     ##############
############################################################################

### written by Steffen Oppel in December 2018
### steffen.oppel@rspb.org.uk
### uses existing CMR and breeding database to extract data
## updated on 31 March 2020

## completely revised on 5 Jan 2021 to incorporate m-array survival estimation in IPM

library(tidyverse)
library(lubridate)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select


#############################################################################
##   1. SPECIFY THE SPECIES AND START YEAR FOR WHICH YOU WANT A SUMMARY ####
#############################################################################
## SPECIFY THE SPECIES AND START YEAR FOR SURVIVAL MODEL
SP<-"TRAL"
start<-1978  ## for CMR data
IPMstart<-2000 ## for count and breeding success data


###################################################################################
##   2. READ IN DATA FROM DATABASE AND FILTER DATA FOR SPECIES OF INTEREST ####
###################################################################################

## run the RODBC import of CMR data in a 32-bit version of R
system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM\\RODBC_CMR_import_TRAL.R")), wait = TRUE, invisible = FALSE, intern = T)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM"), silent=T)
load("GOUGH_seabird_CMR_data.RData")

## run the RODBC import of nest and count data in a 32-bit version of R
system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database\\RODBC_count_import.r")), wait = TRUE, invisible = FALSE, intern = T)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database"), silent=T)
load("GOUGH_seabird_data.RData")

## filter data for the selected species
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
FECUND<-nests %>% filter(Year<2021) %>% filter(Species==SP) %>% mutate(count=1) %>%
  filter(!NestID %in% exclude$NestID) %>%
  group_by(Year,Colony) %>%
  summarise(n_nests=sum(count),BREED_SUCC=mean(SUCCESS, na.rm=T))

### inspect the dubious summaries
FECUND %>% filter(BREED_SUCC==0)		## Tafelkop 2011 appears to have no monitoring data, data from 2019 incomplete
FECUND %>% filter(BREED_SUCC==1)		## Green Hill North 2007 - only one check in August - worthless

FECUND<-FECUND %>% filter(Colony=="Gonydale")


### PLOT TO SPOT ANY OUTLIERS OF BREEDING SUCCESS
ggplot(FECUND, aes(x=Year,y=BREED_SUCC)) +geom_point(size=2, color='darkred')+geom_smooth(method='lm') 


###   PREPARE A METRIC OF ANNUAL MEAN NEST FAILURE DATE
nestfail<-nests %>% filter(Year<2021) %>% filter(Species==SP) %>% mutate(count=1) %>%
  filter(!NestID %in% exclude$NestID) %>%
  filter(SUCCESS==0) %>%
  mutate(LastDay=yday(DateLastAlive)) %>%
  mutate(LastStage=if_else(is.na(LastStage),ifelse(month(DateLastAlive)>4,"CHIC","INCU"), as.character(LastStage))) %>%
  group_by(Year,LastStage) %>%
  summarise(n=sum(count)) %>% #,LDmean=mean(LastDay, na.rm=T), LDmedian=median(LastDay, na.rm=T)) %>%
  filter(!is.na(LastStage)) %>%
  filter(!(LastStage=="FAIL")) %>%
  spread(key=LastStage, value=n) %>%
  mutate(tot=sum(CHIC,INCU,na.rm=T)) %>%
  mutate(propINCU=INCU/tot)
  



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
  filter(Year>IPMstart) %>%
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
  filter(Year>IPMstart) %>%
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





#############################################################################
##   5. PREPARE THE MARK-RECAPTURE DATA FOR SURVIVAL ANALYSIS ###############
#############################################################################

### COPIED FROM C:\STEFFEN\RSPB\UKOT\Gough\ANALYSIS\SeabirdSurvival\TRAL_survival_marray.r

## filter data for the selected species
contacts<-contacts %>% filter(SpeciesCode==SP) ## %>% filter(Location %in% c("Hummocks","Gonydale","Tafelkop","Not Specified")) - this removes age info for chicks ringed elsewhere!
ages<-ages %>% filter(SpeciesCode==SP)
bands<-bands %>% filter(SpeciesCode==SP)

head(contacts)  ## CMR data
dim(contacts)



#############################################################################
##   6. AGE ASSIGNMENT OF BIRDS FOR SURVIVAL ANALYSIS ###############
#############################################################################

### EXTRACT AGE AT DEPLOYMENT FROM DATABASE
deploy_age<-contacts %>% arrange(BirdID, Date_Time,Contact_Year) %>%
  mutate(AGE=ifelse(Age %in% c("Chick","Fledgling"),0,1)) %>%
  group_by(BirdID) %>%
  summarise(MIN_AGE=min(AGE), MAX_AGE=max(AGE), FIRST_AGE=first(Age), FIRST_Date=first(Date_Time), FIRST_YEAR=min(Contact_Year)) %>%
  mutate(FIRST_AGE=ifelse(FIRST_AGE=="Unknown" & month(FIRST_Date)<5,"Adult", as.character(FIRST_AGE)))  ### unknowns marked before May were not chicks

head(deploy_age)
dim(deploy_age)

MISSAGE<-deploy_age %>% filter(is.na(FIRST_AGE)) %>%   left_join(bands, by="BirdID") %>%
  select(BirdID, Band_Number,MIN_AGE,FIRST_Date,FIRST_YEAR)
dim(MISSAGE)

### FIND MISSING DATA FOR SEASON AND REPLACE BASED ON DATE

contacts<-contacts %>%
  mutate(Contact_Year=if_else(is.na(Contact_Year),as.integer(year(Date_Time)),Contact_Year)) %>%
  mutate(Contact_Year=if_else(as.integer(month(Date_Time))==12 & !(Age %in% c("Chick","Fledgling")),Contact_Year+1,as.numeric(Contact_Year))) %>%
  mutate(Contact_Year=if_else(as.integer(month(Date_Time))==1 & (Age %in% c("Chick","Fledgling")),Contact_Year-1,as.numeric(Contact_Year)))
dim(contacts)
head(contacts)


### ASSIGN AGE TO BIRDS WHERE THIS IS NOT SPECIFIED
## include a column with continuous age 

contacts<-contacts %>%
  left_join(deploy_age, by="BirdID") %>%
  mutate(AGE=ifelse(Age=="Adult",1,ifelse(Age %in% c("Chick","Fledgling"),0,NA))) %>%    ### certain assignments based on provided age
  mutate(AGE=ifelse(is.na(AGE), ifelse(Sex %in% c("Male","Female"),1,NA),AGE)) %>%       ### inferred assignment from sex info - only adults can be sexed
  mutate(ContAge=ifelse(FIRST_AGE %in% c("Chick","Fledgling"),Contact_Year-FIRST_YEAR,Contact_Year-FIRST_YEAR+7)) #%>%      ### continuous age since first deployment, at least 7 years for birds marked as 'adult'

contacts %>% filter(is.na(AGE))
contacts %>% filter(is.na(ContAge))




#############################################################################
##   7. REMOVE BIRDS FROM OUTSIDE THE STUDY AREAS ###############
#############################################################################

########## CREATE A LOOP OVER EVERY BIRD TO CHECK WHETHER THEY WERE EVER RECORDED IN STUDY AREAS
STUDY_AREAS<- c("Hummocks","Gonydale","Tafelkop")
allbirds<-unique(contacts$BirdID)
fixed_contacts<-data.frame()

for (xid in allbirds){
  xcont<-contacts %>% filter(BirdID==xid) %>% arrange(Date_Time)
  xcont$INSIDE<-ifelse(xcont$Location %in% STUDY_AREAS,1,0)
  xcont$INSIDE<-ifelse(xcont$Location == "Not Specified" & xcont$Contact_Year>2014,1,xcont$INSIDE) 
  if(sum(xcont$INSIDE, na.rm=T)>0){fixed_contacts<-bind_rows(fixed_contacts ,xcont)}
}
dim(contacts)
dim(fixed_contacts)
length(unique(fixed_contacts$BirdID))
length(allbirds)

### REMOVE RECORDS FROM BEFORE THE SET START YEAR AND BIRDS FIRST MARKED IN LAST YEAR
contacts<-fixed_contacts %>%
  filter(year(Date_Time)>start) %>%
  filter(ContAge!=1)    ## remove 5 records of unfledged chicks within a few weeks/months of ringing
dim(contacts)
unique(contacts$FIRST_AGE)




#############################################################################
##   8. CREATE MATRIX OF ENCOUNTERS AND AGES ###############
#############################################################################
head(contacts)

### SIMPLE BINARY ENCOUNTER HISTORY FOR CHICKS AND ADULTS
TRAL_CHICK<- contacts %>% mutate(count=1) %>%
  filter(FIRST_AGE %in% c("Chick","Fledgling")) %>%
  group_by(BirdID,Contact_Year) %>%
  summarise(STATE=max(count)) %>%
  spread(key=Contact_Year, value=STATE, fill=0) %>%
  arrange(BirdID)
dim(TRAL_CHICK)

TRAL_AD<- contacts %>% mutate(count=1) %>%
  group_by(BirdID,FIRST_AGE,Contact_Year) %>%
  summarise(STATE=max(count)) %>%
  spread(key=Contact_Year, value=STATE, fill=0) %>%
  filter(FIRST_AGE %in% c("Adult")) %>%    ### filter after spread to ensure that years without any adult contacts (1984, 2003, 2005) are included in matrix
  ungroup() %>%
  select(-FIRST_AGE) %>%
  arrange(BirdID)
dim(TRAL_AD)


### CONVERT TO SIMPLE MATRICES WITHOUT BIRD ID COLUMN
CH.J<-as.matrix(TRAL_CHICK[,2:dim(TRAL_CHICK)[2]])
CH.A<-as.matrix(TRAL_AD[,2:dim(TRAL_CHICK)[2]])


### IDENTIFY WHICH CHICKS WERE EVER RECAPTURED AND PUT THOSE IN A SEPARATE ENCOUNTER HISTORY
cap <- apply(CH.J, 1, sum)
ind <- which(cap >= 2)
CH.J.R <- CH.J[ind,]    # Juvenile CH recaptured at least once
CH.J.N <- CH.J[-ind,]   # Juvenile CH never recaptured


# FOR THOSE CHICKS THAT WERE RECAPTURED, Remove first capture and add the rest to the adult encounter history
first <- numeric()
for (i in 1:dim(CH.J.R)[1]){
  first[i] <- min(which(CH.J.R[i,]==1))
}
CH.J.R1 <- CH.J.R
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R1[i,first[i]] <- 0
}

# Add grown-up juveniles to adults FOR COMPLETE ENCOUNTER HISTORY TO BUILD ADULT MARRAY
CH.A.m <- rbind(CH.A, CH.J.R1)


# Create ENCOUNTER HISTORY matrix for juveniles, ignoring subsequent recaptures
second <- numeric()
for (i in 1:dim(CH.J.R1)[1]){
  second[i] <- min(which(CH.J.R1[i,]==1))
}
CH.J.R2 <- matrix(0, nrow = dim(CH.J.R)[1], ncol = dim(CH.J.R)[2])
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R2[i,first[i]] <- 1
  CH.J.R2[i,second[i]] <- 1
}



#############################################################################
##   9. CONVERT TO MARRAY ###############
#############################################################################

# Function to create a m-array based on capture-histories (CH)
marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}


# Create TWO MARRAYS from the capture-histories: one for chicks and one for adults (including those chicks that were ever recaptured)

### ADULT MARRAY
adult.marray <- marray(CH.A.m)

# Create m-array for the chicks that were ultimately recaptured, but only up to the first recapture (the rest is included in the adult.marray)
CH.J.R.marray <- marray(CH.J.R2)
# The last column ought to show the number of juveniles not recaptured again and should all be zeros, since all of them are released as adults
CH.J.R.marray[,dim(CH.J)[2]] <- 0

# Create the m-array for chicks never recaptured and add it to the recaptured chicks m-array
CH.J.N.marray <- marray(CH.J.N)
chick.marray <- CH.J.R.marray + CH.J.N.marray 



### CALCULATE THE PROPORTION OF CHICKS RECRUITING AT A CERTAIN AGE

recruit.age<-function(CH){
  n.occasions <- dim(CH)[2]-1
  total<-as.numeric()
  recruit.mat<- matrix(data = 0, ncol = n.occasions, nrow = n.occasions-1)
  # Calculate the age proportion of returned individuals at each time period
  for (t in 1:(n.occasions-1)){
    total[t] <- sum(CH[t,])
    for (col in (t+1):n.occasions){
      recruit.mat[t,col-t]<- CH[t,col]
    }
  }
  return(list(REC=recruit.mat,TOT=total))
}

RECRUIT.AGE.MAT<-recruit.age(CH.J.R.marray)
RECRUIT.AGE<-data.frame(age=seq(1,dim(RECRUIT.AGE.MAT$REC)[2],1),N=apply(RECRUIT.AGE.MAT$REC,2,sum)) %>%
  mutate(prop=N/sum(RECRUIT.AGE.MAT$TOT))

sum(RECRUIT.AGE.MAT$TOT)
sum(RECRUIT.AGE$N)

ggplot(RECRUIT.AGE) + geom_bar(aes(x=age,y=prop),stat='identity', fill='cornflowerblue') + 
  ylab("Proportion of TRAL first seen on Gough") + 
  xlab("Age in years") + 
  theme(panel.background=element_rect(fill="white", colour="black"),  
        axis.text=element_text(size=16, color="black"),  
        axis.title=element_text(size=18),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.border = element_blank()) 



#################################################################################################################
##   10. ATTEMPT TO QUANTIFY THE LAST NEST ALIVE DATE FOR BIRDS TO STILL RETURN IN FOLLOWING YEAR ###############
#################################################################################################################

### EXTRACT UNSUCCESSFUL NESTS
failures<- nests %>% filter(SUCCESS==0) %>% filter(Species==SP) %>%
  mutate(DateLastAlive=if_else(year(DateLastAlive)<Year,DateLastChecked,DateLastAlive)) ## avoid days >300 for nests that failed in December

### SELECT INDIVIDUALS THAT WERE RECORDED AT THESE UNSUCCESSFUL NESTS
fail.ind<-contacts %>% filter(Nest_Description %in% unique(failures$Nest_label))

### EXTRACT ALL CONTACTS OF INDIVIDUALS ASSOCIATED WITH ANY FAILED NEST
nestfail.contacts<-contacts %>% filter(BirdID %in% unique(fail.ind$BirdID)) %>% group_by(BirdID, Contact_Year) %>%
  summarise(SEEN=n()) %>%
  arrange(BirdID, Contact_Year) %>%
  rename(Year=Contact_Year) 

### EXTRACT ALL INDIVIDUALS ASSOCIATED WITH ALL FAILED NESTS
fail.nest.ind<-fail.ind %>% rename(Nest_label=Nest_Description, Year=Contact_Year) %>%
  left_join(failures, by=c("Year","Nest_label")) %>%
  group_by(Nest_label,Year,LastStage,DateLastAlive,DateLastChecked, BirdID) %>%
  summarise(nest.visits=n())


### COMBINE THE ABOVE IN A FULL MATRIX AND EXTRACT RECORDS BASED ON PREVIOUS YEAR NEST FAILURE

TRAL_NESTFAIL_MATRIX<- expand.grid(BirdID=unique(fail.ind$BirdID), Year=seq(2007,2021,1)) %>%
  left_join(fail.nest.ind, by=c("Year","BirdID")) %>%
  left_join(nestfail.contacts, by=c("Year","BirdID")) %>%
  mutate(SEEN=ifelse(is.na(SEEN),0,SEEN)) %>%
  mutate(SEEN=ifelse(SEEN>1,1,SEEN)) %>%
  arrange(BirdID,Year) %>%
  mutate(prevNest=dplyr::lag(Nest_label),prevStage=dplyr::lag(LastStage),prevAlive=yday(dplyr::lag(DateLastAlive)),prevCheck=yday(dplyr::lag(DateLastChecked))) %>%
  filter(!is.na(prevNest)) %>%
  arrange(BirdID,Year) %>%
  select(BirdID,Year,SEEN,Nest_label,LastStage,DateLastAlive,prevNest,prevStage,prevAlive,prevCheck)
dim(TRAL_NESTFAIL_MATRIX)


### CALCULATE RETURN PROBABILITY BASED ON DATE LAST SEEN ALIVE
ret.mod<-glm(SEEN~prevAlive,TRAL_NESTFAIL_MATRIX, family=binomial)
TRAL_NESTFAIL_MATRIX$pred<-predict(ret.mod,newdat=TRAL_NESTFAIL_MATRIX, type="response")
ggplot(TRAL_NESTFAIL_MATRIX) + geom_point(aes(x=prevAlive,y=SEEN)) + 
  geom_line(aes(x=prevAlive,y=pred), colour='darkred', size=1.5) +
  ylab("Returning to Gough") + 
  xlab("Day when last alive in previous year") + 
  theme(panel.background=element_rect(fill="white", colour="black"),  
        axis.text.y=element_text(size=18, color="black"), 
        axis.text.x=element_text(size=14, color="black", angle=45, vjust=0.5),  
        axis.title=element_text(size=20),  
        legend.position=c(0.15,0.9), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.border = element_blank()) 


### TROUBLESHOOT QUESTIONABLE RECORDS
TRAL_NESTFAIL_MATRIX %>% filter(SEEN==1) %>% filter(prevAlive>250)


### EXTRACT THE DATE WHEN RETURN PROBABILITY DROPS BELOW 50%
CUTOFF <- approx(x = ret.mod$fitted.values, y = TRAL_NESTFAIL_MATRIX$prevAlive[!is.na(TRAL_NESTFAIL_MATRIX$prevAlive)], xout=0.5)$y


### CALCULATE PROPORTION OF NESTS THAT FAIL AFTER THAT DATE
head(nests)

FAIL.PROP<-nests %>% #filter(SUCCESS==0) %>% #filter(Year==2011) %>%
  mutate(DateLastAlive=if_else(SUCCESS==1,ymd_hms(paste(Year,"12","01 12:00:00",sep="-")),DateLastAlive)) %>%
  mutate(DateLastAlive=if_else(year(DateLastAlive)<Year,DateLastChecked,DateLastAlive)) %>% ## avoid days >300 for nests that failed in December
  mutate(DateLastAlive=if_else(is.na(DateLastAlive),DateLastChecked,DateLastAlive)) %>% ## avoid missing data for 2006-2009 and 2011
  mutate(LATEFAIL=if_else(yday(DateLastAlive)>CUTOFF,1,0)) %>%
  group_by(Year) %>%
  summarise(prop.late=mean(LATEFAIL, na.rm=T), mean.last.alive=mean(yday(DateLastAlive), na.rm=T), median.last.alive=median(yday(DateLastAlive), na.rm=T)) %>%
  mutate(prop.late=ifelse(prop.late==1,NA,prop.late)) %>%
  bind_rows(data.frame(Year=2005, prop.late=NA)) %>%
  arrange(Year)



### SIMPLE EXPLORATION WHETHER FAIL PROP EXPLAINS SEQUENCE OF COUNT DATA
### PREPARE COUNT DATA AS IN IPM POPULATION TREND ######
TRAL.pop<-POPSIZE
TRAL.pop[14,4]<-136   ### number of nests monitored in Gonydale that year
TRAL.pop<-TRAL.pop %>% gather(key='Site', value='Count',-Year) %>%
  filter(Year>2003) %>%
  mutate(Site=ifelse(Site %in% c('GP Valley','West Point'),'GP Valley',Site)) %>%
  mutate(Site=ifelse(Site %in% c('Gonydale','Green Hill','Hummocks'),'Gonydale',Site)) %>%
  group_by(Year,Site) %>%
  summarise(Count=sum(Count, na.rm=T)) %>%
  mutate(Count=ifelse(Count==0,NA,Count)) %>%
  spread(key=Site, value=Count)
TRAL.props<-prop.table(as.matrix(TRAL.pop[,2:9]),1)
mean.props<-apply(TRAL.props[c(1,3:7,9:13,15:17),],2,mean) ## for start in 2004
TRAL.pop$prop.counted<-0
for (l in 1:length(TRAL.pop$Year)){
  TRAL.pop$prop.counted[l]<-sum(mean.props[which(!is.na(TRAL.pop[l,2:9]))])
}
TRAL.pop$tot<-rowSums(TRAL.pop[,2:9], na.rm=T)

### CALCULATE CHANE FROM ONE YEAR TO NEXT
TRAL.pop %>% ungroup() %>%
  mutate(tot=ifelse(Year==2011,NA,tot)) %>%
  mutate(change=((tot/prop.counted)-(dplyr::lag(tot)/dplyr::lag(prop.counted)))/(dplyr::lag(tot)/dplyr::lag(prop.counted))) %>%
  select(Year, tot,prop.counted,change) %>%
  left_join(FAIL.PROP, by="Year") %>%
  gather(key="Prev.breed.measure",value="Metric",-Year,-tot,-change,-prop.counted) %>%

  ggplot() + geom_point(aes(x=Metric,y=change), size=2, col='darkred') + 
  geom_smooth(aes(x=Metric,y=change), col="darkblue",method='lm') +
    facet_wrap(~Prev.breed.measure, ncol=1, scales="free_x") +
  ylab("Proportional change in N pairs") + 
  xlab("Measure of previous breeding season end") + 
  theme(panel.background=element_rect(fill="white", colour="black"),  
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),  
        strip.text.x=element_text(size=18, color="black"),  
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.border = element_blank()) 


#############################################################################
##   11. SAVE WORKSPACE ###############
#############################################################################
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
save.image("TRAL_IPM_input.marray.RData")
