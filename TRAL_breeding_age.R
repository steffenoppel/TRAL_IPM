#####################################################################################
######## EXAMINATION OF BREEDING AGE OF TRISTAN ALBATROSS ON GOUGH     ##############
#####################################################################################

### written by Steffen Oppel in August 2021
### steffen.oppel@rspb.org.uk
### examines 3 question re TRAL breeding age:
## 1. does average age of all TRAL first recruiters in a given year change over time? [This is potentially effort-dependent]
## 2. does average age of all TRAL breeders in a given year change over time?
## 3. does average age of all TRAL individuals contacted in a given year change over time? - not used in manuscript as age increases due to length of time since ringing was started

## update on 2 Oct 2021: after chat with Cat Horswill try to examine whether mean age at first breeding 2004-2009 is higher than 2015-2021
## switched questions and figures for manuscript, because prop old breeders is more logical and should be key evidence

## updated 18 February 2022: repeated questions that increase in age proportion could simply be due to when ringing began - examine in more detail

library(tidyverse)
library(lubridate)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select

# LOAD AND MANIPULATE ICONS
library(grid)
library(magick)
imgTRAL<-image_read("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\PR_Comms\\Icons\\alby 4.jpg") %>% image_transparent("white", fuzz=5)
TRALicon <- rasterGrob(imgTRAL, interpolate=TRUE)


#############################################################################
##   DATA IMPORT AND PREPARATION COPIED FROM IPM_DATA_PREPARATION.r     ####
#############################################################################
## SPECIFY THE SPECIES AND START YEAR FOR SURVIVAL MODEL
SP<-"TRAL"
start<-1950  ## for CMR data
IPMstart<-2004 ## for count and breeding success data

## run the RODBC import of CMR data in a 32-bit version of R
#system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM\\RODBC_CMR_import_TRAL.R")), wait = TRUE, invisible = FALSE, intern = T)
#system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM\\RODBC_CMR_import_TRAL.R")), wait = TRUE, invisible = FALSE, intern = T)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM"), silent=T)
load("GOUGH_seabird_CMR_data.RData")

## run the RODBC import of nest and count data in a 32-bit version of R
#system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database\\RODBC_count_import.r")), wait = TRUE, invisible = FALSE, intern = T)
#system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database\\RODBC_count_import.r")), wait = TRUE, invisible = FALSE, intern = T)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Breeding_Database"), silent=T)
load("GOUGH_seabird_data.RData")

## filter data for the selected species
nests<-nests %>% filter(Species==SP)

## filter data for the selected species
contacts<-contacts %>% filter(SpeciesCode==SP) ## %>% filter(Location %in% c("Hummocks","Gonydale","Tafelkop","Not Specified")) - this removes age info for chicks ringed elsewhere!
ages<-ages %>% filter(SpeciesCode==SP)
bands<-bands %>% filter(SpeciesCode==SP)

head(contacts)  ## CMR data
dim(contacts)
dim(contactsbreed)


### EXTRACT AGE AT DEPLOYMENT FROM DATABASE
deploy_age<-contacts %>% arrange(BirdID, Date_Time,Contact_Year) %>%
  mutate(AGE=ifelse(Age %in% c("Chick","Fledgling"),0,1)) %>%
  group_by(BirdID) %>%
  summarise(MIN_AGE=min(AGE), MAX_AGE=max(AGE), FIRST_AGE=first(Age), FIRST_Date=first(Date_Time), FIRST_YEAR=min(Contact_Year)) %>%
  mutate(FIRST_AGE=ifelse(FIRST_AGE=="Unknown" & month(FIRST_Date)<5,"Adult", as.character(FIRST_AGE))) %>%  ### unknowns marked before May were not chicks
  mutate(FIRST_AGE=ifelse(is.na(FIRST_AGE), ifelse(month(FIRST_Date)<5,"Adult","Chick"), as.character(FIRST_AGE)))  ### unknowns marked before May were not chicks

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


### CHECK WHAT BIRDS WERE RINGED AS CHICKS BEFORE 1978
oldchicks<-fixed_contacts %>% filter(Contact_Year<=start) %>% filter(ContAge<2)
fixed_contacts %>% filter(BirdID %in% oldchicks$BirdID)


### REMOVE RECORDS FROM BEFORE THE SET START YEAR AND BIRDS FIRST MARKED IN LAST YEAR
contacts<-fixed_contacts %>%
  filter(year(Date_Time)>start) %>%
  filter(ContAge!=1)    ## remove 5 records of unfledged chicks within a few weeks/months of ringing
dim(contacts)
unique(contacts$FIRST_AGE)


## try to determine years with high and low detection probability

contacts %>% mutate(count=1) %>% group_by(Contact_Year) %>% summarise(n=sum(count)) %>%
  ggplot() + geom_bar(aes(x=Contact_Year,y=n), stat="identity")

n_exist<-deploy_age %>% mutate(count=1) %>% rename(Contact_Year=FIRST_YEAR) %>%
  mutate(FIRST_AGE=if_else(FIRST_AGE=="Fledgling","Chick",FIRST_AGE)) %>%
  group_by(Contact_Year,FIRST_AGE) %>%
  summarise(N_marked=sum(count)) %>%
  arrange(Contact_Year) %>%
  spread(key=FIRST_AGE, value=N_marked, fill=0) %>%
  ungroup() %>%
  mutate(N_marked=Adult+Chick) %>%
  mutate(N_ever_ad = cumsum(Adult),N_ever_ch = cumsum(Chick),N_all = cumsum(N_marked)) %>%
  bind_rows(tibble(Contact_Year=2022,Adult=0,Chick=0,N_marked=0,N_all=0, N_ever_ad=1404, N_ever_ch=3635)) %>%
  mutate(N_all=if_else(Contact_Year==2022,dplyr::lag(N_all),N_all)) %>%
  arrange(Contact_Year)  
n_exist$N_all[1] = n_exist$N_marked[1]
n_exist$N_all[2] = (n_exist$Adult[1]*0.96) + (n_exist$Chick[1]*0.85) + n_exist$N_marked[2]
n_exist$N_all[3] = (n_exist$N_all[2]*(0.96^19)) + n_exist$N_marked[3]
n_exist$N_all[4] = (n_exist$Adult[3]*0.96) + (n_exist$Chick[3]*0.85) + (n_exist$N_all[2]*(0.96^19)) + n_exist$N_marked[4]
for (y in 5:dim(n_exist)[1]) {
  n_exist$N_all[y] = ((n_exist$Adult[y-1]+n_exist$N_all[y-2])*0.96) + (n_exist$Chick[y-1]*0.85) + n_exist$N_marked[y]
}
tail(n_exist)


### ADDED - A COLUMN THAT SAYS WHAT PROPORTION OF BIRDS WERE MARKED > 30 YEARS AGO

goodyears<-contacts %>% group_by(Contact_Year) %>% summarise(n=length(unique(BirdID))) %>%
  filter(!(Contact_Year==1959)) %>%   ## remove single re-sighting from 1959
  left_join(n_exist, by='Contact_Year') %>%
  mutate(prop.seen=n/N_all) %>%
  mutate(old.ad=ifelse(Contact_Year<1981,0,
                             ifelse(Contact_Year<2002,
                                    142,
                                    dplyr::lag(N_ever_ad,n=26)))) %>% ## specifies adults that were already 4 years old when ringed

  mutate(old.ch=ifelse(Contact_Year<1985,0,
                             ifelse(Contact_Year<2006,
                                    21,
                                    dplyr::lag(N_ever_ch,n=30)))) %>%  ## specifies chicks that were already 4 years old when ringed
  mutate(all.pot.old=(old.ad+old.ch)) %>%
  mutate(all.pot.ad=cumsum(Adult), all.pot.ch=cumsum(Chick)) %>%
  mutate(all.pot.breed=dplyr::lag(all.pot.ch,n=4)+all.pot.ad) %>%
  mutate(prop.pot.old=all.pot.old/(all.pot.breed-all.pot.old)) %>%   ## changed denominator to all.not.old on advice from Adam Butler
  select(Contact_Year,n,Adult,Chick,N_marked,N_all,prop.seen,prop.pot.old)

#filter(Contact_Year>1978)
tail(goodyears)





#############################################################################
##   QUESTION 1: does age of breeders change over time? ###############
#############################################################################

### THIS QUESTION IS MASSIVELY AFFECTED BY WHEN THE RINGING BEGAN AS THERE IS AN INCREASING AGE OVER TIME

### this requires the breeding status ID to be added to the access query

breeders<-contacts %>%
  filter(!(FIRST_YEAR==Contact_Year)) %>%   ## potentially change this to remove only ringed chicks? Age %in% c("Chick","Fledgling")
  #filter(!(is.na(Nest_Description))) %>%
  mutate(Breeding_StatusID=ifelse(is.na(Nest_Description),Breeding_StatusID,1)) %>%
  filter(Breeding_StatusID %in% c(1,-1525788936,105568723,1899636611,1899636612,1899636618)) %>%
  group_by(BirdID,Contact_Year) %>%
  summarise(ContAge=mean(ContAge)) %>%
  left_join(goodyears, by="Contact_Year") %>%
  mutate(detrend=Contact_Year-1980) %>%
  filter(Contact_Year>1990) 
dim(breeders)
min(breeders$ContAge)
breeders %>% filter(ContAge<6)

### exploratory plots
ggplot(breeders) +
  geom_point(aes(x=Contact_Year, y=ContAge)) +
  geom_smooth(aes(x=Contact_Year, y=ContAge),method="lm")

# ggplot(breeders) +
#   geom_histogram(aes(x=ContAge)) +
#   facet_wrap(~Contact_Year)

### analysis
m2eff<-glm(ContAge~Contact_Year, data=breeders, family="poisson",weights=prop.seen)
summary(m2eff)

predict(m2eff, newdat=data.frame(Contact_Year=seq(2000,2020,1),detrend=1), type='response')

### does the proportion of young breeders change over time?
## removed on 2 Oct 2021 because it is too confusing and will lead to questions of age at first breeding which we cannot answer

### does the proportion of OLD breeders change over time?

## first, calculate the proportion of birds in each year that could have plausibly been >30 years
head(goodyears)

oldbreeders<-contacts %>%
  filter(!(FIRST_YEAR==Contact_Year)) %>%   ## potentially change this to remove only ringed chicks? Age %in% c("Chick","Fledgling")
  mutate(Breeding_StatusID=ifelse(is.na(Nest_Description),Breeding_StatusID,1)) %>%
  filter(Breeding_StatusID %in% c(1,-1525788936,105568723,1899636611,1899636612,1899636618)) %>%
  group_by(BirdID,Contact_Year) %>%
  summarise(ContAge=mean(ContAge)) %>%
  left_join(goodyears, by="Contact_Year") %>%
  ungroup() %>%
  mutate(YOUNG=ifelse(ContAge<30,1,0)) %>%
  mutate(OLD=ifelse(ContAge>29,1,0)) %>%
  group_by(Contact_Year) %>%
  summarise(prop.old=mean(OLD),n.young=sum(YOUNG),n.old=sum(OLD)) %>%
  left_join(goodyears, by="Contact_Year") %>%
  filter(Contact_Year>2003) 
dim(oldbreeders)
#fwrite(oldbreeders,"TRAL_old_breeders_2004_2021.csv")

### analysis of trend over time

m2oleff<-glm(cbind(n.old,n.young)~Contact_Year+offset(log(prop.pot.old)), data=oldbreeders, family=binomial(link="cloglog"),weights=prop.seen)
summary(m2oleff)
str(m2oleff)
m2oleff$fitted.values
### prediction of effect size
olddat<-data.frame(Contact_Year=seq(2004,2021,1),prop.seen=1, prop.pot.old=0.01)

## grad the inverse link function
ilink <- family(m2oleff)$linkinv
## add fit and se.fit on the **link** scale
olddat <- bind_cols(olddat, setNames(as_tibble(predict(m2oleff, olddat, se.fit = TRUE)[1:2]),
                                     c('fit_link','se_link')))
## create the interval and backtransform
olddat <- mutate(olddat,
                 pred.prop  = ilink(fit_link),
                 ucl = ilink(fit_link + (1.96 * se_link)),
                 lcl = ilink(fit_link - (1.96 * se_link)))



### COMBINE PROPORTION OF YOUNG AND OLD BREEDERS IN ONE PLOT
ggplot(olddat) +

  geom_line(data=oldbreeders,aes(x=Contact_Year, y=prop.pot.old),colour = "darkgrey",linetype="dashed") +
  geom_line(aes(x=Contact_Year, y=pred.prop),colour = "indianred") +
  geom_ribbon(aes(x=Contact_Year, ymin=lcl, ymax=ucl), fill= "indianred",alpha = 0.2) +
  geom_point(data=oldbreeders,aes(x=Contact_Year, y=prop.old), size=3) +
  labs(x = "Year",
       y = "Annual proportion of old breeders") +
  scale_y_continuous(breaks=seq(0,0.40,0.10), limits=c(0,0.40))+
  scale_x_continuous(breaks=seq(2005,2021,2), limits=c(2004,2021))+
  
  
  ### add the bird icons
  annotation_custom(TRALicon, xmin=2004.5, xmax=2006.5, ymin=0.06, ymax=0.10) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        panel.grid.minor = element_blank())

ggsave("C:\\STEFFEN\\MANUSCRIPTS\\Submitted\\TRAL_IPM\\Fig3_rev.jpg", width=9, height=6)




################### CHECK WHETHER THE NUMBER OF OLD BIRDS CHANGES #####################################################
## ####

### does the proportion of OLD breeders change over time?

oldbirds<-contacts %>%
  filter(!(FIRST_YEAR==Contact_Year)) %>%   ## potentially change this to remove only ringed chicks? Age %in% c("Chick","Fledgling")
  group_by(BirdID,Contact_Year) %>%
  summarise(ContAge=mean(ContAge)) %>%
  left_join(goodyears, by="Contact_Year") %>%
  ungroup() %>%
  mutate(YOUNG=ifelse(ContAge<30,1,0)) %>%
  mutate(OLD=ifelse(ContAge>29,1,0)) %>%
  group_by(Contact_Year) %>%
  summarise(prop.old=mean(OLD),n.young=sum(YOUNG),n.old=sum(OLD)) %>%
  left_join(goodyears, by="Contact_Year") %>%
  filter(Contact_Year>2003) 


### analysis of trend over time

m2oleff<-glm(cbind(n.old,n.young)~Contact_Year, data=oldbirds, family="binomial", weights=prop.seen)
summary(m2oleff)

### prediction of effect size
olddat<-data.frame(Contact_Year=seq(2004,2021,1),prop.seen=1)

## grad the inverse link function
ilink <- family(m2oleff)$linkinv
## add fit and se.fit on the **link** scale
olddat <- bind_cols(olddat, setNames(as_tibble(predict(m2oleff, olddat, se.fit = TRUE)[1:2]),
                                     c('fit_link','se_link')))
## create the interval and backtransform
olddat <- mutate(olddat,
                 pred.prop  = ilink(fit_link),
                 ucl = ilink(fit_link + (1.96 * se_link)),
                 lcl = ilink(fit_link - (1.96 * se_link)))








#############################################################################
##   QUESTION 2: does age of first returners change over time? ###############
#############################################################################
firstreturnyear<-contacts %>%
  filter(FIRST_AGE %in% c("Chick","Fledgling")) %>%
  filter(!(FIRST_YEAR==Contact_Year)) %>%
  group_by(BirdID) %>%
  summarise(FirstReturn=min(Contact_Year))

firstreturns<-contacts %>%
  filter(FIRST_AGE %in% c("Chick","Fledgling")) %>%
  left_join(goodyears, by="Contact_Year") %>%
  left_join(firstreturnyear, by="BirdID") %>%
  filter((Contact_Year==FirstReturn)) %>%
  select(ContactID,BirdID,Contact_Year,FIRST_YEAR,ContAge,n,N_marked,prop.seen,FirstReturn) %>%
  rename(effort=n) %>%
  filter(Contact_Year>2003) %>%
  filter(Contact_Year!=2008) %>%  ### remove the THREE recruiters observed in 2008
  filter(Contact_Year!=2005)   ### remove the ONLY recruiter observed in 2005!
dim(firstreturns)


### exploratory plots
ggplot(firstreturns) +
  geom_point(aes(x=Contact_Year, y=ContAge), position=position_jitter()) +
  geom_smooth(aes(x=Contact_Year, y=ContAge),method="lm")

ggplot(firstreturns) +
  geom_histogram(aes(x=ContAge))

### analysis
m3eff<-glm(ContAge~Contact_Year, data=firstreturns, family="poisson",weights=prop.seen)
summary(m3eff)
m3n<-glm(ContAge~Contact_Year, data=firstreturns, family="poisson",weights=effort)
summary(m3n)
m3<-glm(ContAge~Contact_Year, data=firstreturns, family="poisson")
summary(m3)


### prediction of effect size
newdat<-data.frame(Contact_Year=seq(2004,2021,1),prop.seen=1)
#newdat$pred.age<-predict(m3eff, newdat=newdat, type="response", se=T)$fit
#newdat$se.age<-predict(m3eff, newdat=newdat, type="response", se=T)$se.fit
#newdat<-newdat %>% mutate(lcl=pred.age-1.96 * se.age,ucl=pred.age+1.96 * se.age)

## grad the inverse link function
ilink <- family(m3eff)$linkinv
## add fit and se.fit on the **link** scale
newdat <- bind_cols(newdat, setNames(as_tibble(predict(m3eff, newdat, se.fit = TRUE)[1:2]),
                                     c('fit_link','se_link')))
## create the interval and backtransform
newdat <- mutate(newdat,
                 pred.age  = ilink(fit_link),
                 ucl = ilink(fit_link + (1.96 * se_link)),
                 lcl = ilink(fit_link - (1.96 * se_link)))






### plot predicted effect size
ggplot(newdat) +
  geom_line(aes(x=Contact_Year, y=pred.age),colour = "blue") +
  geom_ribbon(aes(x=Contact_Year, ymin=lcl, ymax=ucl), alpha = 0.2,fill = "blue") +
  geom_point(data=firstreturns,aes(x=Contact_Year, y=ContAge), colour="grey45",size=0.7,position=position_jitter()) +
  
  ylab("Age (years) of first return to Gough") +
  xlab("Year") +
  scale_y_continuous(breaks=seq(0,30,5), limits=c(0,32))+
  scale_x_continuous(breaks=seq(2005,2021,2), limits=c(2004,2021))+
  
  ### add the bird icons
  annotation_custom(TRALicon, xmin=2004, xmax=2006, ymin=25, ymax=30) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        panel.grid.minor = element_blank())

ggsave("C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\TRAL_IPM\\Fig3.jpg", width=9, height=6)




######### ANALYSIS WITH QUANTILE REGRESSION ######
## makes no big difference, hence discarded
#install.packages("quantreg")
# library(quantreg)
# 
# firstreturns %>% group_by(Contact_Year) %>% summarise(med=median(ContAge),mean=mean(ContAge))
# firstreturns %>% filter(Contact_Year==2006)
# 
# ### analysis
# m3eff_rq<-rq(ContAge~Contact_Year, tau=c(0.025, 0.5, 0.975),data=firstreturns, weights=prop.seen)
# summary(m3eff_rq)
# plot(m3eff_rq,mfrow = c(1,2))
# 
# 
# plot(firstreturns$Contact_Year,firstreturns$ContAge,cex=.25,type="n",xlab="Year", ylab="Age")
# points(firstreturns$Contact_Year,firstreturns$ContAge,cex=.5,col="blue")
# abline(rq(ContAge~Contact_Year, data=firstreturns, weights=prop.seen),col="blue")
# abline(lm(ContAge~Contact_Year, data=firstreturns, weights=prop.seen),lty=2,col="red") #the mean regression line
# taus <- c(0.025, 0.975)
# for(i in 1:length(taus)){
#   abline(rq(ContAge~Contact_Year, tau=taus[i],data=firstreturns, weights=prop.seen),col="gray")
#   }







#############################################################################
##   QUESTION 3: does age of first breeding change over time? ###############
#############################################################################
unique(contacts$Breeding_StatusID)
firstbreedyear<-contacts %>%
  filter(FIRST_AGE %in% c("Chick","Fledgling")) %>%
  filter(!(FIRST_YEAR==Contact_Year)) %>%
  filter(Breeding_StatusID %in% c(1,1899636611,1899636612,1899636615,1899636613,105568723,1899636616,-1525788936)) %>%
  filter(!(is.na(Nest_Description))) %>%
  group_by(BirdID,FIRST_YEAR) %>%
  summarise(FirstBreed=min(Contact_Year),AgeFirstBreed=min(ContAge))
dim(firstbreedyear)

### exploratory plots
ggplot(firstbreedyear) +
  geom_point(aes(x=FirstBreed, y=AgeFirstBreed), position=position_jitter()) +
  geom_smooth(aes(x=FirstBreed, y=AgeFirstBreed),method="lm")


### decadal analysis
firstbreedyear %>%
  mutate(decade=if_else(FirstBreed<2000,"1990s",if_else(FirstBreed<2015,"2000-2015","2015-2021"))) %>%
  group_by(decade) %>%
  summarise(med=median(AgeFirstBreed),lcl=quantile(AgeFirstBreed,0.05),ucl=quantile(AgeFirstBreed,0.95))







#################################################################################################################
##   QUESTION 4: does age of all recorded bird change over time? ABANDONED BECAUSE NOT MEANINGFUL ###############
#################################################################################################################
head(contacts)
dim(contacts)
returns<-contacts %>%
  filter(!(FIRST_YEAR==Contact_Year)) %>%   ## potentially change this to remove only ringed chicks? Age %in% c("Chick","Fledgling")
  group_by(BirdID,Contact_Year) %>%
  summarise(ContAge=mean(ContAge)) %>%
  left_join(goodyears, by="Contact_Year") %>%
  filter(Contact_Year>2003) 
dim(returns)


### exploratory plots
ggplot(returns) +
  geom_point(aes(x=Contact_Year, y=ContAge)) +
  geom_smooth(aes(x=Contact_Year, y=ContAge),method="lm")

ggplot(returns) +
  geom_histogram(aes(x=ContAge))

### analysis
m1eff<-glm(ContAge~Contact_Year, data=returns, family="poisson",weights=prop.seen)
#m1n<-glm(ContAge~Contact_Year, data=returns, family="poisson",weights=n)
#m1<-glm(ContAge~Contact_Year, data=returns, family="poisson")
summary(m1eff)
#summary(m1n)
#summary(m1)



