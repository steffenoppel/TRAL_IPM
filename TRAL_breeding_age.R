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
start<-1978  ## for CMR data
IPMstart<-2000 ## for count and breeding success data

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

## calculate number of individuals that had been marked by a given year
n_exist<-deploy_age %>% mutate(count=1) %>% rename(Contact_Year=FIRST_YEAR) %>%
  group_by(Contact_Year) %>%
  summarise(N_marked=sum(count)) %>%
  arrange(Contact_Year) %>%
  mutate(N_all = cumsum(N_marked))

goodyears<-contacts %>% mutate(count=1) %>% group_by(Contact_Year) %>% summarise(n=sum(count)) %>%
  left_join(n_exist, by='Contact_Year') %>%
  mutate(prop.seen=n/N_all) %>%
  mutate(p.sel=if_else(prop.seen>0.1,2,1))
tail(goodyears)





#############################################################################
##   QUESTION 1: does age of breeders change over time? ###############
#############################################################################

### this requires the breeding status ID to be added to the access query

breeders<-contacts %>%
  filter(!(FIRST_YEAR==Contact_Year)) %>%   ## potentially change this to remove only ringed chicks? Age %in% c("Chick","Fledgling")
  #filter(!(is.na(Nest_Description))) %>%
  mutate(Breeding_StatusID=ifelse(is.na(Nest_Description),Breeding_StatusID,1)) %>%
  filter(Breeding_StatusID %in% c(1,-1525788936,105568723,1899636611,1899636612,1899636618)) %>%
  group_by(BirdID,Contact_Year) %>%
  summarise(ContAge=mean(ContAge)) %>%
  left_join(goodyears, by="Contact_Year") %>%
  filter(Contact_Year>2009) 
dim(breeders)


### exploratory plots
ggplot(breeders) +
  geom_point(aes(x=Contact_Year, y=ContAge)) +
  geom_smooth(aes(x=Contact_Year, y=ContAge),method="lm")

ggplot(breeders) +
  geom_histogram(aes(x=ContAge)) +
  facet_wrap(~Contact_Year)

### analysis
m2eff<-glm(ContAge~Contact_Year, data=breeders, family="poisson",weights=prop.seen)
summary(m2eff)



### does the proportion of young breeders change over time?
## removed on 2 Oct 2021 because it is too confusing and will lead to questions of age at first breeding which we cannot answer
# 
# youngbreeders<-contacts %>%
#   filter(!(FIRST_YEAR==Contact_Year)) %>%   ## potentially change this to remove only ringed chicks? Age %in% c("Chick","Fledgling")
#   mutate(Breeding_StatusID=ifelse(is.na(Nest_Description),Breeding_StatusID,1)) %>%
#   filter(Breeding_StatusID %in% c(1,-1525788936,105568723,1899636611,1899636612,1899636618)) %>%
#   group_by(BirdID,Contact_Year) %>%
#   summarise(ContAge=mean(ContAge)) %>%
#   left_join(goodyears, by="Contact_Year") %>%
#   ungroup() %>%
#   mutate(YOUNG=ifelse(ContAge<10,1,0)) %>%
#   mutate(OLD=ifelse(ContAge>9,1,0)) %>%
#   group_by(Contact_Year) %>%
#   summarise(prop.young=mean(YOUNG),n.young=sum(YOUNG),n.old=sum(OLD)) %>%
#   left_join(goodyears, by="Contact_Year") %>%
#   filter(Contact_Year>2003) 
# dim(youngbreeders)
# 
# 
# ### analysis of trend over time
# 
# # m2qeff<-glm(cbind(n.young,n.old)~Contact_Year+I(Contact_Year^2), data=youngbreeders, family="binomial", weights=prop.seen)
# # summary(m2qeff)
# 
# m2leff<-glm(cbind(n.young,n.old)~Contact_Year, data=youngbreeders, family="binomial", weights=prop.seen)
# summary(m2leff)
# 
# 
# ### prediction of effect size
# newdat<-data.frame(Contact_Year=seq(2004,2021,1),prop.seen=1)
# 
# ## grad the inverse link function
# ilink <- family(m2leff)$linkinv
# ## add fit and se.fit on the **link** scale
# newdat <- bind_cols(newdat, setNames(as_tibble(predict(m2leff, newdat, se.fit = TRUE)[1:2]),
#                                      c('fit_link','se_link')))
# ## create the interval and backtransform
# newdat <- mutate(newdat,
#                  pred.prop  = ilink(fit_link),
#                  ucl = ilink(fit_link + (1.96 * se_link)),
#                  lcl = ilink(fit_link - (1.96 * se_link)))




### does the proportion of OLD breeders change over time?

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


### analysis of trend over time

m2oleff<-glm(cbind(n.old,n.young)~Contact_Year, data=oldbreeders, family="binomial", weights=prop.seen)
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

### COMBINE PROPORTION OF YOUNG AND OLD BREEDERS IN ONE PLOT
## abandoned on 2 Oct 2021 because only old breeder proportion shown in manuscript
# colors <- c("young (< 10 years old)" = "steelblue", "old (>30 years old)" = "indianred")
# ggplot(newdat) +
#   geom_line(aes(x=Contact_Year, y=pred.prop,colour = "young (< 10 years old)")) +
#   geom_ribbon(aes(x=Contact_Year, ymin=lcl, ymax=ucl),fill = "steelblue", alpha = 0.2) +
#   geom_point(data=youngbreeders,aes(x=Contact_Year, y=prop.young,colour = "young (< 10 years old)"), size=3) +
#   
#   geom_line(data=olddat,aes(x=Contact_Year, y=pred.prop,colour = "old (>30 years old)")) +
#   geom_ribbon(data=olddat,aes(x=Contact_Year, ymin=lcl, ymax=ucl), fill= "indianred",alpha = 0.2) +
#   geom_point(data=oldbreeders,aes(x=Contact_Year, y=prop.old,colour = "old (>30 years old)"), size=3) +
#   labs(x = "Year",
#        y = "Annual proportion of breeders",
#        color = "Age group") +
#   scale_color_manual(values = colors) +
#   scale_y_continuous(breaks=seq(0,0.4,0.1), limits=c(0,0.4))+
#   scale_x_continuous(breaks=seq(2005,2021,2), limits=c(2004,2021))+
#   
#   
#   ### add the bird icons
#   annotation_custom(TRALicon, xmin=2013.5, xmax=2015.5, ymin=0.32, ymax=0.42) +
#   
#   theme(panel.background=element_rect(fill="white", colour="black"), 
#         axis.text=element_text(size=18, color="black"), 
#         axis.title=element_text(size=20),
#         legend.title=element_text(size=18),
#         legend.text=element_text(size=16),
#         legend.background=element_blank(),
#         legend.key=element_blank(),
#         legend.position=c(0.80, 0.90),
#         panel.grid.minor = element_blank())
# 
# ggsave("C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\TRAL_IPM\\Fig2.jpg", width=9, height=6)


### COMBINE PROPORTION OF YOUNG AND OLD BREEDERS IN ONE PLOT
ggplot(olddat) +

  geom_line(aes(x=Contact_Year, y=pred.prop),colour = "indianred") +
  geom_ribbon(aes(x=Contact_Year, ymin=lcl, ymax=ucl), fill= "indianred",alpha = 0.2) +
  geom_point(data=oldbreeders,aes(x=Contact_Year, y=prop.old), size=3) +
  labs(x = "Year",
       y = "Annual proportion of old breeders") +
  scale_y_continuous(breaks=seq(0,0.15,0.03), limits=c(0,0.15))+
  scale_x_continuous(breaks=seq(2005,2021,2), limits=c(2004,2021))+
  
  
  ### add the bird icons
  annotation_custom(TRALicon, xmin=2004.5, xmax=2006.5, ymin=0.12, ymax=0.15) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        panel.grid.minor = element_blank())

ggsave("C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\TRAL_IPM\\Fig2.jpg", width=9, height=6)







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



