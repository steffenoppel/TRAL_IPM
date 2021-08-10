#####################################################################################
######## EXAMINATION OF BREEDING AGE OF TRISTAN ALBATROSS ON GOUGH     ##############
#####################################################################################

### written by Steffen Oppel in August 2021
### steffen.oppel@rspb.org.uk
### examines 3 question re TRAL breeding age:
## 1. does average age of all TRAL individuals contacted in a given year change over time?
## 2. does average age of all TRAL breeders in a given year change over time?
## 3. does average age of all TRAL first recruiters in a given year change over time? [This is potentially effort-dependent]


library(tidyverse)
library(lubridate)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select


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
##   QUESTION 1: does age of all recorded bird change over time? ###############
#############################################################################
head(contacts)
dim(contacts)
returns<-contacts %>%
  filter(!(FIRST_YEAR==Contact_Year)) %>%   ## potentially change this to remove only ringed chicks? Age %in% c("Chick","Fledgling")
  left_join(goodyears, by="Contact_Year") %>%
  select(ContactID,BirdID,Contact_Year,Age,Sex,Nest_Description,FIRST_AGE,FIRST_YEAR,ContAge,n,N_marked,prop.seen) %>%
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



#############################################################################
##   QUESTION 2: does age of breeders change over time? ###############
#############################################################################

breeders<-contacts %>%
  filter(!(FIRST_YEAR==Contact_Year)) %>%   ## potentially change this to remove only ringed chicks? Age %in% c("Chick","Fledgling")
  filter(!(is.na(Nest_Description))) %>%
  left_join(goodyears, by="Contact_Year") %>%
  select(ContactID,BirdID,Contact_Year,Age,Sex,Nest_Description,FIRST_AGE,FIRST_YEAR,ContAge,n,N_marked,prop.seen) %>%
  filter(Contact_Year>2009) 
dim(breeders)


### exploratory plots
ggplot(breeders) +
  geom_point(aes(x=Contact_Year, y=ContAge)) +
  geom_smooth(aes(x=Contact_Year, y=ContAge),method="lm")

ggplot(breeders) +
  geom_histogram(aes(x=ContAge))

### analysis
m2eff<-glm(ContAge~Contact_Year, data=breeders, family="poisson",weights=prop.seen)
summary(m2eff)





#############################################################################
##   QUESTION 3: does age of first returners change over time? ###############
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
  filter(Contact_Year>2003) 
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




# LOAD AND MANIPULATE ICONS
library(grid)
library(magick)
imgTRAL<-image_read("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\PR_Comms\\Icons\\alby 4.jpg") %>% image_transparent("white", fuzz=5)
TRALicon <- rasterGrob(imgTRAL, interpolate=TRUE)

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

ggsave("TRAL_recruit_age_time.jpg", width=9, height=6)
