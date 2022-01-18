##########################################################################
#
# TRISTAN ALBATROSS INTEGRATED POPULATION MODEL 2001-2050
#
##########################################################################
# based on Kery and Schaub 2012, Chapter 11
# modified by Steffen oppel, JULY 2019
# based on AYNA IPM code
# major changes to AYNA: more gappy count data series
# adult encounter occasions span almost entire year
# breeding success is area specific

# v2 initiated on 17 July 2019 after advice from Adam Butler to include multinomial probability for N.est
# dmulti did not work with 'cannot find sampler' error, so changed to obs model instead
# v3 initiated on 17 July 2019 after advice from Chris Jones that we cannot exclude sites - lumped all data together and used a single input vector for count data, thus reducing the problem of dmulti to include multinomial probability for N.est
# seems to work but when niter>500 error occurs 'slicer stuck at infinite density'
# updated 22 July 2019: inserted truncations T(0.001,0.999) to prevent slice sampler getting stuck at infinite density
## UPDATED on 24 July 2019 to include future projection
## UPDATED ON 6 August 2019 to revise SSM part of the model after nonconvergence of basic SSM
## UPDATED ON 2 April 2020 to fix juvenile survival to constant and fix first pop count data to a lower number (based on Peter Ryan's suggestion via email)
## UPDATED ON 21 April 2020 to include future scenarios - mouse eradication (fec increases) or mouse predation (surv decreases)
## TRAL_IPM_fut_v3.jags runs ok and produces sensible output
## TRAL_IPM_fut_v4.jags intended to include 3 future scenarios
## UPDATED ON 20 May 2020 to include carrying capacity to prevent overshoot to 10,000 pairs
## TRAL_IPM_fut_v5.jags includes carrying capacity at 2500 birds

## v5 completely revised in January 2021 to update survival estimation
## abandoned individual survival via SSM and instead use m-array approach with two age groups
## REMOVED 2001 data and start process in 2004!

## 9 Jan 2021: tried to relate pop change (skip.prob) to proportion of nests failing at egg vs chick stage (or day of failure)
## none of the breed.fail metrics explains pop change from one year to the next - need to keep this totally random

## 9 JANUARY 2021: major revision to population process trying to simplify the quantities that explain which birds return to island
## 10-12 Jan 2021: simplified initial values and adjusted priors - MODEL CONVERGED BUT has very high juvenile survival bounded by priors on p and phi

## 13 January 2021: re-inserted carrying capacity into future projections to prevent pop sizes >5000
## 15 January 2021: removed data from 2011 and included biennial breeding schedule into future projections

## 20 January 2021: preliminary output suggests slower increase with than without eradication - examined code to identify problem
## 21 January 2021: added Ntot into model and as output

## NEED TO DO: calculate lambda for total pop - not just breeding pop (include imm and at sea birds)
## revise transition to future to avoid the massive spike in 2022 in scenarios 1 and 3 (due to succ breeders remaining at sea?)
## calculate total pop at each time step and include in output

## 24 JANUARY 2021: MAJOR revision to data preparation and survival part of the model - included two p-intercepts and set p.juv to 0 for first year
## 25 January 2021: improved calculation of mean.p.ad to be used in future projection

## FINAL IS BASED ON V5 but includes data from 2021 and eliminates all unnecessary code
## updated 7 Feb 2021 to include age-dependent recapture probability

## updated 10 February to inform initial population status (Ntot[1]) based on counts of adults and chicks in 1999-2003

## MODEL CONVERGED ON 11 FEBRUARY 2021 with sensible estimates

## survival model requires 'agebeta' prior to be set to 0.5 or higher, otherwise it will not converge

## short version initiated 5 June 2021: reduce data series by another 4 years and start in 2008 rather than 2004 to improve data quality
## short version redacted 18 June 2021 - did not improve convergence

## 1 July 2021: revised input data (eliminated duplicate chick counts in 2011 and 2012)

## 4 AUGUST 2021: revised priors in model based on Beth Clark's suggestions

## 22 DECEMBER 2021: revision of manuscript started and included 2021 fecundity data

## 07 January 2022: re-run model to report IM, and fix Ntot.breed to never fall below actually counted value

## 17 January 2022: revised model to 3 detection parameters and included a max for breed.prop

library(tidyverse)
library(lubridate)
library(data.table)
library(jagsUI)
library(runjags)   ## added by Beth in July 2021 because jagsUI would not converge
filter<-dplyr::filter
select<-dplyr::select



#########################################################################
# LOAD PRE-PREPARED DATA ON COUNTS AND BREEDING SUCCESS
#########################################################################
### see 'IPM_DATA_PREPARATION_TRAL.R' for details on how data are aggregated

### NOTE THAT SORT ORDER OF GONYDALE AND GP VALLEY HAS SHIFTED ON 15 Jan 2021 (due to switch to if_else on R4.0.2)

## LOAD PREPARED M-ARRAY FOR SURVIVAL ESTIMATION
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
load("TRAL_IPM_input.marray.REV2021.RData")

## BOTH ARRAYS MUST HAVE EXACT SAME DIMENSIONS
dim(chick.marray)
dim(adult.marray)


### COUNT DATA FOR POPULATION TREND ######
TRAL.pop<-POPSIZE ##fread("TRAL_INCU_counts_2001_2020.csv")
head(TRAL.pop)
names(TRAL.pop)
TRAL.pop[14,4]<-136   ### number of nests monitored in Gonydale that year

## COMBINE SITES THAT WERE AMBIGUOUSLY DEFINED OVER TIME 
TRAL.pop<-TRAL.pop %>% gather(key='Site', value='Count',-Year) %>%
  filter(Year>2003) %>%   ### changed from 2003
  mutate(Site=if_else(Site %in% c('GP Valley','West Point'),'GP Valley',as.character(Site))) %>%
  mutate(Site=if_else(Site %in% c('Gonydale','Green Hill','Hummocks'),'Gonydale',as.character(Site))) %>%
  group_by(Year,Site) %>%
  summarise(Count=sum(Count, na.rm=T)) %>%
  mutate(Count=ifelse(Count==0,NA,Count)) %>%
  spread(key=Site, value=Count)

## to account for data gaps, need to quantify what proportion of the population was counted

TRAL.props<-prop.table(as.matrix(TRAL.pop[,2:9]),1)
#mean.props<-apply(TRAL.props[c(1,4,6:10,12:16,18:19),],2,mean) ## for start in 2001
mean.props<-apply(TRAL.props[c(1,3:7,9:13,15:17),],2,mean) ## for start in 2004
#mean.props<-apply(TRAL.props[c(1:3,5:9,11:14),],2,mean) ## for start in 2008

TRAL.pop$prop.counted<-0
for (l in 1:length(TRAL.pop$Year)){
  TRAL.pop$prop.counted[l]<-sum(mean.props[which(!is.na(TRAL.pop[l,2:9]))])
}


## CALCULATE SUM PER YEAR 
TRAL.pop$tot<-rowSums(TRAL.pop[,2:9], na.rm=T)
#TRAL.pop$tot[c(2,3,11)]<-NA
R<- as.matrix(TRAL.pop[,2:9])
n.years<-nrow(R)
n.sites<-ncol(R)



### PLOT TO SPOT ANY OUTLIERS OF BCOUNTS
ggplot(TRAL.pop, aes(x=Year,y=tot)) +geom_point(size=2, color='darkred')+geom_smooth(method='lm') 




#### BREEDING SUCCESS DATA FOR FECUNDITY ######
TRAL.chick<-CHICKCOUNT  ##fread("TRAL_CHIC_counts_2001_2020.csv")
## detailed nest data from Gonydale only available for some years - we use this to replace missing count data
TRAL.bs<-FECUND ##fread("TRAL_breed_success_2006_2019.csv")
TRAL.chick[TRAL.chick$Year==2013,4]<-as.integer(TRAL.pop[TRAL.pop$Year==2013,4]*TRAL.bs$BREED_SUCC[TRAL.bs$Year==2013])
TRAL.chick[TRAL.chick$Year==2014,4]<-as.integer(TRAL.pop[TRAL.pop$Year==2014,4]*TRAL.bs$BREED_SUCC[TRAL.bs$Year==2014])

TRAL.chick<-TRAL.chick %>% gather(key='Site', value='Count',-Year) %>%
  filter(Year>2003) %>%  ### changed from 2003
  mutate(Site=if_else(Site %in% c('GP Valley','West Point'),'GP Valley',as.character(Site))) %>%
  mutate(Site=if_else(Site %in% c('Gonydale','Green Hill','Hummocks'),'Gonydale',as.character(Site))) %>%
  group_by(Year,Site) %>%
  summarise(Count=sum(Count, na.rm=T)) %>%
  mutate(Count=ifelse(Count==0,NA,Count)) %>%
  spread(key=Site, value=Count)

### NOTE THAT USE OF 'if_else' switches Gonydale and GP_Valley columns around [only happens in R4.0.2!!]
TRAL.chick[4,4]<-CHICKCOUNT[11,4] ### in 2011 no adults were counted in Green hill and Hummocks, so we cannot add up the chicks across those 3 sites


TRAL.chick$tot<-rowSums(TRAL.chick[,2:9], na.rm=T)
#TRAL.chick$tot[c(3,13,19)]<-NA   ## when start in 2001
TRAL.chick$tot[10]<-NA   ## when start in 2004
#TRAL.chick$tot[6]<-NA   ## when start in 2008

J<- as.matrix(TRAL.chick[,2:9])

### specify constants for JAGS
n.years<-dim(R)[1]		## defines the number of years
n.sites<-dim(R)[2]    ## defines the number of study areas


### UPDATE 10 January 2021 - reduce R and J to vectors of sum across the study areas for which we have data

Jlong<-TRAL.chick %>% gather(key='Site', value="chicks",-Year)
PROD.DAT<-TRAL.pop %>% select(-prop.counted,-tot) %>% gather(key='Site', value="adults",-Year) %>%
  left_join(Jlong, by=c("Year","Site")) %>%
  mutate(include=ifelse(is.na(adults+chicks),0,1)) %>%
  filter(include==1) %>%
  group_by(Year) %>%
  summarise(J=sum(chicks),R=sum(adults))





### DIMENSION MISMATCH IN DATA
# IPM runs from 2004-2021 # changed on 7 February - and again on 4 June 2021
# survival analysis runs from 1979-2021, but recapture index refers to columns, which represent year 1980-2021 plus the ones never recaptured (last column)
# very difficult
names(TRAL_CHICK)
TRAL.pop$Year

OFFSET<-min(which(!is.na(match(as.numeric(names(TRAL_CHICK)[2:44]),TRAL.pop$Year))))
names(TRAL_CHICK)[OFFSET+1]
TRAL.pop$Year[1]



#########################################################################
# SPECIFY FUTURE DECREASE IN SURVIVAL
#########################################################################

dec.surv=0.9  ## we assume that adult survival will decrease by 10%
lag.time=10    ## the decrease will take 10 years to materialise
PROJECTION.years<-seq(1,30,1)  ## we specify the relative survival decrease for all 30 years in the projection

fut.surv.change<- expand.grid(PROJECTION.years,dec.surv,lag.time) %>%
  rename(Year=Var1,SURV3=Var2,LAG=Var3) %>%
  mutate(ann.offset=(SURV3-1)/LAG) %>%
  mutate(SURV3=ifelse(Year<LAG,1+(Year*ann.offset),SURV3)) %>%
  mutate(SURV1=1,SURV2=1) %>%
  select(Year, SURV1,SURV2,SURV3)
  


#########################################################################
# SPECIFY DETECTION PROBABILITY
#########################################################################
# use 3 rather than just 2 parameters, with higher probability after 2009
goodyears$p.sel<-ifelse(goodyears$p.sel==2 & goodyears$Contact_Year>2008,3,goodyears$p.sel)



#########################################################################
# SPECIFY MODEL IN JAGS
#########################################################################
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
sink("TRAL_IPM_FINAL_REV2022.jags")
cat("

model {
    #-------------------------------------------------
    # integrated population model for the Gough TRAL population
    # - age structured model with 10 age classes 
    # - adult survival based on CMR ringing data
    # - pre breeding census, female-based assuming equal sex ratio & survival
    # - productivity based on all areas incu and chick counts
    # - simplified population process with assumed age at recruiting = 10 (Wanless et al. 2009)
    # - adult breeders skipping when unsuccessful at rate of 22-32%, all successful breeders will skip
    # - linked population process with SUM OF count data
    # - v4 includes 3 scenarios of future projection: no change, improved fecundity, reduced adult survival
    # - marray_v1 uses marray for survival estimation to speed up computation time
    # - marray_simplified adjusts population process to make number of bird returning completely random
    # - marray_simplified_v2 further simplifies return process to tie to recapture probability
    # - marray_simplified_v3 includes carrying capacity for future projections and re-inserts breeding success
    # - marray_simplified_v4 calculates Ntot and lambda based on Ntot (rather than just breeding pop).
    # - marray_simplified_v5 includes two values for p (high and low monitoring effort years) and sets p.juv to 0 for first year.
    # - marray_simplified_v5 also calculates mean propensity and recruitment from annual values for good monitoring years.
    # - marray_age_recruit allows for varying juvenile recapture probability with age.
    # - marray_age_recruit_immat creates a loop over immature age classes rather than specify each age separately
    # - marray_age_recruit_immat_2008 sets juv surv to mean and starts in 2008 rather than 2004
    # -------------------------------------------------
    
#-------------------------------------------------  
# 1. PRIORS FOR ALL DATA SETS
#-------------------------------------------------
    
    
    # -------------------------------------------------        
    # 1.1. Priors and constraints FOR FECUNDITY
    # -------------------------------------------------
    
    for (t in 1:T){  
      ann.fec[t] ~ dbeta(32,68) ## dnorm(0.32,10) T(0.001,0.999)        ## Informative Priors on fecundity based on Wanless et al 2009
    } #t
    
    
    # -------------------------------------------------        
    # 1.2. Priors and constraints FOR POPULATION COUNTS
    # -------------------------------------------------
    for (s in 1:n.sites){			### start loop over every study area
      for (t in 1:T){			### start loop over every year
        sigma.obs[s,t] ~ dexp(0.1)	#Prior for SD of observation process (variation in detectability)
        tau.obs[s,t]<-pow(sigma.obs[s,t],-2)
      }
    }
    
    
    # -------------------------------------------------        
    # 1.3. Priors and constraints FOR SURVIVAL
    # -------------------------------------------------
    
    ### RECAPTURE PROBABILITY
    for (gy in 1:3){    ## for good and poor monitoring years
      mean.p.juv[gy] ~ dunif(0, 1)	           # Prior for mean juvenile recapture - should be higher than 20% if they survive!
      mean.p.ad[gy] ~ dunif(0, 1)	           # Prior for mean adult recapture - should be higher than 20%
      mu.p.juv[gy] <- log(mean.p.juv[gy] / (1-mean.p.juv[gy])) # Logit transformation
      mu.p.ad[gy] <- log(mean.p.ad[gy] / (1-mean.p.ad[gy])) # Logit transformation
    }
    agebeta ~ dunif(0,1)    # Prior for shape of increase in juvenile recapture probability with age
    
    ## RANDOM TIME EFFECT ON RESIGHTING PROBABILITY OF JUVENILES
    for (t in 1:(n.occasions-1)){
      for (j in 1:t){ ## zero by definition (these are never actually used)
        p.juv[t,j] <- 0
      }
      for (j in (t+1):(n.occasions-1)){
        logit(p.juv[t,j])  <- mu.p.juv[goodyear[j]] + agebeta*(j - t) + eps.p[j]
      }
    }
    
    ## PRIORS FOR RANDOM EFFECTS
    sigma.p ~ dexp(1)                # Prior for standard deviation
    tau.p <- pow(sigma.p, -2)
    
    
    ### SURVIVAL PROBABILITY
    mean.phi.juv ~ dbeta(75.7,24.3)             # Prior for mean juvenile survival first year 0.757, second year 0.973 in Laysan albatross
    mean.phi.ad ~ dbeta(91,9)              # Prior for mean adult survival - should be higher than 70%
    mu.juv <- log(mean.phi.juv / (1-mean.phi.juv)) # Logit transformation
    mu.ad <- log(mean.phi.ad / (1-mean.phi.ad)) # Logit transformation
    
    ## PRIORS FOR RANDOM EFFECTS
    sigma.phi ~ dexp(1)                # Prior for standard deviation
    tau.phi <- pow(sigma.phi, -2)
    
    ## RANDOM TIME EFFECT ON SURVIVAL AND ADULT RECAPTURE
    for (j in 1:(n.occasions-1)){
      logit(phi.juv[j]) <- mu.juv + eps.phi[j]*juv.poss[j]
      logit(phi.ad[j]) <- mu.ad + eps.phi[j]
      eps.phi[j] ~ dnorm(0, tau.phi) 
      logit(p.ad[j])  <- mu.p.ad[goodyear[j]] + eps.p[j]    #### CAT HORSWILL SUGGESTED TO HAVE A CONTINUOUS EFFORT CORRECTION: mu.p.ad + beta.p.eff*goodyear[j] + eps.p[j]
      eps.p[j] ~ dnorm(0, tau.p)
    }
    
    
    
#-------------------------------------------------  
# 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
#-------------------------------------------------
    
    # -------------------------------------------------        
    # 2.1. System process: female based matrix model
    # -------------------------------------------------
    
    ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1 - based on deterministic multiplications
    ## ADJUSTED BASED ON PAST POPULATION SIZES WITH CHICK COUNTS SINCE 1999
    
    IM[1,1,1] ~ dnorm(324,20) T(0,)                                 ### number of 1-year old survivors is low because few chicks hatched in 2003 - CAN BE MANIPULATED
    IM[1,1,2] <- 0
    IM[1,1,3] <- IM[1,1,1] - IM[1,1,2]
    
    IM[1,2,1] ~ dnorm(257,20) T(0,)                                  ### number of 2-year old survivors is very low because very few chicks hatched in 2002 - CAN BE MANIPULATED
    IM[1,2,2] <- IM[1,2,1]*p.juv.recruit.f[2]
    IM[1,2,3] <- IM[1,1,1] - IM[1,1,2]
    
    IM[1,3,1] ~ dnorm(462,20) T(0,)                                 ### number of 3-year old survivors is higher because many chicks hatched in 2001 - CAN BE MANIPULATED
    IM[1,3,2] <- IM[1,3,1]*p.juv.recruit.f[3]
    IM[1,3,3] <- IM[1,1,1] - IM[1,1,2]
    
    IM[1,4,1] ~ dnorm(207,20) T(0,)                                 ### number of 4-year old survivors is very low because few chicks hatched in 2000 - CAN BE MANIPULATED
    IM[1,4,2] <- IM[1,4,1]*p.juv.recruit.f[4]
    IM[1,4,3] <- IM[1,1,1] - IM[1,1,2]
    
    IM[1,5,1] ~ dnorm(700,10) T(0,)                                  ### number of 5-year old survivors is huge because a lot of chicks hatched in 1999  - CAN BE MANIPULATED
    IM[1,5,2] <- IM[1,5,1]*p.juv.recruit.f[5]
    IM[1,5,3] <- IM[1,1,1] - IM[1,1,2]
    
    IM[1,6,1] ~ dnorm(225,20) T(0,)                                 ### very uncertain number of of 6-year old survivors because no data from 1998 or previously  - CAN BE MANIPULATED
    IM[1,6,2] <- IM[1,6,1]*p.juv.recruit.f[6]
    IM[1,6,3] <- IM[1,1,1] - IM[1,1,2]
    
    for(age in 7:30) {
    IM[1,age,1] ~ dbin(pow(mean.phi.ad,(age-1)), round(IM[1,age-1,3]))
    IM[1,age,2] <- IM[1,age,1]*p.juv.recruit.f[age]
    IM[1,age,3] <- IM[1,age,1] - IM[1,age,2]
    }
    N.recruits[1] <- sum(IM[1,,2])  ### number of this years recruiters - irrelevant in year 1 as already included in Ntot.breed prior
    
    Ntot.breed[1] ~ dnorm(1869,100) T(0,)  ### sum of counts is 1869
    JUV[1] ~ dnorm(510,100) T(0,)          ### sum of chicks is 510
    N.atsea[1] ~ dnorm(530,20) T(0,)    ### unknown number - CAN BE MANIPULATED
    Ntot[1]<-sum(IM[1,,3]) + Ntot.breed[1]+N.atsea[1]  ## total population size is all the immatures plus adult breeders and adults at sea - does not include recruits in Year 1
    breed.prop[1] <- Ntot.breed[1]/Ntot[1]
    
    ### FOR EVERY SUBSEQUENT YEAR POPULATION PROCESS
    
    for (tt in 2:T){
    
    ## THE PRE-BREEDING YEARS ##
    
    ## define recruit probability for various ages ##
    for (age in 1:30) {
    logit(p.juv.recruit[age,tt])<-mu.p.juv[2] + eps.p[tt+24] + (agebeta * age)
    }
    
    
    ## IMMATURE MATRIX WITH 3 columns:
    # 1: survivors from previous year
    # 2: recruits in current year
    # 3: unrecruited in current year (available for recruitment next year)
    
    nestlings[tt] <- ann.fec[tt] * 0.5 * Ntot.breed[tt]                                                     ### number of locally produced FEMALE chicks
    JUV[tt] ~ dpois(nestlings[tt])                                                                     ### need a discrete number otherwise dbin will fail, dpois must be >0
    IM[tt,1,1] ~ dbin(phi.juv[tt+24], max(1,round(JUV[tt-1])))                                  ### number of 1-year old survivors 
    IM[tt,1,2] <- 0
    IM[tt,1,3] <- IM[tt,1,1] - IM[tt,1,2]
    
    for(age in 2:30) {
    IM[tt,age,1] ~ dbin(phi.ad[tt+24], max(1,round(IM[tt-1,age-1,3])))
    IM[tt,age,2] <- min(round(IM[tt,age-1,3]),IM[tt,age,1])*p.juv.recruit[age,tt]
    IM[tt,age,3] <- IM[tt,age,1] - IM[tt,age,2]
    }
    N.recruits[tt] <- sum(IM[tt,,2])  ### number of this years recruiters
    
    
    ## THE BREEDING POPULATION ##
    # Ntot.breed comprised of first-time breeders, previous skippers, and previous unsuccessful breeders
    # simplified in simplified_v2 to just adult survivors with p.ad as proportion returning
    
    N.ad.surv[tt] ~ dbin(phi.ad[tt+24], round(Ntot.breed[tt-1]+N.atsea[tt-1]))           ### previous year's adults that survive
    breed.prop[tt] <- max(p.ad[tt+24],mean.p.ad[3])						### p.ad is only equivalent to breeding propensity in years with high effort
    N.breed.ready[tt] ~ dbin(breed.prop[tt], N.ad.surv[tt])                  ### number of available breeders is proportion of survivors that returns
    Ntot.breed[tt]<- round(N.breed.ready[tt]+N.recruits[tt])              ### number of counted breeders is sum of old breeders returning and first recruits
    N.atsea[tt] <- round(N.ad.surv[tt]-(Ntot.breed[tt]-N.recruits[tt]))                     ### potential breeders that remain at sea    
    
    ### THE TOTAL TRAL POPULATION ###
    Ntot[tt]<-sum(IM[tt,,3]) + Ntot.breed[tt]+N.atsea[tt]  ## total population size is all the immatures plus adult breeders and adults at sea
    
    } # tt
    
    


    
    # -------------------------------------------------        
    # 2.2. Observation process for population counts: state-space model of annual counts
    # -------------------------------------------------
    
    for (s in 1:n.sites){			### start loop over every study area
    
      ## Observation process
    
      for (t in 1:T){
        y.count[t,s] ~ dnorm(Ntot.breed[t]*prop.sites[s], tau.obs[s,t])								# Distribution for random error in observed numbers (counts)
      }														# run this loop over t= nyears
    }		## end site loop
    
    
    # -------------------------------------------------        
    # 2.3. Likelihood for fecundity: Logistic regression from the number of surveyed broods
    # -------------------------------------------------

    for (t in 1:(T-1)){
      J[t] ~ dbin(ann.fec[t], R[t])
    } #	close loop over every year in which we have fecundity data
    
    
    
    # -------------------------------------------------        
    # 2.4. Likelihood for adult and juvenile survival from CMR
    # -------------------------------------------------
    
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
      marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], r.j[t])
      marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], r.a[t])
    }
    
    
    # Define the cell probabilities of the m-arrays
    # Main diagonal
    for (t in 1:(n.occasions-1)){
      q.ad[t] <- 1-p.ad[t]            # Probability of non-recapture
    
      for(j in 1:(n.occasions-1)){
        q.juv[t,j] <- 1 - p.juv[t,j]
      }    
    
      pr.j[t,t] <- 0
      pr.a[t,t] <- phi.ad[t]*p.ad[t]
    
      # Above main diagonal
      for (j in (t+1):(n.occasions-1)){
        pr.j[t,j] <- phi.juv[t]*prod(phi.ad[(t+1):j])*prod(q.juv[t,t:(j-1)])*p.juv[t,j]
        pr.a[t,j] <- prod(phi.ad[t:j])*prod(q.ad[t:(j-1)])*p.ad[j]
      } #j
    
      # Below main diagonal
      for (j in 1:(t-1)){
        pr.j[t,j] <- 0
        pr.a[t,j] <- 0
      } #j
    } #t
    
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
      pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
      pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
    } #t
    
    
#-------------------------------------------------  
# 3. DERIVED PARAMETERS FOR OUTPUT REPORTING
#-------------------------------------------------
    
    
    ## DERIVED POPULATION GROWTH RATE PER YEAR
    for (t in 1:(T-1)){
      lambda[t]<-Ntot[t+1]/max(1,Ntot[t])  ## division by 0 creates invalid parent value
    }		## end year loop
    
    ## DERIVED MEAN FECUNDITY 
    mean.fec <- mean(ann.fec)
    pop.growth.rate <- exp((1/(T-1))*sum(log(lambda[1:(T-1)])))   # Geometric mean
    
    
#-------------------------------------------------  
# 4. PROJECTION INTO FUTURE
#-------------------------------------------------
  ## includes 3 scenarios
  ## scenario 1: projection with no changes in demography
  ## scenario 2: successful mouse eradication in 2021 - fecundity doubles
  ## scenario 3: increasing mouse impacts on adult survival (adult survival decreases by 10%)
    
    ## recruit probability
    for (age in 1:30) {
      logit(p.juv.recruit.f[age])<-mu.p.juv[2] + (agebeta * age)
    }


  # -------------------------------------------------        
  # 4.1. System process for future
  # -------------------------------------------------
        
  ## LOOP OVER EACH SCENARIO  
  for(scen in 1:n.scenarios){
    
    ### ~~~~~~~~~~ COPY POPULATIONS FROM LAST YEAR OF DATA SERIES FOR FIRST FUTURE YEAR ~~~~~~~~~###
    
    ## IMMATURE MATRIX WITH 3 columns:
    # 1: survivors from previous year
    # 2: recruits in current year
    # 3: unrecruited in current year (available for recruitment next year)

    nestlings.f[scen,1] ~ dbin(fut.fec.change[scen]*mean.fec*0.5,round(Ntot.breed.f[scen,1]))                      ### number of locally produced FEMALE chicks based on average fecundity 
    IM.f[scen,1,1,1] ~ dbin(mean.phi.juv, max(1,round(JUV[T])))                                  ### number of 1-year old survivors 
    IM.f[scen,1,1,2] <- 0
    IM.f[scen,1,1,3] <- IM.f[scen,1,1,1] - IM.f[scen,1,1,2]
    
    for(age in 2:30) {
      IM.f[scen,1,age,1] ~ dbin(mean.phi.ad, max(1,round(IM[T,age-1,3])))
      IM.f[scen,1,age,2] <- min(round(IM[T,age-1,3]),IM.f[scen,1,age,1])*p.juv.recruit.f[age]
      IM.f[scen,1,age,3]   <- IM.f[scen,1,age,1] - IM.f[scen,1,age,2]
    }
    N.recruits.f[scen,1] <- sum(IM.f[scen,1,,2])  ### number of this years recruiters
    
    N.ad.surv.f[scen,1] ~ dbin(mean.phi.ad, round(Ntot.breed[T]+N.atsea[T]))              ### previous year's adults that survive
    N.breed.ready.f[scen,1] ~ dbin(mean.p.ad[2], round(N.ad.surv.f[scen,1]))              ### number of available breeders is proportion of survivors that returns, with fecundity INCLUDED in return probability
    Ntot.breed.f[scen,1]<- round(N.breed.ready.f[scen,1]+N.recruits.f[scen,1])            ### number of counted breeders is sum of old breeders returning and first recruits
    N.atsea.f[scen,1] <- round(N.ad.surv.f[scen,1]-N.breed.ready.f[scen,1])               ### potential breeders that remain at sea
    N.succ.breed.f[scen,1] ~ dbin(mean.fec, round(Ntot.breed.f[scen,1]))                  ### these birds will  remain at sea because they bred successfully
    
    ### THE TOTAL TRAL POPULATION ###
    Ntot.f[scen,1]<-sum(IM.f[scen,1,,3])+Ntot.breed.f[scen,1]+N.atsea.f[scen,1]  ## total population size is all the immatures plus adult breeders and adults at sea


    
    ### ~~~~~~~~~~ LOOP OVER ALL SUBSEQUENT FUTURE YEARS ~~~~~~~~~###

    for (tt in 2:FUT.YEAR){
    
  
      ## INCLUDE CARRYING CAPACITY OF 2500 breeding pairs (slightly more than maximum ever counted)
      carr.capacity[scen,tt] ~ dnorm(2500,5) T(0,)
    
      ## THE PRE-BREEDING YEARS ##
      ## because it goes for 30 years, all pops must be safeguarded to not become 0 because that leads to invald parent error
    
      ## IMMATURE MATRIX WITH 3 columns:
      # 1: survivors from previous year
      # 2: recruits in current year
      # 3: unrecruited in current year (available for recruitment next year)
      nestlings.f[scen,tt] ~ dbin(fut.fec.change[scen]*mean.fec*0.5,round(Ntot.breed.f[scen,tt]))                       ### number of locally produced FEMALE chicks based on average fecundity
      IM.f[scen,tt,1,1] ~ dbin(mean.phi.juv, max(1,round(nestlings.f[scen,tt-1])))                                  ### number of 1-year old survivors 
      IM.f[scen,tt,1,2] <- 0
      IM.f[scen,tt,1,3] <- IM.f[scen,tt,1,1] - IM.f[scen,tt,1,2]
    
      for(age in 2:30) {
        IM.f[scen,tt,age,1] ~ dbin(mean.phi.ad, max(1,round(IM.f[scen,tt-1,age-1,3])))
        IM.f[scen,tt,age,2] <- min(round(IM.f[scen,tt-1,age-1,3]),IM.f[scen,tt,age,1])*p.juv.recruit.f[age]
        IM.f[scen,tt,age,3] <- IM.f[scen,tt,age,1] - IM.f[scen,tt,age,2]
      }
      N.recruits.f[scen,tt] <- sum(IM.f[scen,tt,,2])  ### number of this years recruiters
    
      ## THE BREEDING POPULATION ##
      N.ad.surv.f[scen,tt] ~ dbin(fut.surv.change[tt,scen]*mean.phi.ad, round((Ntot.breed.f[scen,tt-1]-N.succ.breed.f[scen,tt-1])+N.atsea.f[scen,tt-1]))           ### previous year's adults that survive
      N.prev.succ.f[scen,tt] ~ dbin(fut.surv.change[tt,scen]*mean.phi.ad, round(N.succ.breed.f[scen,tt-1]))                  ### these birds will  remain at sea because tey bred successfully
      N.breed.ready.f[scen,tt] ~ dbin(min(0.99,(mean.p.ad[2])), max(1,round(N.ad.surv.f[scen,tt])))                  ### number of available breeders is proportion of unsuccessful or non-breeding survivors that returns
      Ntot.breed.f[scen,tt]<- min(carr.capacity[scen,tt],round(N.breed.ready.f[scen,tt]+N.recruits.f[scen,tt]))              ### number of counted breeders is sum of old breeders returning and first recruits
      N.succ.breed.f[scen,tt] ~ dbin(fut.fec.change[scen]*mean.fec, round(Ntot.breed.f[scen,tt]))                  ### these birds will  remain at sea because they bred successfully
      N.atsea.f[scen,tt] <- round(N.ad.surv.f[scen,tt]-N.breed.ready.f[scen,tt]+N.prev.succ.f[scen,tt])                     ### potential breeders that remain at sea    
    
      ### THE TOTAL TRAL POPULATION ###
      Ntot.f[scen,tt]<-sum(IM.f[scen,tt,,3])+Ntot.breed.f[scen,tt]+N.atsea.f[scen,tt]  ## total population size is all the immatures plus adult breeders and adults at sea
    
    
    } ### end future loop
    
    ## CALCULATE ANNUAL POP GROWTH RATE ##
    for (fut2 in 1:(FUT.YEAR-1)){
      fut.lambda[scen,fut2] <- Ntot.f[scen,fut2+1]/max(1,Ntot.f[scen,fut2])                                 ### inserted safety to prevent denominator being 0
    } # fut2
    
    
    ## DERIVED MEAN FUTURE GROWTH RATE 
    fut.growth.rate[scen] <- exp((1/(FUT.YEAR-1))*sum(log(fut.lambda[scen,1:(FUT.YEAR-1)])))   # Geometric mean
    
  } # end future projection scenarios
    
}  ## end model loop

",fill = TRUE)
sink()





#########################################################################
# PREPARE DATA FOR MODEL
#########################################################################

# Bundle data
jags.data <- list(marr.j = chick.marray,
                  marr.a = adult.marray,
                  n.occasions = dim(chick.marray)[2],
                  r.j=apply(chick.marray,1,sum),
                  r.a=apply(adult.marray,1,sum),
                  goodyear=goodyears$p.sel,
                  #goodyear=goodyears$prop.seen,   ### if using a continuous effort correction
                  juv.poss=phi.juv.possible$JuvSurv, ### sets the annual survival of juveniles to the mean if <70 were ringed
                  
                  ### count data
                  n.sites=n.sites,
                  T = n.years,
                  prop.sites=mean.props,
                  y.count=R,    ### use log(R) here if using the logscale model
                  
                  ### breeding success data
                  J=PROD.DAT$J,
                  R=PROD.DAT$R,
                  
                  ### longline effort data
                  #longline=longlineICCAT,
                  
                  # ### FUTURE PROJECTION
                  FUT.YEAR=30,  ### for different scenarios future starts at 1
                  n.scenarios=3,
                  fut.surv.change=as.matrix(fut.surv.change[,2:4]),  ## future survival rate change - matrix that adjusts gradual decrease in survival
                  fut.fec.change=c(1,2,1)     ## future fecundity change - vector with one element for each scenario
                  )


# Initial values 
inits <- function(){list(mean.phi.ad = runif(1, 0.7, 0.97),
                         mean.phi.juv = runif(1, 0.5, 0.9),
                         mean.p.ad = runif(3, 0.2, 1),
                         mean.p.juv = runif(3, 0, 1),
                         Ntot.breed= c(runif(1, 4950, 5050),rep(NA,n.years-1)),
                         JUV= c(rnorm(1, 246, 0.1),rep(NA,n.years-1)),
                         N.atsea= c(rnorm(1, 530, 0.1),rep(NA,n.years-1)),
                         # IM[,1,1]= c(rnorm(1, 324, 0.1),rep(NA,n.years-1)),
                         # IM[,2,1]= c(rnorm(1, 257, 0.1),rep(NA,n.years-1)),
                         # IM[,3,1]= c(rnorm(1, 462, 0.1),rep(NA,n.years-1)),
                         # IM[,4,1]= c(rnorm(1, 207, 0.1),rep(NA,n.years-1)),
                         # IM[,5,1]= c(rnorm(1, 700, 0.1),rep(NA,n.years-1)),
                         # IM[,6,1]= c(runif(1, 150, 300),rep(NA,n.years-1)),
                         sigma.obs=matrix(runif(n.sites*n.years,1,20),ncol=n.years))}

 

# Parameters monitored
parameters <- c("mean.phi.ad","mean.phi.juv","mean.fec",
                "pop.growth.rate","fut.growth.rate","lambda",
                "agebeta","Ntot","Ntot.f","phi.ad","phi.juv","Ntot.breed",   ## added Ntot.breed to provide better contrast with Ntot?
                #new
                "ann.fec", "mean.p.juv","mean.p.ad",#"p.juv","p.ad",
                "IM") ## added IM and JUV to facilitate LTRE analysis

### REDUCE WORKSPACE SIZE
rm(list=setdiff(ls(), c("parameters","jags.data","inits","n.years","n.sites")))
gc()



# RUN THE MODEL {took 3 hours for niter=125000)
## _logscale model requires log(R) as input for count data
## THIS DOES NOT CONVERGE
# TRALipm <- autojags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM\\TRAL_IPM_marray_age_recruit_immat_FINAL.jags",
#                     n.chains = nc, n.thin = nt, n.burnin = nb,parallel=T, #n.iter = ni)
#                     Rhat.limit=1.2, max.iter=200000)  

nt <- 20
nb <- 25000
nad <- 10000
nc <- 3
ns <- 10000 #longest

TRALipm <- run.jags(data=jags.data, inits=inits, parameters, 
                    model="C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM\\TRAL_IPM_FINAL_REV2022.jags",
                    n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                    method = "rjparallel") 

str(TRALipm)



#########################################################################
# SAVE OUTPUT - RESULT PROCESSING in TRAL_IPM_result_summaries.r
#########################################################################
### DO NOT UPLOAD THIS TO GITHUB - IT WILL CORRUPT THE REPOSITORY
### monitoring IM and JUV makes this workspace way too large to save
### writing output files one by one

## updated script for 'runjags' output

for(p in parameters){
summary_tralipm <- summary(TRALipm, vars=p)   ### ("mean.phi.ad","mean.phi.juv","mean.fec","mean.propensity", "mean.recruit",)does not run for all variables if IM is included
predictions <- data.frame(summary_tralipm,
                          parameter = row.names(summary_tralipm))
#summary_tralipm_df <- as.data.frame(summary_tralipm)
fwrite(predictions, sprintf("IPM_output_%s.csv",p))
}



#### FOR IM the above does not work

### need the following parameters from model
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="IM[18,30,1]") # IM: 2210-3829
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="IM[1,1,1]") 
IM <- as.matrix(TRALipm$mcmc[[1]][,c(406:945)])    #### modify sequence of dimnames to only include ,,1] 

### cannot include all chains as memory allocation error
#for(ch in 2:nc){
#   IM<-rbind(IM,as.matrix(TRALipm$mcmc[[ch]][,c(406:2025)]))  ## only for the years coinciding with the fecundity and pop size
# }

IMpredictions<-as_tibble(IM) %>% gather(key="parameter", value="value") %>%
	group_by(parameter) %>%
	summarise(mean=mean(value), median=median(value), lcl=quantile(value,0.025), ucl=quantile(value,0.975))
fwrite(IMpredictions, sprintf("IPM_output_%s.csv",p))







#### old code that worked before inclusion of IM and JUV - not possible with memory limitations

View(summary_tralipm_df)
head(summary_tralipm_df)
min(summary_tralipm_df$SSeff) #Ntot[1]
max(summary_tralipm_df$psrf) #Ntot[1]

addsummary_tralipm <- add.summary(TRALipm,plots = runjags.getOption("predraw.plots"))
addsummary_tralipm #18 min

plot(addsummary_tralipm, layout=c(2,2))

predictions <- data.frame(summary(addsummary_tralipm),
                          parameter = row.names(summary(addsummary_tralipm)))
head(predictions)
row.names(predictions) <- 1:nrow(predictions)

predictions <- predictions[1:253,]   ### 200 cuts off ann.fec
#predictions[1:5,]
tail(predictions)

predictions$Mode <- NULL
np <- names(predictions) 
names(predictions) <- c("lcl",np[2],"ucl",np[4:9],"Rhat",np[11])

max(predictions$Rhat)



setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
save.image("TRAL_IPM_output_FINAL_REV2022b.RData")


