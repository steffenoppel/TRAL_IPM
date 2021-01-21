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

library(tidyverse)
library(lubridate)
library(data.table)
library(jagsUI)
filter<-dplyr::filter
select<-dplyr::select



#########################################################################
# LOAD PRE-PREPARED DATA ON COUNTS AND BREEDING SUCCESS
#########################################################################
### see 'IPM_DATA_PREPARATION.R' for details on how data are aggregated

### NOTE THAT SORT ORDER OF GONYDALE AND GP VALLEY HAS SHIFTED ON 15 Jan 2021 (due to switch to if_else on R4.0.2)

## LOAD PREPARED M-ARRAY FOR SURVIVAL ESTIMATION
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
load("TRAL_IPM_input.marray.RData")

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
  filter(Year>2003) %>%
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
  filter(Year>2003) %>%
  mutate(Site=if_else(Site %in% c('GP Valley','West Point'),'GP Valley',as.character(Site))) %>%
  mutate(Site=if_else(Site %in% c('Gonydale','Green Hill','Hummocks'),'Gonydale',as.character(Site))) %>%
  group_by(Year,Site) %>%
  summarise(Count=sum(Count, na.rm=T)) %>%
  mutate(Count=ifelse(Count==0,NA,Count)) %>%
  spread(key=Site, value=Count)

### NOTE THAT USE OF 'if_else' switches Gonydale and GP_Valley columns around [only happens in R4.0.2!!]
TRAL.chick[8,4]<-CHICKCOUNT[11,4] ### in 2011 no adults were counted in Green hill and Hummocks, so we cannot add up the chicks across those 3 sites


TRAL.chick$tot<-rowSums(TRAL.chick[,2:9], na.rm=T)
#TRAL.chick$tot[c(3,13,19)]<-NA   ## when start in 2001
TRAL.chick$tot[10]<-NA   ## when start in 2004

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
# IPM runs from 2001-2020
# survival analysis runs from 1979-2021
names(TRAL_CHICK)
TRAL.pop$Year

OFFSET<-min(which(!is.na(match(as.numeric(names(TRAL_CHICK)[2:44]),TRAL.pop$Year))))
names(TRAL_CHICK)[OFFSET+1]
TRAL.pop$Year[1]

#########################################################################
# SPECIFY MODEL IN JAGS
#########################################################################
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
sink("TRAL_IPM_marray_simplified_v4.jags")
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
    # -------------------------------------------------
    
    #-------------------------------------------------  
    # 1. PRIORS FOR ALL DATA SETS
    #-------------------------------------------------
    
    
    # -------------------------------------------------        
    # 1.1. Priors and constraints FOR FECUNDITY
    # -------------------------------------------------
    
    for (t in 1:T){  
      ann.fec[t] ~ dunif(0.05,0.65) ## dnorm(0.32,10) T(0.001,0.999)        ## Informative Priors on fecundity based on Wanless et al 2009
    } #t

    
    # -------------------------------------------------        
    # 1.2. Priors and constraints FOR POPULATION COUNTS
    # -------------------------------------------------
    for (s in 1:n.sites){			### start loop over every study area
      for (t in 1:T){			### start loop over every year
        sigma.obs[s,t] ~ dunif(0,20)	#Prior for SD of observation process (variation in detectability)
        tau.obs[s,t]<-pow(sigma.obs[s,t],-2)
      }
    }


    # -------------------------------------------------        
    # 1.3. Priors and constraints FOR SURVIVAL
    # -------------------------------------------------
    
    ### RECAPTURE PROBABILITY
    mean.p.juv ~ dunif(0.2, 1)	           # Prior for mean juvenile recapture - should be higher than 20% if they survive!
    mean.p.ad ~ dunif(0.2, 1)	           # Prior for mean adult recapture - should be higher than 20%
    mu.p.juv <- log(mean.p.juv / (1-mean.p.juv)) # Logit transformation
    mu.p.ad <- log(mean.p.ad / (1-mean.p.ad)) # Logit transformation

    ## PRIORS FOR RANDOM EFFECTS
    sigma.p ~ dunif(0, 1)                # Prior for standard deviation
    tau.p <- pow(sigma.p, -2)

    
    ### SURVIVAL PROBABILITY
    mean.phi.juv ~ dunif(0.5, 0.9)             # Prior for mean juvenile survival first year 0.757, second year 0.973 in Laysan albatross
    mean.phi.ad ~ dunif(0.7, 1)             # Prior for mean adult survival - should be higher than 70%
    mu.juv <- log(mean.phi.juv / (1-mean.phi.juv)) # Logit transformation
    mu.ad <- log(mean.phi.ad / (1-mean.phi.ad)) # Logit transformation

    ## PRIORS FOR RANDOM EFFECTS
    sigma.phi ~ dunif(0, 1)                # Prior for standard deviation
    tau.phi <- pow(sigma.phi, -2)
    
    ## RANDOM TIME EFFECT ON SURVIVAL ONLY FOR ADULTS (age group = 2)
    for (t in 1:(n.occasions-1)){
      logit(phi.juv[t]) <- mu.juv + eps.phi[t]
      logit(phi.ad[t]) <- mu.ad + eps.phi[t]
      eps.phi[t] ~ dnorm(0, tau.phi)
    }
    for (t in 1:(n.occasions)){
      logit(p.ad[t])  <- mu.p.ad + eps.p[t]
      logit(p.juv[t])  <- mu.p.juv + eps.p[t]
      eps.p[t] ~ dnorm(0, tau.p) 
    }


    #-------------------------------------------------  
    # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
    #-------------------------------------------------
    
    # -------------------------------------------------        
    # 2.1. System process: female based matrix model
    # -------------------------------------------------
    
    for (tt in 2:T){
    
      ## THE PRE-BREEDING YEARS ##
    
      nestlings[tt] <- ann.fec[tt] * 0.5 * Ntot.breed[tt]                                                     ### number of locally produced FEMALE chicks
      JUV[tt] ~ dpois(nestlings[tt])                                                                     ### need a discrete number otherwise dbin will fail, dpois must be >0
      N1[tt]  ~ dbin(mean.phi.juv, round(JUV[tt-1]))                                                    ### number of 1-year old survivors 
      N2[tt] ~ dbin(mean.phi.juv, round(N1[tt-1]))                                                      ### number of 2-year old survivors
      N3[tt] ~ dbin(mean.phi.juv, round(N2[tt-1]))                                                       ### number of 3-year old survivors
      N4[tt] ~ dbin(mean.phi.juv, round(N3[tt-1]))                                                       ### number of 4-year old survivors
      N.notrecruited.surv[tt] ~ dbin(phi.ad[tt+25], round(N.notrecruited[tt-1]))                        ### number of older immature survivors that have not recruited 
      N.recruits[tt] ~ dbin(p.juv[tt+26], round(N.notrecruited.surv[tt]+N4[tt]))                             ### number of this years recruiters
      N.notrecruited[tt] <-round(N.notrecruited.surv[tt]+N4[tt]-N.recruits[tt])                         ### number of birds not yet recruited at the end of this year
 
      ## THE BREEDING POPULATION ##
      # Ntot.breed comprised of first-time breeders, previous skippers, and previous unsuccessful breeders
      # simplified in simplified_v2 to just adult survivors with p.ad as proportion returning
      #N.prev.succ[tt] ~ dbin(ann.fec[tt-1], round(Ntot.breed[tt-1]))                    ### number of successful breeders from previous year
      #N.prev.succ.atsea[tt] ~ dbin(phi.ad[tt+25], N.prev.succ[tt])                     ### previous year's successful breeders at sea
      #N.prev.unsucc.breed[tt] <- round(Ntot.breed[tt-1]-N.prev.succ[tt])
      #N.prev.unsucc.return[tt] ~ dbin(skip.prob[tt], N.prev.unsucc.breed[tt])             ### previous year's unsuccessful breeders that are returning to breed
      #N.prev.unsucc.atsea[tt] ~ dbin(phi.ad[tt+25], (N.prev.unsucc.breed[tt]-N.prev.unsucc.return[tt])) ### previous year's unsuccessful breeders that are staying at sea
      #N.breed.ready[tt]<- round(N.prev.unsucc.return[tt]+ N.prev.unsucc.atsea[tt-1]+N.prev.succ.atsea[tt-1]+N.recruits[tt-1])       ### number of available breeders is failed breeders from previous year plus failed and successful breeders from 2 years ago plus recruits (N9)
      #Ntot.breed[tt] ~ dbin(phi.ad[tt+25], N.breed.ready[tt])                         ### number of potential old breeders is the number of survivors from previous year breeders and nonbreeders
      
      N.ad.surv[tt] ~ dbin(phi.ad[tt+25], round(Ntot.breed[tt-1]+N.atsea[tt-1]))           ### previous year's adults that survive
      N.breed.ready[tt] ~ dbin(p.ad[tt+26], N.ad.surv[tt])                  ### number of available breeders is proportion of survivors that returns
      Ntot.breed[tt]<- round(N.breed.ready[tt]+N.recruits[tt])              ### number of counted breeders is sum of old breeders returning and first recruits
      N.atsea[tt] <- round(N.ad.surv[tt]-N.breed.ready[tt])                     ### potential breeders that remain at sea    
      
      ### THE TOTAL TRAL POPULATION ###
      Ntot[tt]<-N1[tt]+N2[tt]+N3[tt]+N.notrecruited[tt]+Ntot.breed[tt]+N.atsea[tt]  ## total population size is all the immatures plus adult breeders and adults at sea

    } # tt
    
    
    
    ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1 - based on deterministic multiplications
    ## ADJUSTED BASED ON PAST POPULATION SIZES
    Ntot.breed[1] ~ dunif(1500,2000) ###   ### sum of counts is 2400, but we take average over 3 years because 2001 was an outlier year
    JUV[1]<-round(Ntot.breed[1]*0.5*ann.fec[1])
    N1[1] ~ dunif(200,700)  ##<-round(JUV[1]*mean.phi.juv)
    N2[1] ~ dunif(100,600)  ##<-round(N1[1]*mean.phi.juv)
    N3[1] ~ dunif(300,800)  ## very large number of birds due to 2400 breeders in 2001 <-round(N2[1]*mean.phi.juv)
    N4[1] ~ dunif(100,400)  ##<-round(N3[1]*mean.phi.juv)
    N.notrecruited[1] ~ dunif(50,800)  ### all birds not recruited, including those from 1999 bumper year <-round(N4[1]*mean.phi.juv*0.9)                         ### number of birds not yet recruited ~90% at age 4
    N.atsea[1] ~ dunif(100,1000)
    Ntot[1]<-N1[1]+N2[1]+N3[1]+N.notrecruited[1]+Ntot.breed[1]+N.atsea[1]  ## total population size is all the immatures plus adult breeders and adults at sea


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
    for (t in 1:(T)){
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
      q.juv[t] <- 1-p.juv[t]            # Probability of non-recapture
      q.ad[t] <- 1-p.ad[t]            # Probability of non-recapture
      pr.j[t,t] <- phi.juv[t]*p.juv[t]
      pr.a[t,t] <- phi.ad[t]*p.ad[t]
    
      # Above main diagonal
      for (j in (t+1):(n.occasions-1)){
        pr.j[t,j] <- phi.juv[t]*prod(phi.ad[(t+1):j])*prod(q.juv[t:(j-1)])*p.juv[j]
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


for(scen in 1:n.scenarios){


    ### COPY POPULATIONS FROM LAST YEAR OF DATA SERIES

    nestlings.f[scen,1] <- round(fut.fec.change[scen]*mean.fec* 0.5 * Ntot.breed[T])                                             ### number of locally produced FEMALE chicks based on average fecundity - to use just one take ann.fec[FUT.int[tt]] 
    N1.f[scen,1]  ~ dbin(mean.phi.juv, max(1,round(nestlings[T])))                                                    ### number of 1-year old survivors 
    N2.f[scen,1] ~ dbin(mean.phi.juv, max(1,round(N1[T])))                                                      ### number of 2-year old survivors
    N3.f[scen,1] ~ dbin(mean.phi.juv, max(1,round(N2[T])))                                                       ### number of 3-year old survivors
    N4.f[scen,1] ~ dbin(mean.phi.juv, max(1,round(N3[T])))                                                       ### number of 4-year old survivors
    N.notrecruited.surv.f[scen,1] ~ dbin(fut.surv.change[scen]*mean.phi.ad, round(N.notrecruited[T]))                        ### number of older immature survivors that have not recruited 
    N.recruits.f[scen,1] ~ dbin(mean.p.juv, round(N.notrecruited.surv.f[scen,1]+N4.f[scen,1]))                             ### number of this years recruiters
    N.notrecruited.f[scen,1] <-round(N.notrecruited.surv.f[scen,1]+N4.f[scen,1]-N.recruits.f[scen,1])                         ### number of birds not yet recruited at the end of this year

    N.ad.surv.f[scen,1] ~ dbin(fut.surv.change[scen]*mean.phi.ad, round(Ntot.breed[T]+N.atsea[T]))           ### previous year's adults that survive
    N.breed.ready.f[scen,1] ~ dbin(p.ad[T+26], round(N.ad.surv.f[scen,1]))                  ### number of available breeders is proportion of survivors that returns, with fecundity INCLUDED in return probability
    Ntot.breed.f[scen,1]<- round(N.breed.ready.f[scen,1]+N.recruits.f[scen,1])              ### number of counted breeders is sum of old breeders returning and first recruits
    N.atsea.f[scen,1] <- round(N.ad.surv.f[scen,1]-N.breed.ready.f[scen,1])                     ### potential breeders that remain at sea
    N.succ.breed.f[scen,1] ~ dbin(fut.fec.change[scen]*mean.fec, round(Ntot.breed.f[scen,1]))                  ### these birds will  remain at sea because tey bred successfully
      
    ### THE TOTAL TRAL POPULATION ###
    Ntot.f[scen,1]<-N1.f[scen,1]+N2.f[scen,1]+N3.f[scen,1]+N.notrecruited.f[scen,1]+Ntot.breed.f[scen,1]+N.atsea.f[scen,1]  ## total population size is all the immatures plus adult breeders and adults at sea

   
    for (tt in 2:FUT.YEAR){
    

    # -------------------------------------------------        
    # 4.1. System process for future
    # -------------------------------------------------

    ## INCLUDE CARRYING CAPACITY OF 2500 breeding pairs (slightly more than maximum ever counted)
    carr.capacity[scen,tt] ~ dnorm(2500,0.05)
    
    ## THE PRE-BREEDING YEARS ##
    ## because it goes for 30 years, all pops must be safeguarded to not become 0 because that leads to invald parent error
    
    nestlings.f[scen,tt] <- round(fut.fec.change[scen]*mean.fec* 0.5 * Ntot.breed.f[scen,tt])                       ### number of locally produced FEMALE chicks based on average fecundity - to use just one take ann.fec[FUT.int[tt]] 
    N1.f[scen,tt]  ~ dbin(mean.phi.juv, max(1,round(nestlings.f[scen,tt-1])))                                              ### number of 1-year old survivors 
    N2.f[scen,tt] ~ dbin(mean.phi.juv, max(1,round(N1.f[scen,tt-1])))                                                      ### number of 2-year old survivors
    N3.f[scen,tt] ~ dbin(mean.phi.juv, max(1,round(N2.f[scen,tt-1])))                                                       ### number of 3-year old survivors
    N4.f[scen,tt] ~ dbin(mean.phi.juv, max(1,round(N3.f[scen,tt-1])))                                                       ### number of 4-year old survivors
    N.notrecruited.surv.f[scen,tt] ~ dbin(fut.surv.change[scen]*mean.phi.ad, round(N.notrecruited.f[scen,tt-1]))                        ### number of older immature survivors that have not recruited 
    N.recruits.f[scen,tt] ~ dbin(mean.p.juv, round(N.notrecruited.surv.f[scen,tt]+N4.f[scen,tt]))                             ### number of this years recruiters
    N.notrecruited.f[scen,tt] <-round(N.notrecruited.surv.f[scen,tt]+N4.f[scen,tt]-N.recruits.f[scen,tt])                         ### number of birds not yet recruited at the end of this year

    ## THE BREEDING POPULATION ##
      N.ad.surv.f[scen,tt] ~ dbin(fut.surv.change[scen]*mean.phi.ad, round((Ntot.breed.f[scen,tt-1]-N.succ.breed.f[scen,tt-1])+N.atsea.f[scen,tt-1]))           ### previous year's adults that survive
      N.prev.succ.f[scen,tt] ~ dbin(fut.surv.change[scen]*mean.phi.ad, round(N.succ.breed.f[scen,tt-1]))                  ### these birds will  remain at sea because tey bred successfully
      N.breed.ready.f[scen,tt] ~ dbin(min(0.99,(mean.p.ad/(1-mean.fec))), max(1,round(N.ad.surv.f[scen,tt])))                  ### number of available breeders is proportion of survivors that returns, with fecundity partialled out of return probability
      Ntot.breed.f[scen,tt]<- min(carr.capacity[scen,tt],round(N.breed.ready.f[scen,tt]+N.recruits.f[scen,tt]))              ### number of counted breeders is sum of old breeders returning and first recruits
      N.succ.breed.f[scen,tt] ~ dbin(fut.fec.change[scen]*mean.fec, round(Ntot.breed.f[scen,tt]))                  ### these birds will  remain at sea because tey bred successfully
      N.atsea.f[scen,tt] <- round(N.ad.surv.f[scen,tt]-N.breed.ready.f[scen,tt]+N.prev.succ.f[scen,tt])                     ### potential breeders that remain at sea    

    ### THE TOTAL TRAL POPULATION ###
    Ntot.f[scen,tt]<-N1.f[scen,tt]+N2.f[scen,tt]+N3.f[scen,tt]+N.notrecruited.f[scen,tt]+Ntot.breed.f[scen,tt]+N.atsea.f[scen,tt]  ## total population size is all the immatures plus adult breeders and adults at sea


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

### FIX FIRST POPULATION COUNT FOR 2001 (too high due to exceptional breeding success in 1999)
# R[1,]<-round(colMeans(R[4:6,], na.rm=T),0)
# sum(R[1,])   ## reduces count from 2400 to 1638

# Bundle data
jags.data <- list(marr.j = chick.marray,
                  marr.a = adult.marray,
                  n.occasions = dim(chick.marray)[2],
                  r.j=apply(chick.marray,1,sum),
                  r.a=apply(adult.marray,1,sum),
                  
                  ### count data
                  n.sites=n.sites,
                  T = n.years,
                  prop.sites=mean.props,
                  y.count=R,
                  
                  ### breeding success data
                  J=PROD.DAT$J,
                  R=PROD.DAT$R,
                  
                  ### longline effort data
                  #longline=longlineICCAT,
                  
                  # ### FUTURE PROJECTION
                  FUT.YEAR=30,  ### for different scenarios future starts at 1
                  n.scenarios=3,
                  fut.surv.change=c(1,1,0.9),     ## future survival rate change - vector with one element for each scenario
                  fut.fec.change=c(1,2,1)     ## future fecundity change - vector with one element for each scenario
                  )


# Initial values 
inits <- function(){list(mean.phi.ad = runif(1, 0.7, 1),
                         mean.phi.juv = runif(1, 0.5, 0.9),
                         mean.p.ad = runif(1, 0.2, 1),
                         mean.p.juv = runif(1, 0.2, 1),
                         Ntot.breed= c(runif(1, 1500, 2000),rep(NA,n.years-1)),
                         #bycatch = rnorm(1,0,0.01),
                         #hookpod = rnorm(1,0,0.01),
                         
                         ### count data
                         #sigma.proc=runif(n.sites,0,10),
                         #mean.lambda=runif(n.sites,0.1,2),
                         #N.est=N.init,  ## results in 'inconsistent with parents' error
                         sigma.obs=matrix(runif(n.sites*n.years,1,20),ncol=n.years))}

 

# Parameters monitored
parameters <- c("ann.fec","mean.phi.ad","mean.phi.juv","mean.fec","mean.skip","mean.p.ad","mean.p.juv","pop.growth.rate","fut.growth.rate","Ntot","Ntot.f","lambda","p.ad")  

# MCMC settings
ni <- 250000
nt <- 4
nb <- 50000
nc <- 3



# RUN THE FOUR SCENARIOS {took 2 hours for niter=150000)
TRALipm <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM\\TRAL_IPM_marray_simplified_v4.jags",
                    n.chains = nc, n.thin = nt, n.burnin = nb,parallel=T, n.iter = ni)
                    #Rhat.limit=1.5, max.iter=150000)  





#########################################################################
# SAVE OUTPUT - RESULT PROCESSING in TRAL_IPM_result_summaries.r
#########################################################################
### DO NOT UPLOAD THIS TO GITHUB - IT WILL CORRUPT THE REPOSITORY

setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
save.image("TRAL_IPM_output_v5_Ntot.RData")











