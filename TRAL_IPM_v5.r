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
  mutate(Site=ifelse(Site %in% c('GP Valley','West Point'),'GP Valley',Site)) %>%
  mutate(Site=ifelse(Site %in% c('Gonydale','Green Hill','Hummocks'),'Gonydale',Site)) %>%
  group_by(Year,Site) %>%
  summarise(Count=sum(Count, na.rm=T)) %>%
  mutate(Count=ifelse(Count==0,NA,Count)) %>%
  spread(key=Site, value=Count)

## to account for data gaps, need to quantify what proportion of the population was counted

TRAL.props<-prop.table(as.matrix(TRAL.pop[,2:9]),1)
mean.props<-apply(TRAL.props[c(1,4,6:10,12:16,18:19),],2,mean)

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



#### BREEDING SUCCESS DATA FOR FECUNDITY ######
TRAL.chick<-CHICKCOUNT  ##fread("TRAL_CHIC_counts_2001_2020.csv")
## detailed nest data from Gonydale only available for some years - we replace those years in the chick count matrix
TRAL.bs<-FECUND ##fread("TRAL_breed_success_2006_2019.csv")
TRAL.chick[match(TRAL.bs$Year,TRAL.chick$Year),4]<-as.integer(R[match(TRAL.bs$Year,TRAL.chick$Year),3]*TRAL.bs$BREED_SUCC)
TRAL.chick[14,4]<-as.integer(TRAL.bs$n_nests[7]*TRAL.bs$BREED_SUCC[7])  ## use Gonydale nest monitoring to fill in data gaps
TRAL.chick<-TRAL.chick %>% gather(key='Site', value='Count',-Year) %>%
  mutate(Site=ifelse(Site %in% c('GP Valley','West Point'),'GP Valley',Site)) %>%
  mutate(Site=ifelse(Site %in% c('Gonydale','Green Hill','Hummocks'),'Gonydale',Site)) %>%
  group_by(Year,Site) %>%
  summarise(Count=sum(Count, na.rm=T)) %>%
  mutate(Count=ifelse(Count==0,NA,Count)) %>%
  spread(key=Site, value=Count)

TRAL.chick$tot<-rowSums(TRAL.chick[,2:9], na.rm=T)
TRAL.chick$tot[c(3,13,19)]<-NA
#J<-TRAL.chick$tot
J<- as.matrix(TRAL.chick[,2:9])

### specify constants for JAGS
n.years<-dim(R)[1]		## defines the number of years
n.sites<-dim(R)[2]    ## defines the number of study areas



### DIMENSION MISMATCH IN DATA
# IPM runs from 2001-2020
# survival analysis runs from 1979-2021
min(which(!is.na(match(as.numeric(names(TRAL_CHICK)[2:44]),POPSIZE$Year))))



#########################################################################
# SPECIFY MODEL IN JAGS
#########################################################################
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
sink("TRAL_IPM_fut_marray_v1.jags")
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
    # -------------------------------------------------
    
    #-------------------------------------------------  
    # 1. PRIORS FOR ALL DATA SETS
    #-------------------------------------------------
    
    
    # -------------------------------------------------        
    # 1.1. Priors and constraints FOR FECUNDITY
    # -------------------------------------------------
    
    for (t in 1:T){  
      ann.fec[t] ~ dunif(0.05,0.65) ## dnorm(0.32,10) T(0.001,0.999)        ## Informative Priors on fecundity based on Wanless et al 2009
      skip.prob[t] ~ dunif(0.05,0.95) ## dunif(0.22,0.32) T(0.001,0.999)           ## PRIOR FOR ADULT FAILED BREEDER SKIPPING PROBABILITY from Wanless et al 2009
    } #t

    
    # -------------------------------------------------        
    # 1.2. Priors and constraints FOR POPULATION COUNTS
    # -------------------------------------------------
    for (s in 1:n.sites){			### start loop over every study area
      #N.est[1,s] ~ dunif(0,400)   ## draw random value from a uniform distribution between 0 and 200 for initial population size
      sigma.proc[s] ~ dunif(0,10)	#Prior for SD of state process (annual variation in pop size)
      tau.proc[s]<-pow(sigma.proc[s],-2)
      sigma.obs[s] ~ dunif(0,100)	#Prior for SD of observation process (variation in detectability)
      tau.obs[s]<-pow(sigma.obs[s],-2)
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
    mean.phi.juv ~ dunif(0, 1)             # Prior for mean juvenile survival
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
      N5[tt] ~ dbin(phi.ad[tt+22], round(N4[tt-1]))                                                       ### number of 5-year old survivors
      N6[tt] ~ dbin(phi.ad[tt+22], round(N5[tt-1]))                                                       ### number of 6-year old survivors
      N7[tt] ~ dbin(phi.ad[tt+22], round(N6[tt-1]))                                                       ### number of 7-year old survivors
      N8[tt] ~ dbin(phi.ad[tt+22], round(N7[tt-1]))                                                       ### number of 8-year old survivors
      N9[tt] ~ dbin(phi.ad[tt+22], round(N8[tt-1]))                                                       ### number of 9-year old survivors
    

      ## THE BREEDING POPULATION ##

      # Ntot.breed comprised of first-time breeders, previous skippers, and previous unsuccessful breeders
      N.succ.breed[tt] ~ dbin(ann.fec[tt], round(Ntot.breed[tt]))                     ### number of successful breeders not breeding the following year
      N.unsucc.breed[tt] <- round(Ntot.breed[tt]-N.succ.breed[tt])
      N.fail.skip[tt] ~ dbin(skip.prob[tt], max(1,N.unsucc.breed[tt]))   ### number of unsuccessful breeders not breeding the following year
      N.nonbreed[tt] <- round(sum(N.succ.breed[tt-1],N.fail.skip[tt-1])) ### number of nonbreeders as the sum of previous year successful and skipping unsuccessful breeders
      N.non.breed[tt] ~ dbin(phi.ad[tt+22], N.nonbreed[tt])  ### number of non-breeders in a given year is composed of the successful skippers and failed skippers from previous year       

      N.breed.ready[tt]<- round((Ntot.breed[tt-1]-N.succ.breed[tt-1]-N.fail.skip[tt-1])+N9[tt-1]+N.non.breed[tt-1])       ### number of available breeders is failed breeders from previous year plus non-breeders from previous year plus recruits (N9)
      Ntot.breed[tt] ~ dbin(phi.ad[tt+22], N.breed.ready[tt])   ### number of potential old breeders is the number of survivors from previous year breeders and nonbreeders

 
    } # tt
    
    
    
    ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1 - based on deterministic multiplications
    Ntot.breed[1] ~ dnorm(1638,1000) #1638   ##~ dunif(1000,2000) ###   ### sum of counts is 2400, but we take average over 3 years because 2001 was an outlier year
    JUV[1]<-round(Ntot.breed[1]*0.5*ann.fec[1])
    N1[1]<-round(JUV[1]*mean.phi.juv)
    N2[1]<-round(N1[1]*mean.phi.juv)
    N3[1]<-round(N2[1]*mean.phi.juv)
    N4[1]<-round(N3[1]*mean.phi.juv)
    N5[1]<-round(N4[1]*mean.phi.ad)
    N6[1]<-round(N5[1]*mean.phi.ad)
    N7[1]<-round(N6[1]*mean.phi.ad)
    N8[1]<-round(N7[1]*mean.phi.ad)
    N9[1]<-round(N8[1]*mean.phi.ad)
    N.nonbreed[1]<-round(sum(N.succ.breed[1],N.fail.skip[1]))
    N.non.breed[1]~ dbin(mean.phi.ad, N.nonbreed[1])
    N.succ.breed[1]<- round(Ntot.breed[1]*ann.fec[1])
    N.unsucc.breed[1] <- (Ntot.breed[1]-N.succ.breed[1])
    N.fail.skip[1]<- round(N.unsucc.breed[1]*skip.prob[1])


    # -------------------------------------------------        
    # 2.2. Observation process for population counts: state-space model of annual counts
    # -------------------------------------------------
    
    for (s in 1:n.sites){			### start loop over every study area
    
      ## Observation process
    
      for (t in 1:T){
        y.count[t,s] ~ dnorm(Ntot.breed[t]*prop.sites[s], tau.obs[s])								# Distribution for random error in observed numbers (counts)
      }														# run this loop over t= nyears
    }		## end site loop

    
    # -------------------------------------------------        
    # 2.3. Likelihood for fecundity: Logistic regression from the number of surveyed broods
    # -------------------------------------------------
    for (s in 1:n.sites){			### start loop over every study area
      for (t in 1:(T-1)){
        J[t,s] ~ dpois(rho.fec[t,s])
        rho.fec[t,s] <- y.count[t,s]*ann.fec[t]
      } #	close loop over every year in which we have fecundity data
    }		## end site loop
    
    
    
    
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
    

    
    #-------------------------------------------------  
    # 3. DERIVED PARAMETERS FOR OUTPUT REPORTING
    #-------------------------------------------------
    

    ## DERIVED POPULATION GROWTH RATE PER YEAR
    for (t in 1:(T-1)){
      lambda[t]<-Ntot.breed[t+1]/max(1,Ntot.breed[t])  ## division by 0 creates invalid parent value
    }		## end year loop
    

    ## DERIVED MEAN FECUNDITY 
    mean.fec <- mean(ann.fec)
    mean.skip <- mean(skip.prob)
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
    N5.f[scen,1] ~ dbin(fut.surv.change[scen]*mean.phi.ad, max(1,round(N4[T])))                                                       ### number of 5-year old survivors
    N6.f[scen,1] ~ dbin(fut.surv.change[scen]*mean.phi.ad, max(1,round(N5[T])))                                                       ### number of 6-year old survivors
    N7.f[scen,1] ~ dbin(fut.surv.change[scen]*mean.phi.ad, max(1,round(N6[T])))                                                       ### number of 7-year old survivors
    N8.f[scen,1] ~ dbin(fut.surv.change[scen]*mean.phi.ad, max(1,round(N7[T])))                                                       ### number of 8-year old survivors
    N9.f[scen,1] ~ dbin(fut.surv.change[scen]*mean.phi.ad, max(1,round(N8[T])))                                                       ### number of 9-year old survivors
    N.succ.breed.f[scen,1] ~ dbin(fut.fec.change[scen]*mean.fec, max(1,round(Ntot.breed[T])))                     ### number of successful breeders not breeding the following year
    N.unsucc.breed.f[scen,1] <- round(Ntot.breed[T]-N.succ.breed[T])
    N.fail.skip.f[scen,1] ~ dbin(mean.skip, max(1,N.unsucc.breed[T]))   ### number of unsuccessful breeders not breeding the following year
    N.nonbreed.f[scen,1] <- round(sum(N.succ.breed[T-1],N.fail.skip[T-1])) ### number of nonbreeders as the sum of previous year successful and skipping unsuccessful breeders
    N.non.breed.f[scen,1] ~ dbin(fut.surv.change[scen]*mean.phi.ad, max(1,N.nonbreed[T]))  ### number of non-breeders in a given year is composed of the successful skippers and failed skippers from previous year       
    N.breed.ready.f[scen,1]<- round((Ntot.breed[T-1]-N.succ.breed[T-1]-N.fail.skip[T-1])+N9[T-1]+N.non.breed[T-1])       ### number of available breeders is failed breeders from previous year plus non-breeders from previous year plus recruits (N9)
    Ntot.breed.f[scen,1] ~ dbin(fut.surv.change[scen]*mean.phi.ad, max(1,N.breed.ready[T]))   ### number of potential old breeders is the number of survivors from previous year breeders and nonbreeders
    carr.capacity[scen,1]<-2500
    Ntot.breed.raw[scen,1]<-Ntot.breed.f[scen,1]

   
    for (tt in 2:FUT.YEAR){
    

    # -------------------------------------------------        
    # 4.1. System process for future
    # -------------------------------------------------

    ## INCLUDE CARRYING CAPACITY OF 2500 breeding pairs (slightly more than maximum ever counted)
    carr.capacity[scen,tt] ~ dnorm(2500,100)
    
    ## THE PRE-BREEDING YEARS ##
    ## because it goes for 30 years, all pops must be safeguarded to not become 0 because that leads to invald parent error
    
    nestlings.f[scen,tt] <- round(fut.fec.change[scen]*mean.fec* 0.5 * Ntot.breed.f[scen,tt])                       ### number of locally produced FEMALE chicks based on average fecundity - to use just one take ann.fec[FUT.int[tt]] 
    N1.f[scen,tt]  ~ dbin(mean.phi.juv, max(1,round(nestlings.f[scen,tt])))                                              ### number of 1-year old survivors 
    N2.f[scen,tt] ~ dbin(mean.phi.juv, max(1,round(N1.f[scen,tt])))                                                      ### number of 2-year old survivors
    N3.f[scen,tt] ~ dbin(mean.phi.juv, max(1,round(N2.f[scen,tt])))                                                       ### number of 3-year old survivors
    N4.f[scen,tt] ~ dbin(mean.phi.juv, max(1,round(N3.f[scen,tt])))                                                       ### number of 4-year old survivors
    N5.f[scen,tt] ~ dbin(fut.surv.change[scen]*mean.phi.ad, max(1,round(N4.f[scen,tt])))                                                       ### number of 5-year old survivors
    N6.f[scen,tt] ~ dbin(fut.surv.change[scen]*mean.phi.ad, max(1,round(N5.f[scen,tt])))                                 ### number of 6-year old survivors
    N7.f[scen,tt] ~ dbin(fut.surv.change[scen]*mean.phi.ad, max(1,round(N6.f[scen,tt])))                                 ### number of 7-year old survivors
    N8.f[scen,tt] ~ dbin(fut.surv.change[scen]*mean.phi.ad, max(1,round(N7.f[scen,tt])))                                 ### number of 8-year old survivors
    N9.f[scen,tt] ~ dbin(fut.surv.change[scen]*mean.phi.ad, max(1,round(N8.f[scen,tt-1])))                               ### number of 9-year old survivors
    

    ## THE BREEDING POPULATION ##

    # Ntot.breed comprised of first-time breeders, previous skippers, and previous unsuccessful breeders
    N.succ.breed.f[scen,tt] ~ dbin(fut.fec.change[scen]*mean.fec, max(1,round(Ntot.breed.f[scen,tt])))                     ### number of successful breeders not breeding the following year
    N.unsucc.breed.f[scen,tt] <- round(Ntot.breed.f[scen,tt]-N.succ.breed.f[scen,tt])
    N.fail.skip.f[scen,tt] ~ dbin(mean.skip, max(1,N.unsucc.breed.f[scen,tt]))   ### number of unsuccessful breeders not breeding the following year
    N.nonbreed.f[scen,tt] <- round(sum(N.succ.breed.f[scen,tt-1],N.fail.skip.f[scen,tt-1])) ### number of nonbreeders as the sum of previous year successful and skipping unsuccessful breeders
    N.non.breed.f[scen,tt] ~ dbin(fut.surv.change[scen]*mean.phi.ad, max(1,N.nonbreed.f[scen,tt]))  ### number of non-breeders in a given year is composed of the successful skippers and failed skippers from previous year       
    
    N.breed.ready.f[scen,tt]<- round((Ntot.breed.f[scen,tt-1]-N.succ.breed.f[scen,tt-1]-N.fail.skip.f[scen,tt-1])+N9.f[scen,tt-1]+N.non.breed.f[scen,tt-1])       ### number of available breeders is failed breeders from previous year plus non-breeders from previous year plus recruits (N9)
    Ntot.breed.raw[scen,tt] ~ dbin(fut.surv.change[scen]*mean.phi.ad, max(1,N.breed.ready.f[scen,tt]))   ### number of potential old breeders is the number of survivors from previous year breeders and nonbreeders
    Ntot.breed.f[scen,tt] <- min(Ntot.breed.raw[scen,tt],carr.capacity[scen,tt])  ## simply reduces number of breeders to carrying capacity

    } ### end future loop
    
    ## CALCULATE ANNUAL POP GROWTH RATE ##
      for (fut2 in 1:(FUT.YEAR-1)){
        fut.lambda[scen,fut2] <- Ntot.breed.f[scen,fut2+1]/max(1,Ntot.breed.f[scen,fut2])                                 ### inserted safety to prevent denominator being 0
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
R[1,]<-round(colMeans(R[4:6,], na.rm=T),0)
sum(R[1,])   ## reduces count from 2400 to 1638

# Bundle data
jags.data <- list(marr.j = chick.marray,
                  marr.a = adult.marray,
                  n.occasions = dim(chick.marray)[2],
                  r.j=apply(chick.marray,1,sum),
                  r.a=apply(adult.marray,1,sum),
                  
                  ### count data
                  T = n.years,
                  prop.sites=mean.props,
                  y.count=R,
                  
                  ### breeding success data
                  n.sites=n.sites,
                  J=J,
                  R=R,
                  
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
                         mean.phi.juv = runif(1, 0, 1),
                         mean.p = runif(1, 0, 1),
                         Ntot.breed= c(1638,rep(NA,n.years-1)),
                         #bycatch = rnorm(1,0,0.01),
                         #hookpod = rnorm(1,0,0.01),
                         
                         ### count data
                         #sigma.proc=runif(n.sites,0,10),
                         #mean.lambda=runif(n.sites,0.1,2),
                         #N.est=N.init,  ## results in 'inconsistent with parents' error
                         sigma.obs=runif(n.sites,1,100))}

 

# Parameters monitored
parameters <- c("Ntot.breed","Ntot.breed.f","ann.fec","mean.phi.ad","mean.phi.juv","mean.fec","mean.skip","mean.p","pop.growth.rate","fut.growth.rate")  

# MCMC settings
ni <- 150000
nt <- 2
nb <- 500
nc <- 3



# RUN THE FOUR SCENARIOS {took 20 hours for niter=150000)
TRALipm <- autojags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM\\TRAL_IPM_fut_marray_v1.jags",
                    n.chains = nc, n.thin = nt, n.burnin = nb,parallel=T,#n.iter = ni)
                    Rhat.limit=1.5, max.iter=25000)  



#########################################################################
# SAVE OUTPUT - RESULT PROCESSING in TRAL_IPM_result_summaries.r
#########################################################################
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
save.image("TRAL_IPM_output_v5.RData")





















#########################################################################
# ESTIMATING STABLE STAGE DISTRIBUTION WITH SIMPLE DETERMINISTIC MODEL ##
#########################################################################

#################### CREATING THE POPULATION MATRIX #####################
## transition probabilities are FROM (col) TO (row)
library(popbio)
seabird.matrix<-expression(
  0,0,0,0,0,0,0,0,0,0,(F*0.5*S2),
  S1,0,0,0,0,0,0,0,0,0,0,
  0,S1,0,0,0,0,0,0,0,0,0,
  0,0,S1,0,0,0,0,0,0,0,0,
  0,0,0,S1,0,0,0,0,0,0,0,
  0,0,0,0,S1,0,0,0,0,0,0,
  0,0,0,0,0,S2,0,0,0,0,0,
  0,0,0,0,0,0,S2,0,0,0,0,
  0,0,0,0,0,0,0,S2,0,0,0,
  0,0,0,0,0,0,0,0,S2,0,0,
  0,0,0,0,0,0,0,0,0,(S2*F+S2*(1-F)*0.3),S2*(1-F)*0.7
)



### CREATE LESLIE MATRIX WITH SUBSET OF VITAL RATES

seabird.vr<-list(F=0.32, S1=0.80,S2=0.93)
A<-matrix(sapply(seabird.matrix, eval,seabird.vr , NULL), nrow=sqrt(length(seabird.matrix)), byrow=TRUE)
stable.stage(A)



#### TROUBLESHOOT INCONSISTENT ERROR PROBLEM ####
rmultinom(n=rep(100,5), size=runif(100,600,700), prob=jags.data$prop)

