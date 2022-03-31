##########################################################################
#
# TRISTAN ALBATROSS INTEGRATED POPULATION MODEL 2001-2050
#
##########################################################################
## written by steffen.oppel@rspb.org.uk
## code to replicate analysis published in:
## Oppel, Clark, Risi, Horswill, Converse, Jones, Osborne, Stevens, Perold, Bond, Wanless, Cuthbert, Cooper, and Ryan. 2022."Demographic consequences of invasive species predation and management on the population of a long-lived seabird species"

library(tidyverse)
library(lubridate)
library(data.table)
library(runjags)
filter<-dplyr::filter
select<-dplyr::select



#########################################################################
# LOAD ALL PREPARED DATA
#########################################################################
#setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
#load("Oppel_etal_TristanAlbatross_IPM_Input.RData")
file_url <- "https://github.com/steffenoppel/TRAL_IPM/blob/main/Oppel_etal_TristanAlbatross_IPM_Input.RData?raw=true"
load(url(file_url))
ls()

## Description of data:
# R: matrix with count data of Tristan Albatross pairs from 2004-2021 (rows) in 8 study areas (columns)
# PROD.DAT: dataframe with breeding success data of Tristan Albatross from 2004-2021, where column R denotes the total number of breeding pairs, and J the total number of fledglings
# n.sites: vector of length 1 that specifies the number of study areas in matrix R (8)
# n.years: vector of length 1 that specifies the number of years in matrix R (18)
# adult.marray: m-array matrix for adult mark recapture data with 42 rows and 43 columns specifying the number of adult birds 'released' in occasion[row] that were next observed in occasion[column]. Last column indicates adults never seen again
# chick.marray: m-array matrix for chick mark recapture data with 42 rows and 43 columns specifying the number of chicks ringed in occasion[row] that were first observed in occasion[column]. Last column indicates chicks never seen again
# fut.surv.change: dataframe that specifies the annual survival change multiplier for all future 30 years (rows) under 3 scenarios (columns)
# goodyears: dataframe that specifies the annual observation effort as the proportion of marked animals observed in that year, column 'p.sel' indicates the effort intercept for low effort (1) or high effort (2) years
# mean.props: vector of length 8 specifying the mean proportion of the annual breeding population that is counted in each of the 8 study areas in matrix R - used to relate population process to count data
# phi.juv.possible: dataframe that specifies the number of chicks ringed in each year (row) of the study, indicating whether it was possible to estimate juvenile survival (JuvSurv) for that year (=1) or not (=0, when N too low)




#########################################################################
# PREPARE DATA FOR MODEL
#########################################################################
goodyears$p.sel<-ifelse(goodyears$Contact_Year<2004,1,ifelse(goodyears$Contact_Year==2001,1,2)) ### re-calibrated which year has good and poor coverage


# Create list with input data to pass to JAGS
jags.data <- list(marr.j = chick.marray,
                  marr.a = adult.marray,
                  n.occasions = dim(chick.marray)[2],
                  r.j=apply(chick.marray,1,sum),
                  r.a=apply(adult.marray,1,sum),
                  goodyear=goodyears$p.sel,
                  juv.poss=phi.juv.possible$JuvSurv, ### sets the annual survival of juveniles to the mean if <70 were ringed
                  
                  ### count data
                  n.sites=n.sites,
                  T = n.years,
                  prop.sites=mean.props,
                  y.count=R,
                  
                  ### breeding success data
                  J=PROD.DAT$J,
                  R=PROD.DAT$R,
                  
                  ### FUTURE PROJECTION
                  FUT.YEAR=30, ## number of years for future projection
                  n.scenarios=3, ## number of future scenarios explored
                  fut.surv.change=as.matrix(fut.surv.change[,2:4]),  ## future survival rate change - matrix that adjusts gradual decrease in survival
                  fut.fec.change=c(1,2,1)     ## future fecundity change - vector with one element for each scenario
)





#########################################################################
# SPECIFY MODEL IN JAGS
#########################################################################
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
sink("Oppel_etal_TRAL_IPM.jags")
cat("

model {
    #-------------------------------------------------
    # integrated population model for the Gough TRAL population
    # - age structured model with 31 age classes 
    # - adult and juvenile survival based on CMR ringing data
    # - pre breeding census, female-based assuming equal sex ratio & survival
    # - productivity based on all areas incubator and chick counts
    # - adult breeders assumed to be detected perfectly and p in CJS model is taken as return and breed probability in population process model
    # - linked population process with SUM OF count data
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
    for (gy in 1:2){    ## for good and poor monitoring years
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
      logit(p.ad[j])  <- mu.p.ad[goodyear[j]] + eps.p[j]    #### ALTERNATIVE CONTINUOUS EFFORT CORRECTION: mu.p.ad + beta.p.eff*goodyear[j] + eps.p[j] - made no difference
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
    
    IM[1,1,1] ~ dnorm(324,20) T(0,)                                 ### number of 1-year old survivors is low because few chicks hatched in 2003
    IM[1,1,2] <- 0
    IM[1,1,3] <- IM[1,1,1] - IM[1,1,2]
    
    IM[1,2,1] ~ dnorm(257,20) T(0,)                                  ### number of 2-year old survivors is very low because very few chicks hatched in 2002
    IM[1,2,2] <- 0
    IM[1,2,3] <- IM[1,1,1] - IM[1,1,2]
    
    IM[1,3,1] ~ dnorm(462,20) T(0,)                                 ### number of 3-year old survivors is higher because many chicks hatched in 2001
    IM[1,3,2] <- 0
    IM[1,3,3] <- IM[1,1,1] - IM[1,1,2]
    
    IM[1,4,1] ~ dnorm(207,20) T(0,)                                 ### number of 4-year old survivors is very low because few chicks hatched in 2000
    IM[1,4,2] <- IM[1,4,1]*p.juv.recruit.f[4]
    IM[1,4,3] <- IM[1,1,1] - IM[1,1,2]
    
    IM[1,5,1] ~ dnorm(700,10) T(0,)                                  ### number of 5-year old survivors is huge because a lot of chicks hatched in 1999
    IM[1,5,2] <- IM[1,5,1]*p.juv.recruit.f[5]
    IM[1,5,3] <- IM[1,1,1] - IM[1,1,2]
    
    IM[1,6,1] ~ dnorm(225,20) T(0,)                                 ### very uncertain number of of 6-year old survivors because no data from 1998 or previously
    IM[1,6,2] <- IM[1,6,1]*p.juv.recruit.f[6]
    IM[1,6,3] <- IM[1,1,1] - IM[1,1,2]
    
    for(age in 7:30) {
      IM[1,age,1] ~ dbin(pow(mean.phi.ad,(age-1)), round(IM[1,age-1,3]))
      IM[1,age,2] ~ dbin(p.juv.recruit.f[age],round(IM[1,age,1]))
      IM[1,age,3] <- IM[1,age,1] - IM[1,age,2]
    }
    N.recruits[1] <- sum(IM[1,,2])  ### number of this years recruiters - irrelevant in year 1 as already included in Ntot.breed prior
    
    Ntot.breed[1] ~ dnorm(1869,100) T(0,)  ### sum of counts is 1869
    JUV[1] ~ dnorm(510,100) T(0,)          ### sum of chicks is 510
    N.atsea[1] ~ dnorm(530,20) T(0,)    ### unknown number
    Ntot[1]<-sum(IM[1,,3]) + Ntot.breed[1]+N.atsea[1]  ## total population size is all the immatures plus adult breeders and adults at sea - does not include recruits in Year 1
    
    
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
    

    for(age in 2:3) {
      IM[tt,age,1] ~ dbin(phi.ad[tt+24], max(1,round(IM[tt-1,age-1,3])))
      IM[tt,age,2] <- 0 #min(round(IM[tt,age-1,3]),IM[tt,age,1])*p.juv.recruit[age,tt]            ### prevent 2-3 year old birds from recruiting
      IM[tt,age,3] <- IM[tt,age,1] - IM[tt,age,2]
    }
    
    for(age in 4:30) {
      IM[tt,age,1] ~ dbin(phi.ad[tt+24], max(1,round(IM[tt-1,age-1,3])))
      IM[tt,age,2] ~ dbin(p.juv.recruit[age,tt],round(IM[tt,age,1]))
      IM[tt,age,3] <- IM[tt,age,1] - IM[tt,age,2]
    }
    N.recruits[tt] <- sum(IM[tt,,2])  ### number of this years recruiters
    
    
    ## THE BREEDING POPULATION ##
    # Ntot.breed comprised of first-time breeders, previous skippers, and previous unsuccessful breeders
    
    N.ad.surv[tt] ~ dbin(phi.ad[tt+24], round(Ntot.breed[tt-1]+N.atsea[tt-1]))           ### previous year's adults that survive
    breed.prop[tt] <- max(p.ad[tt+24],0.35)						### p.ad is only equivalent to breeding propensity in years with high effort
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
    
    for(age in 2:3) {
      IM.f[scen,1,age,1] ~ dbin(mean.phi.ad, max(1,round(IM[T,age-1,3])))
      IM.f[scen,1,age,2] <- 0 #min(round(IM[T,age-1,3]),IM.f[scen,1,age,1])*p.juv.recruit.f[age]  ## prevent 2-3 year olds from recruiting
      IM.f[scen,1,age,3]   <- IM.f[scen,1,age,1] - IM.f[scen,1,age,2]
    }

    for(age in 4:30) {
      IM.f[scen,1,age,1] ~ dbin(mean.phi.ad, max(1,round(IM[T,age-1,3])))
      IM.f[scen,1,age,2] ~ dbin(p.juv.recruit.f[age],round(IM.f[scen,1,age,1]))
      IM.f[scen,1,age,3]   <- IM.f[scen,1,age,1] - IM.f[scen,1,age,2]
    }
    N.recruits.f[scen,1] <- sum(IM.f[scen,1,,2])  ### number of this years recruiters


    ## THE BREEDING POPULATION ##
    
    N.ad.surv.f[scen,1] ~ dbin(mean.phi.ad, round(Ntot.breed[T]+N.atsea[T]))              ### previous year's adults that survive
    N.breed.ready.f[scen,1] ~ dbin(min(0.95,(mean.p.ad[2]/(1-mean.fec))), round(N.ad.surv.f[scen,1]))              ### number of available breeders is proportion of survivors that returns, with fecundity INCLUDED in return probability
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
    
    for(age in 2:3) {
      IM.f[scen,tt,age,1] ~ dbin(mean.phi.ad, max(1,round(IM.f[scen,tt-1,age-1,3])))
      IM.f[scen,tt,age,2] <- 0
      IM.f[scen,tt,age,3] <- IM.f[scen,tt,age,1] - IM.f[scen,tt,age,2]
    }
    
    for(age in 4:30) {
      IM.f[scen,tt,age,1] ~ dbin(mean.phi.ad, max(1,round(IM.f[scen,tt-1,age-1,3])))
      IM.f[scen,tt,age,2] ~ dbin(p.juv.recruit.f[age], round(IM.f[scen,tt,age,1]))
      IM.f[scen,tt,age,3] <- IM.f[scen,tt,age,1] - IM.f[scen,tt,age,2]
    }
    N.recruits.f[scen,tt] <- sum(IM.f[scen,tt,,2])  ### number of this years recruiters
    
    ## THE BREEDING POPULATION ##
    N.ad.surv.f[scen,tt] ~ dbin(fut.surv.change[tt,scen]*mean.phi.ad, round((Ntot.breed.f[scen,tt-1]-N.succ.breed.f[scen,tt-1])+N.atsea.f[scen,tt-1]))           ### previous year's adults that survive
    N.prev.succ.f[scen,tt] ~ dbin(fut.surv.change[tt,scen]*mean.phi.ad, round(N.succ.breed.f[scen,tt-1] + IM.f[scen,tt-1,3,2]))     ### add recruits age 2 and 3 here             ### these birds will  remain at sea because tey bred successfully
    N.breed.ready.f[scen,tt] ~ dbin(min(0.95,(mean.p.ad[2]/(1-mean.fec))), max(1,round(N.ad.surv.f[scen,tt])))                  ### number of available breeders is proportion of unsuccessful or non-breeding survivors that returns, subtracting breeding success from return probability in best monitoring years
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
# SET UP MODEL RUN
#########################################################################

# Specify initial values 
inits <- function(){list(mean.phi.ad = runif(1, 0.7, 0.97),
                         mean.phi.juv = runif(1, 0.5, 0.9),
                         mean.p.ad = runif(2, 0.2, 1),
                         mean.p.juv = runif(2, 0, 1),
                         Ntot.breed= c(runif(1, 4950, 5050),rep(NA,n.years-1)),
                         JUV= c(rnorm(1, 246, 0.1),rep(NA,n.years-1)),
                         N.atsea= c(rnorm(1, 530, 0.1),rep(NA,n.years-1)),
                         sigma.obs=matrix(runif(n.sites*n.years,1,20),ncol=n.years))}

# Select the parameters monitored
parameters <- c("mean.phi.ad","mean.phi.juv","mean.fec","mean.propensity",
                "mean.recruit","pop.growth.rate","fut.growth.rate",
                "agebeta","Ntot","Ntot.f","phi.ad","phi.juv","Ntot.breed",
                "ann.fec", "sigma.obs", "mean.p.juv","mean.p.ad",
                "mean.p.sd","sigma.p","sigma.phi")

# MCMC settings
nt <- 10
nb <- 25000
nad <- 20000
nc <- 3
ns <- 200000

# RUN THE MODEL (took 3 days for ns=200000)
TRALipm <- run.jags(data=jags.data, inits=inits, parameters, 
                    model="C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM\\Oppel_etal_TRAL_IPM.jags",
                    n.chains = nc, thin = nt, burnin = nb, adapt = nad,sample = ns, 
                    method = "rjparallel") 







#########################################################################
# CREATE TABLE AND FIGURE TO SUMMARISE POPULATION TRAJECTORY
#########################################################################

## extract summary table from raw run.jags output
summary_tralipm <- summary(TRALipm)
summary_tralipm_df <- as.data.frame(summary_tralipm)
addsummary_tralipm <- add.summary(TRALipm,plots = runjags.getOption("predraw.plots"))
predictions <- data.frame(summary(addsummary_tralipm),
                          parameter = row.names(summary(addsummary_tralipm)))
row.names(predictions) <- 1:nrow(predictions)
predictions <- predictions[1:218,]   ### 200 cuts off ann.fec
predictions$Mode <- NULL
np <- names(predictions) 
names(predictions) <- c("lcl",np[2],"ucl",np[4:9],"Rhat",np[11])
max(predictions$Rhat)

## convert predictions into export table in legible format with years rather than indices
export<-predictions %>%
  mutate(Year=c(
    rep(NA,8),         ## for mean phi, p, and growth rates
    seq(2004,2021,1),   ## for N.tot
    rep(seq(2022,2051,1),each=3), ##  for Ntot.f with 3 scenarios
    #seq(2004,2020,1), ##  for lambda 
    rep(seq(1979,2020,1), 2), ##  for phi.ad and phi.juv
    seq(2004,2021,1)   ## for ann.fec
  )) %>%     ## for deviance and agebeta
  mutate(demographic=parameter) %>%
  mutate(demographic=ifelse(grepl("fec",parameter,perl=T,ignore.case = T)==T,"fecundity",demographic))%>%
  mutate(demographic=ifelse(grepl("phi",parameter,perl=T,ignore.case = T)==T,"survival",demographic))%>%
  mutate(demographic=ifelse(grepl("Ntot",parameter,perl=T,ignore.case = T)==T,"pop.size",demographic)) %>%
  mutate(demographic=ifelse(grepl("growth",parameter,perl=T,ignore.case = T)==T,"growth.rate",demographic)) %>%
  mutate(demographic=ifelse(grepl("agebeta",parameter,perl=T,ignore.case = T)==T,"agebeta",demographic)) %>%
  arrange(demographic,Year)


## PRODUCE TABLE 1 THAT SUMMARISES DEMOGRAPHIC RATES
TABLE1<-export %>% 
  filter(!grepl("Ntot",parameter)) %>%
  filter(parameter %in% c("fut.growth.rate[1]",
                          "fut.growth.rate[2]",
                          "fut.growth.rate[3]",
                          "mean.fec",
                          "pop.growth.rate",
                          "mean.phi.ad",
                          "mean.phi.juv" )) 
TABLE1




## PRODUCE FIGURE 1 THAT SHOWS POPULATION TRAJECTORY

## scenario 1: projection with no changes in demography
## scenario 2: successful mouse eradication in 2021 - fecundity doubles
## scenario 3: increasing mouse impacts on adult survival (adult survival decreases by 10%)

export %>%
  filter(grepl("Ntot",parameter,perl=T,ignore.case = T)) %>%
  filter(!grepl("Ntot.breed",parameter,perl=T,ignore.case = T)) %>%
  arrange(Year) %>%
  mutate(Scenario="no future change") %>%
  mutate(Scenario=if_else(grepl("f\\[2",parameter,perl=T,ignore.case = T), "after successful mouse eradication",if_else(grepl("f\\[3",parameter,perl=T,ignore.case = T),"no mouse eradication and worsening impacts",Scenario))) %>%
  mutate(ucl=if_else(ucl>15000,15000,ucl)) %>%
  filter(!(Median<500 & Year<2020)) %>%

  ggplot() + 
  geom_line(aes(y=Median*2, x=Year, colour=Scenario), size=1)+   #
  geom_ribbon(aes(x=Year, ymin=lcl*2,ymax=ucl*2, fill=Scenario),alpha=0.3)+ #

  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  
  scale_y_continuous(breaks=seq(0,18000,2000), limits=c(0,20000),expand = c(0, 0))+
  scale_x_continuous(breaks=seq(2005,2050,5), limits=c(2004,2050))+
  labs(x="Year", y="\nTristan Albatross Population Size (Individuals)\n",
       col="Total population scenario",
       fill="Total population scenario") +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        legend.text=element_text(size=14),
        legend.title = element_text(size=16),
        legend.position=c(0.26,0.82),
        panel.grid.major = element_line(size=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))


######################################################################################
# FIGURE TO SHOW PROBABILITY OF FUTURE POP GROWTH RATES
######################################################################################

## get all mcmc samples for fut.growth.rate[1-3]
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="fut.growth.rate[1]")
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="fut.growth.rate[2]")
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="fut.growth.rate[3]")

## collate all samples
fut.lam.samp<-data.frame()
for(ch in 1:nc){
  fut.lam.samp<-bind_rows(fut.lam.samp,as.data.frame((TRALipm$mcmc[[1]])[,5:7]))
}
head(fut.lam.samp)
dim(fut.lam.samp)


## plot histograms for future pop growth rate
fut.lam.samp %>% rename(nochange=`fut.growth.rate[1]`,erad=`fut.growth.rate[2]`,worse=`fut.growth.rate[3]`) %>%
  gather(key="Scenario",value="N") %>%
  mutate(Scenario=if_else(Scenario=="erad", "after successful mouse eradication",if_else(Scenario=="worse","no mouse eradication and worsening impacts","no future change"))) %>%
  
  ggplot(aes(x = N, fill = Scenario)) +                       # Draw overlaying histogram
  geom_histogram(position = "identity", alpha = 0.2, bins = 80, aes(y = ..density..), color="black") +
  geom_density(alpha=0.5) +
  #geom_vline(aes(xintercept = 1), colour="indianred3", size=1) +
  geom_segment(aes(x = 1, y = 0, xend = 1, yend = 50), colour="gray15", linetype = "dashed", size=1)+
  
  labs(x="Future population growth rate", y="Probability density",
       fill="Scenario") +
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14),
        legend.position=c(0.14,0.88),
        panel.grid.major = element_line(size=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))


