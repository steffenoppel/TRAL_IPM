##########################################################################
#
# TRISTAN ALBATROSS INTEGRATED POPULATION MODEL 2001-2019
#
##########################################################################
# based on Kery and Schaub 2012, Chapter 11
# modified by Steffen oppel, JULY 2019
# based on AYNA IPM code
# major changes to AYNA: more gappy count data series
# adult encounter occasions span almost entire year
# breeding success is area specific


library(tidyverse)
library(jagsUI)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select



# #########################################################################
# # LOAD FISHERY DATA FROM ICCAT (n hooks 2000 - 2017)
# #########################################################################
#  see AYNA_bycatch_effort_IPM.r script for data preparation and alternatives
#
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\AYNA_IPM\\BycatchData"), silent=T)
longline<-fread("ICCAT_AYNA_overlay_nhooks_2000_2017.csv")

## scale
longlineICCAT<- (longline$EFF-mean(longline$EFF))/sd(longline$EFF)




#########################################################################
# LOAD PRE-PREPARED DATA
#########################################################################
### see 'IPM_DATA_PREPARATION.R' for details on how data are aggregated


#### CMR SURVIVAL DATA ######

try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM"), silent=T)
#try(setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM"), silent=T)
TRAL<-fread("TRAL_simple_encounter_history_2000_2019.csv")
names(TRAL)
CH<-as.matrix(TRAL[,3:21], dimnames=F)
TRAL$AGE[is.na(TRAL$AGE)]<-1    ## set all NA as 'adult'

### check that there are contacts in every season
apply(CH,2,sum)




### COUNT DATA FOR POPULATION TREND ######

TRAL.pop<-fread("TRAL_INCU_counts_2001_2018.csv")
head(TRAL.pop)
names(TRAL.pop)

R<- as.matrix(TRAL.pop[,c(4,9,10)])
#R[1,11]<- 350 ### fill in this so that there is no need to anchor the state-space model to an NA value
TRAL.pop$tot<-rowSums(TRAL.pop[,c(4,9,10)], na.rm=T)
n.years<-nrow(R)
n.sites<-ncol(R)





#### BREEDING SUCCESS DATA FOR FECUNDITY ######
TRAL.chick<-fread("TRAL_CHIC_counts_2001_2018.csv")

## detailed nest data from Gonydale only available for some years - we replace those years in the chick count matrix
TRAL.bs<-fread("TRAL_breed_success_2006_2018.csv")
TRAL.chick[match(TRAL.bs$Year,TRAL.chick$Year),4]<-as.integer(R[match(TRAL.bs$Year,TRAL.chick$Year),3]*TRAL.bs$BREED_SUCC)
TRAL.chick[14,4]<-as.integer(TRAL.bs$n_nests[7]*TRAL.bs$BREED_SUCC[7])  ## use Gonydale nest monitoring to fill in data gaps
R[14,3]<-TRAL.bs$n_nests[7]  ## use Gonydale nest monitoring to fill in data gaps
J<- as.matrix(TRAL.chick[,2:12])
R[1,11]<- 350 ### fill in this so that there is no need to anchor the state-space model to an NA value

### specify constants for JAGS
n.years<-dim(R)[1]		## defines the number of years
n.sites<-dim(R)[2]    ## defines the number of study areas




#########################################################################
# MANIPULATE DATA: CREATE MATRIX OF AGE FOR EACH OCCASION AND INDIVIDUAL
#########################################################################

## this matrix will relate to the survival parameter estimates chosen in the model
## simple model only has 2 survival parameters:
## 1 - juvenile and immature survival (years 1-5)
## 2 - adult survival (birds >5 years old)

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x==1))
f <- apply(CH, 1, get.first)


## REMOVE BIRDS THAT ARE TOO YOUNG TO HAVE HAD A CHANCE TO RETURN
tooyoung<-ifelse(f>(dim(CH)[2]-7),ifelse(TRAL$AGE==0,1,0),0)
CH<-CH[tooyoung==0,]  ## removes individuals that were ringed as chicks <7 years before end of time series
f <- apply(CH, 1, get.first)
toolate<-ifelse(f==dim(CH)[2],1,0)
CH<-CH[toolate==0,]  ## removes individuals ringed in last occasion end of time series
ages<-TRAL$AGE[tooyoung==0]

## CREATE BLANK AGE MATRIX
AGEMAT<-matrix(2,nrow=nrow(CH),ncol=ncol(CH))
n.occ<-ncol(CH)

## LOOP OVER EACH BIRD RINGED AND SET PRE-CAPTURE DATA TO NA AND ADJUST AGE
for (l in 1:nrow(AGEMAT)){
  firstocc<-get.first(CH[l,])
  lastjuv<-firstocc+6
  lastjuv<-ifelse(lastjuv>n.occ,n.occ,lastjuv)
  young<-ages[l]
  if(firstocc>1){AGEMAT[l,1:(firstocc-1)]<-NA}  ## sets everything before first contact to NA
  if(young==0){AGEMAT[l,firstocc:lastjuv]<-1}  ## sets all juvenile years to 1
}

### CHECK WHETHER IT LOOKS OK ###
head(AGEMAT)
head(CH)
tail(AGEMAT)
tail(CH)
dim(CH)


## PREPARE CONSTANTS
n.ind<-dim(CH)[1]		## defines the number of individuals
n.cmr.years<-dim(CH)[2]  ## defines the number of years
f <- apply(CH, 1, get.first)




#########################################################################
# MANIPULATE DATA: INITIAL VALUES
#########################################################################


## CREATE MATRIX for INITIAL STATE Z FOR SURVIVAL MODEL
zinit<-CH
for (l in 1:nrow(zinit)){
  firstocc<-get.first(zinit[l,])
  zinit[l,1:firstocc]<-NA  ## sets everything up to first contact to NA
  zinit[l,(firstocc+1):n.years]<-1  ## alive after first contact
}
dim(zinit)


## CREATE MATRIX for COUNTS

N.init=matrix(NA, nrow=n.years,ncol=n.sites)
N.init[1,]<-as.matrix(R[1,])



#########################################################################
# SPECIFY MODEL IN JAGS
#########################################################################
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
sink("TRAL_IPM_v1.jags")
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
    # -------------------------------------------------
    
    #-------------------------------------------------  
    # 1. PRIORS FOR ALL DATA SETS
    #-------------------------------------------------
    
    
    # -------------------------------------------------        
    # 1.1. Priors and constraints FOR FECUNDITY
    # -------------------------------------------------
    
    for (t in 1:T){  
      ann.fec[t] ~ dnorm(0.32,10) T(0,1)        ## Informative Priors on fecundity based on Wanless et al 2009
      skip.prob[t] ~ dunif(0.22,0.32)            ## PRIOR FOR ADULT FAILED BREEDER SKIPPING PROBABILITY from Wanless et al 2009
    } #t

    
    # -------------------------------------------------        
    # 1.2. Priors and constraints FOR POPULATION COUNTS
    # -------------------------------------------------
    for (s in 1:n.sites){			### start loop over every study area
      N.est[1,s] ~ dunif(0,200)   ## draw random value from a uniform distribution between 0 and 200 for initial population size
      mean.lambda[s] ~ dunif(0,10)	#Prior for mean growth rate
      sigma.proc[s] ~ dunif(0,10)	#Prior for SD of state process (annual variation in pop size)
      tau.proc[s]<-pow(sigma.proc[s],-2)
      sigma.obs[s] ~ dunif(0,100)	#Prior for SD of observation process (variation in detectability)
      tau.obs[s]<-pow(sigma.obs[s],-2)
    }



    # -------------------------------------------------        
    # 1.3. Priors and constraints FOR SURVIVAL
    # -------------------------------------------------
    
    ### RECAPTURE PROBABILITY
    mean.p ~ dunif(0, 1)                          # Prior for mean recapture
    logit.p <- log(mean.p / (1-mean.p))           # Logit transformation
    
    for (t in 1:T){
      logit(p[t]) <- logit.p  + capt.raneff[t]
      capt.raneff[t] ~ dnorm(0, tau.capt)
    }
    
    ### SURVIVAL PROBABILITY
    for (i in 1:nind){
      for (t in f[i]:(T-1)){
        logit(phi[i,t]) <- mu[AGEMAT[i,t]] + surv.raneff[t]       ###+ bycatch*longline[t]
      } #t
    } #i
    
    
    ## AGE-SPECIFIC SURVIVAL 
    for (age in 1:2){
      beta[age] ~ dunif(0.5, 1)                         # Priors for age-specific survival
      mu[age] <- log(beta[age] / (1-beta[age]))       # Logit transformation
    }
    
    ## RANDOM TIME EFFECT ON SURVIVAL 
    for (t in 1:(T-1)){
      surv.raneff[t] ~ dnorm(0, tau.surv)
    }
    
    ### PRIORS FOR RANDOM EFFECTS
    sigma.surv ~ dunif(0, 3)                     # Prior for standard deviation of survival
    tau.surv <- pow(sigma.surv, -2)
    
    sigma.capt ~ dunif(0, 3)                     # Prior for standard deviation of capture
    tau.capt <- pow(sigma.capt, -2)
    
    

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
      N1[tt]  ~ dbin(ann.surv[1,tt-1], round(JUV[tt-1]))                                                    ### number of 1-year old survivors 
      N2[tt] ~ dbin(ann.surv[1,tt-1], round(N1[tt-1]))                                                      ### number of 2-year old survivors
      N3[tt] ~ dbin(ann.surv[1,tt-1], round(N2[tt-1]))                                                       ### number of 3-year old survivors
      N4[tt] ~ dbin(ann.surv[1,tt-1], round(N3[tt-1]))                                                       ### number of 4-year old survivors
      N5[tt] ~ dbin(ann.surv[1,tt-1], round(N4[tt-1]))                                                       ### number of 5-year old survivors
      N6[tt] ~ dbin(ann.surv[2,tt-1], round(N5[tt-1]))                                                       ### number of 6-year old survivors
      N7[tt] ~ dbin(ann.surv[2,tt-1], round(N6[tt-1]))                                                       ### number of 7-year old survivors
      N8[tt] ~ dbin(ann.surv[2,tt-1], round(N7[tt-1]))                                                       ### number of 8-year old survivors
      N9[tt] ~ dbin(ann.surv[2,tt-1], round(N8[tt-1]))                                                       ### number of 9-year old survivors
    
      ## THE RECRUITERS AT AGE 10 YEARS ##
      #ann.recruits[tt] ~ dbin(ann.surv[2,tt-1], round(N9[tt-1]))                            ### number of 9-year old survivors that are ready for recruitment - using adult survival

      ## THE BREEDING POPULATION ##

      Ntot.breed[tt] ~ dpois(max(pop.size[tt],1))                                           ### the annual number of breeding birds is the estimate from the count SSM
      #pop.size[tt] ~ dpois(Ntot.breed[tt])                                           ### the annual number of breeding birds is the estimate from the count SSM
      
      # Ntot.breed comprised of first-time breeders, previous skippers, and previous unsuccessful breeders
      N.succ.breed[tt] ~ dbin(ann.fec[tt], round(Ntot.breed[tt]))                     ### number of successful breeders not breeding the following year
      N.unsucc.breed[tt] <- round(Ntot.breed[tt]-N.succ.breed[tt])
      N.fail.skip[tt] ~ dbin(skip.prob[tt], max(1,N.unsucc.breed[tt]))   ### number of unsuccessful breeders not breeding the following year
      N.nonbreed[tt] <- round(sum(N.succ.breed[tt-1],N.fail.skip[tt-1]))
      N.non.breed[tt] ~ dbin(ann.surv[2,tt-1], N.nonbreed[tt])  ### number of non-breeders in a given year is composed of the successful skippers and failed skippers from previous year       

      N.breed.ready[tt]<- round((Ntot.breed[tt-1]-N.succ.breed[tt-1]-N.fail.skip[tt-1])+N9[tt-1]+N.non.breed[tt-1])       ### number of available breeders is failed breeders from previous year plus non-breeders from previous year plus recruits (N9)
      Npot.breed[tt] ~ dbin(ann.surv[2,tt-1], N.breed.ready[tt])   ### number of potential old breeders is the number of survivors from previous year breeders and nonbreeders
      #sum(N.est[tt,1:n.sites]) ~ dbin(ann.surv[2,tt-1], N.breed.ready[tt])   ### number of potential old breeders is the number of survivors from previous year breeders and nonbreeders
      
 
    } # tt
    
    
    
    ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1 - based on deterministic 
    Ntot.breed[1] ~ dpois(pop.size[1])   
    JUV[1]<-round(Ntot.breed[1]*0.5*ann.fec[1])
    N1[1]<-round(JUV[1]*beta[1])
    N2[1]<-round(N1[1]*beta[1])
    N3[1]<-round(N2[1]*beta[1])
    N4[1]<-round(N3[1]*beta[1])
    N5[1]<-round(N4[1]*beta[1])
    N6[1]<-round(N5[1]*beta[2])
    N7[1]<-round(N6[1]*beta[2])
    N8[1]<-round(N7[1]*beta[2])
    N9[1]<-round(N8[1]*beta[2])
    N.nonbreed[1]<-round(sum(N.succ.breed[1],N.fail.skip[1]))
    N.non.breed[1]~ dbin(beta[2], N.nonbreed[1])
    N.succ.breed[1]<- round(Ntot.breed[1]*ann.fec[1])
    N.unsucc.breed[1] <- (Ntot.breed[1]-N.succ.breed[1])
    N.fail.skip[1]<- round(N.unsucc.breed[1]*skip.prob[1])


    # -------------------------------------------------        
    # 2.2. Observation process for population counts: state-space model of annual counts
    # -------------------------------------------------
    

    for (s in 1:n.sites){			### start loop over every study area
    
      ## State process for entire time series
    
      for (t in 1:(T-1)){
        lambda[t,s] ~ dnorm(mean.lambda[s], tau.proc[s])								# Distribution for random error of growth rate
        N.est[t+1,s]<-N.est[t,s]*lambda[t,s]										        # Linear predictor (population size based on past pop size and change rate)
      }														# run this loop over nyears
    
    
      ## Observation process
    
      for (t in 1:T){
        R[t,s] ~ dnorm(N.est[t,s], tau.obs[s])								# Distribution for random error in observed numbers (counts)
        }														# run this loop over t= nyears
    }		## end site loop
    
    
    
    # -------------------------------------------------        
    # 2.3. Likelihood for fecundity: Poisson regression from the number of surveyed broods
    # -------------------------------------------------
    for (s in 1:n.sites){			### start loop over every study area
      for (t in 1:(T-1)){
        J[t,s] ~ dpois(rho.fec[t,s])
        rho.fec[t,s] <- R[t,s]*ann.fec[t]
      } #	close loop over every year in which we have fecundity data
    }		## end site loop
    
    
    
    
    # -------------------------------------------------        
    # 2.4. Likelihood for adult and juvenile survival from CMR
    # -------------------------------------------------
    
    # Likelihood 
    for (i in 1:nind){

      # Define latent state at first capture
      z[i,f[i]] <- 1

      for (t in (f[i]+1):T){
    
        # State process
        z[i,t] ~ dbern(mu1[i,t])
        mu1[i,t] <- phi[i,t-1] * z[i,t-1]
    
        # Observation process
        y[i,t] ~ dbern(mu2[i,t])
        mu2[i,t] <- p[t] * z[i,t]
      } #t
    } #i
    
    
    
    
    #-------------------------------------------------  
    # 3. DERIVED PARAMETERS FOR OUTPUT REPORTING
    #-------------------------------------------------
    
    ## DERIVED SURVIVAL PROBABILITIES PER YEAR 
    for (t in 1:(T-1)){
      for (age in 1:2){
        logit(ann.surv[age,t]) <- mu[age] + surv.raneff[t]
      }
    }
    
    ## DERIVED POPULATION SIZE PER YEAR
    for (t in 1:T){
      pop.size[t]<-sum(N.est[t,1:n.sites])               ## introduced max to prevent this number from being 0 which leads to invalid parent error on Ntot.breed
    }		## end year loop
    

    ## DERIVED OVERALL POPULATION GROWTH RATE 
    pop.growth.rate<-mean(mean.lambda[]) 				# Arithmetic mean for whole time series
    
    ## DERIVED MEAN FECUNDITY 
    mean.fec <- mean(ann.fec)

    }
    
    
    
    ",fill = TRUE)
sink()





#########################################################################
# PREPARE DATA FOR MODEL
#########################################################################

# Bundle data
jags.data <- list(y = CH,
                  f = f,
                  T = n.years,
                  nind = n.ind,
                  AGEMAT=AGEMAT,
                  
                  ### count data
                  n.sites=n.sites,
                  #y.count=as.matrix(TRAL.pop[,2:12]),
                  
                  ### breeding success data
                  J=J,
                  R=R
                  
                  ### longline effort data
                  #longline=longlineICCAT,
                  
                  # ### FUTURE PROJECTION
                  #FUT.YEAR=n.years+10
                  # FUT.int=c(seq(1,(n.years-1),1),rep(NA,11)),
                  # fut.fec=c(rep(0.5,(n.years)),rep(NA,10))     ## blank vector to hold index for future demographic rates
                  )


# Initial values 
inits <- function(){list(beta = runif(2, 0.5, 1),
                         z = zinit,
                         mean.p = runif(1, 0, 1),
                         #bycatch = rnorm(1,0,0.01),
                         #hookpod = rnorm(1,0,0.01),
                         
                         ### count data
                         sigma.proc=runif(n.sites,0,10),
                         mean.lambda=runif(n.sites,0.1,2),
                         sigma.obs=runif(n.sites,0,10),
                         N.est=N.init)}
 

# Parameters monitored
parameters <- c("Ntot.breed","ann.fec","ann.surv","beta","pop.growth.rate","mean.fec","mean.p")  

# MCMC settings
ni <- 50000
nt <- 1
nb <- 20000
nc <- 4



# RUN THE FOUR SCENARIOS
TRALipm <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM\\TRAL_IPM_v1.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)




#########################################################################
# SAVE OUTPUT - RESULT PROCESSING in TRAL_IPM_result_summaries.r
#########################################################################
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
save.image("TRAL_IPM_output_basic.RData")





















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
