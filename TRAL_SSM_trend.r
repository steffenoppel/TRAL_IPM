##########################################################################
#
# TRISTAN ALBATROSS POPULATION TREND 2001-2019
#
##########################################################################
# based on Kery and Schaub 2012, Chapter 5
# modified by Steffen oppel, JULY 2019
# explored various formulations on 16 July 2019
# using all sites produces rubbish output as some counts are crazy
# used originally 5 sites (Albatross Plain, Gonydale, Tarn Moss, Tafelkop, Triple Peak) produces somewhat sensible output
# included derived variable to quantify proportion of total pop in each site

## updated 17 July 2019 because Chris indicated that we cannot drop sites due to inconsistent counts - as these may have in fact the poorest breeding success
## instead lumped counts into island total and ran model over island total (no site-specific loop) (v2) - never converged!
## reverted to site-specific approach but lumping West Point and GP Valley and Gonydale, Hummocks, Green Hill

## NOTHING CONVERGES THAT USES THE DUBIOUS DATA INCLUDING GREEN HILL

library(tidyverse)
library(jagsUI)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select




#########################################################################
# LOAD PRE-PREPARED DATA
#########################################################################
### see 'IPM_DATA_PREPARATION.R' for details on how data are aggregated
### COUNT DATA FOR POPULATION TREND ######

TRAL.pop<-fread("TRAL_INCU_counts_2001_2018.csv")
head(TRAL.pop)
names(TRAL.pop)
TRAL.pop[14,4]<-136   ### number of nests monitored in Gonydale that year

## CALCULATE SUM PER YEAR 
#TRAL.pop$tot<-rowSums(TRAL.pop[,2:12], na.rm=T)
#TRAL.pop$tot[c(2,3,11)]<-NA

## COMBINE SITES THAT WERE AMBIGUOUSLY DEFINED OVER TIME 
TRAL.pop<-TRAL.pop %>% gather(key='Site', value='Count',-Year) %>%
  mutate(Site=ifelse(Site %in% c('GP Valley','West Point'),'GP Valley',Site)) %>%
  mutate(Site=ifelse(Site %in% c('Gonydale','Green Hill','Hummocks'),'Gonydale',Site)) %>%
  group_by(Year,Site) %>%
  summarise(Count=sum(Count, na.rm=T)) %>%
  mutate(Count=ifelse(Count==0,NA,Count)) %>%
  spread(key=Site, value=Count)

R<- as.matrix(TRAL.pop[4:19,c(2,4,5,7,8,9)])
# #R[1,11]<- 350 ### fill in this so that there is no need to anchor the state-space model to an NA value
# TRAL.pop$tot<-rowSums(TRAL.pop[,c(2,4,9,10,11)], na.rm=T)
n.years<-nrow(R)
n.sites<-ncol(R)


#### INCLUDE PROJECTION INTO FUTURE ####
#fut<- 20  ## number of years to predict into future
#R<-rbind(R, matrix(NA,nrow=fut,ncol=n.sites))
#n.years<-nrow(R)



#########################################################################
# MANIPULATE DATA: INITIAL VALUES
#########################################################################

N.init=matrix(NA, nrow=n.years,ncol=n.sites)
N.init[1,]<-as.matrix(R[1,])




#########################################################################
# SPECIFY MODEL IN JAGS
#########################################################################
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
sink("TRAL_SSM_v1.jags")
cat("

  
    model {
    #-------------------------------------------------
    # state space model for population counts
    # -------------------------------------------------
    
    #-------------------------------------------------  
    # 1. PRIORS FOR POPULATION COUNTS
    #-------------------------------------------------
    
    for (s in 1:n.sites){			### start loop over every study area
      N.est[1,s] ~ dunif(0,400)   ## draw random value from a uniform distribution between 0 and 200 for initial population size
      mean.lambda[s] ~ dunif(0,10)	#Prior for mean growth rate
      sigma.proc[s] ~ dunif(0,10)	#Prior for SD of state process (annual variation in pop size)
      tau.proc[s]<-pow(sigma.proc[s],-2)
      sigma.obs[s] ~ dunif(0,100)	#Prior for SD of observation process (variation in detectability)
      tau.obs[s]<-pow(sigma.obs[s],-2)
    }
    
    
    #-------------------------------------------------  
    # 2. LIKELIHOOD FOR POPULATION COUNT DATA
    #-------------------------------------------------
    
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
    
    
    #-------------------------------------------------  
    # 3. DERIVED DATA FOR OUTPUT
    #-------------------------------------------------
    
    ## DERIVED POPULATION SIZE PER YEAR
    for (t in 1:T){
      pop.size[t]<-sum(N.est[t,1:n.sites])
    
        for (s in 1:n.sites) {
          prop[t,s]<-N.est[t,s]/pop.size[t]
        } ## end site loop
    }		## end year loop
    
    pop.growth.rate<-mean(mean.lambda[])
    
    
}      ## END MODEL 
    
    
    
    ",fill = TRUE)
sink()





#########################################################################
# PREPARE DATA FOR MODEL
#########################################################################

# Bundle data
jags.data <- list(T = n.years,
                  n.sites=n.sites,
                  R=R)


# Initial values 
inits <- function(){list(sigma.proc=runif(n.sites,0,10),
                          mean.lambda=runif(n.sites,0.1,2),
                          sigma.obs=runif(n.sites,0,100))}
 

# Parameters monitored
parameters <- c("pop.size","prop","pop.growth.rate")  

# MCMC settings
ni <- 1000000
nt <- 5
nb <- 500000
nc <- 4



# RUN THE FOUR SCENARIOS
TRALssm <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM\\TRAL_SSM_v1.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)
TRALssm



########################################################################
# PRODUCE OUTPUT TABLE
#########################################################################


## write output into file ##
pop.count<-as.data.frame(TRALipm$summary) %>% select(c(1,5,2,3,7,8)) %>%
  setNames(c('Mean', 'Median','SD','lcl', 'ucl','Rhat')) %>%
  mutate(parameter=row.names(TRALipm$summary)) %>%
  mutate(Year=c(rep(seq(2001,2019,1),6),rep(NA,2))) %>%
  filter(grepl("pop.size",parameter,perl=T,ignore.case = T)) %>%
  arrange(Year)  

pop.props<-as.data.frame(TRALipm$summary) %>% select(c(1,5,2,3,7,8)) %>%
  setNames(c('Mean', 'Median','SD','lcl', 'ucl','Rhat')) %>%
  mutate(parameter=row.names(TRALipm$summary)) %>%
  mutate(Year=c(rep(seq(2001,2019,1),6),rep(NA,2))) %>%
  filter(grepl("prop",parameter,perl=T,ignore.case = T)) %>%
  mutate(site=ifelse(grepl(",1]",parameter,perl=T,ignore.case = T)==T,"AlbatrossPlain",NA)) %>%
  mutate(site=ifelse(grepl(",2]",parameter,perl=T,ignore.case = T)==T,"Gonydale",site)) %>%
  mutate(site=ifelse(grepl(",3]",parameter,perl=T,ignore.case = T)==T,"Tafelkop",site)) %>%
  mutate(site=ifelse(grepl(",4]",parameter,perl=T,ignore.case = T)==T,"TarnMoss",site)) %>%
  mutate(site=ifelse(grepl(",5]",parameter,perl=T,ignore.case = T)==T,"TriplePeak",site)) %>%
  group_by(site) %>%
  summarise(prop=mean(Median))
pop.props


as.data.frame(TRALipm$summary) %>% select(c(1,5,2,3,7,8)) %>%
  setNames(c('Mean', 'Median','SD','lcl', 'ucl','Rhat')) %>%
  mutate(parameter=row.names(TRALipm$summary)) %>%
  filter(grepl("pop.growth",parameter,perl=T,ignore.case = T))


#########################################################################
# PRODUCE OUTPUT GRAPH THAT SHOWS ESTIMATES FOR POPULATION TREND
#########################################################################
#pdf("TRAL_pop_trend_Gough_2001_2019.pdf", width=14, height=8)

  ggplot(pop.count) + 
  geom_line(aes(y=Median, x=Year), size=1)+
  geom_ribbon(aes(x=Year, ymin=lcl,ymax=ucl),alpha=0.2)+
  geom_point(data=TRAL.pop,aes(x=Year, y=tot),col='firebrick')+
  ylab("Number of Tristan Albatross pairs") +
  scale_y_continuous(breaks=seq(0,800,100), limits=c(0,800))+
  scale_x_continuous(breaks=seq(2001,2019,2))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=16, color="black"), 
        axis.title=element_text(size=18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())


dev.off()




#########################################################################
# SAVE OUTPUT - RESULT PROCESSING in TRAL_IPM_result_summaries.r
#########################################################################
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
save.image("TRAL_SSM_output_basic.RData")

