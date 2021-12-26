##########################################################################
#
# TRISTAN ALBATROSS INTEGRATED POPULATION MODEL 2001-2018
#
##########################################################################
# LTRE from Koons et al 2017
# requested by reviewer at J Appl Ecol
# implemented by steffen.oppel@rspb.org.uk on 22 December 2021

library(tidyverse)
library(jagsUI)
library(data.table)
library(runjags)
#library(nimble)
filter<-dplyr::filter
select<-dplyr::select



#########################################################################
# LOAD MODEL OUTPUT FROM IPMs
#########################################################################

setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
#load("TRAL_IPM_output_2020.RData")
#load("TRAL_IPM_output_v5_Ntot_agerecruit.RData")
load("TRAL_IPM_output_FINAL_REV2021.RData")
ls()



#########################################################################
# IMPLEMENT LTRE using code from Koons et al. 2017
#########################################################################
### predictions created in TRAL_IPM_FINAL.r

### need the following parameters from model
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="lambda[17]") # lambda: 8-24
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="phi.ad[26]") # phi.ad: 159-175
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="phi.juv[26]") # phi.juv: 201-217
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="ann.fec[17]") # ann.fec: 236 - 253
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="p.ad[26]") # p.ad: 2191-2208

## # Step 1: Using the JAGS output (named TRALipm$mcmc),         # compute realized population growth rates for Tristan Albatross 
library(matrixStats)
lam <- as.matrix(TRALipm$mcmc[[1]][,c(8:24)])
Sa <- as.matrix(TRALipm$mcmc[[1]][,c(159:175)])
Sj <- as.matrix(TRALipm$mcmc[[1]][,c(201:217)])
Fec <- as.matrix(TRALipm$mcmc[[1]][,c(236:252)])
Bp <- as.matrix(TRALipm$mcmc[[1]][,c(2191:2207)])


# Account for different start of time series for survival and fecundity data
for(ch in 2:nc){
  lam<-rbind(lam,as.matrix(TRALipm$mcmc[[ch]][,c(8:24)]))
  Sa<-rbind(Sa,as.matrix(TRALipm$mcmc[[ch]][,c(159:175)]))  ## only for the years coinciding with the fecundity and pop size
  Sj<-rbind(Sj,as.matrix(TRALipm$mcmc[[ch]][,c(201:217)]))  ## only for the years coinciding with the fecundity and pop size
  Fec<-rbind(Fec,as.matrix(TRALipm$mcmc[[ch]][,c(236:252)]))  ## without the last year for which no lambda is available
  Bp<-rbind(Bp,as.matrix(TRALipm$mcmc[[ch]][,c(2191:2207)]))  ## only for the years coinciding with the fecundity and pop size
  }

tempvar_real <- rowVars(lam)
mean(tempvar_real)
quantile(tempvar_real,0.05)
quantile(tempvar_real,0.95)




# Step 2: Extract population sizes at each time step and for each of the saved MCMC samples. 

### need the following parameters from model
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot[1]") # Ntot: 26-43
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot.breed[1]") # Ntot.breed: 218-235

Ntot <- as.matrix(TRALipm$mcmc[[1]][,c(26:43)])
Nbreed <- as.matrix(TRALipm$mcmc[[1]][,c(218:235)])
for(ch in 2:nc){
  Ntot<-rbind(Ntot,as.matrix(TRALipm$mcmc[[ch]][,c(26:43)]))
  Nbreed<-rbind(Nbreed,as.matrix(TRALipm$mcmc[[ch]][,c(218:235)]))  ## only for the years coinciding with the fecundity and pop size
}

#proportion of "breeders" of the total population size
prop_breed<-(Nbreed)/Ntot



############ PLOT RAW CORRELATIONS BETWEEN DEMOGRAPHIC PARAMETERS #####

dim(lam)
corplot.l<-data.frame(lambda=colMeans(lam),lcl.lam=apply(lam,2,quantile,probs = 0.025, na.rm = TRUE),ucl.lam=apply(lam,2,quantile,probs = 0.975, na.rm = TRUE),Year=c(2004:2020))
corplot.f<-data.frame(parm=dimnames(Fec)[[2]],value=colMeans(Fec),Year=c(2004:2020),dem="Breeding success",lcl=apply(Fec,2,quantile,probs = 0.025, na.rm = TRUE),ucl=apply(Fec,2,quantile,probs = 0.975, na.rm = TRUE))
corplot.Sa<-data.frame(parm=dimnames(Sa)[[2]],value=colMeans(Sa),Year=c(2004:2020),dem="Adult survival",lcl=apply(Sa,2,quantile,probs = 0.025, na.rm = TRUE),ucl=apply(Sa,2,quantile,probs = 0.975, na.rm = TRUE))
corplot.Sj<-data.frame(parm=dimnames(Sj)[[2]],value=colMeans(Sj),Year=c(2004:2020),dem="Juvenile survival",lcl=apply(Sj,2,quantile,probs = 0.025, na.rm = TRUE),ucl=apply(Sj,2,quantile,probs = 0.975, na.rm = TRUE))
corplot.Bp<-data.frame(parm=dimnames(Nbreed)[[2]],value=colMeans(Nbreed),Year=c(2004:2021),dem="N breeding pairs",lcl=apply(Nbreed,2,quantile,probs = 0.025, na.rm = TRUE),ucl=apply(Nbreed,2,quantile,probs = 0.975, na.rm = TRUE))

corplot<-bind_rows(corplot.f,corplot.Sa,corplot.Sj,corplot.Bp) %>%
  left_join(corplot.l, by="Year")

### Pearson correlation tests
test.vals<-corplot %>% group_by(dem) %>% summarise(x=min(lcl)) %>% mutate(text=NA)
ct.f<-cor.test(corplot.l$lambda,corplot.f$value)
test.vals$text[2]<-paste("r = ",round(ct.f$estimate,2)," (",round(ct.f$conf.int[1],2),", ",round(ct.f$conf.int[2],2),")", sep="")

ct.Sa<-cor.test(corplot.l$lambda,corplot.Sa$value)
test.vals$text[1]<-paste("r = ",round(ct.Sa$estimate,2)," (",round(ct.Sa$conf.int[1],2),", ",round(ct.Sa$conf.int[2],2),")", sep="")

ct.Sj<-cor.test(corplot.l$lambda,corplot.Sj$value)
test.vals$text[3]<-paste("r = ",round(ct.Sj$estimate,2)," (",round(ct.Sj$conf.int[1],2),", ",round(ct.Sj$conf.int[2],2),")", sep="")

ct.Bp<-cor.test(corplot.l$lambda,corplot.Bp$value[corplot.Bp$Year<2021])
test.vals$text[4]<-paste("r = ",round(ct.Bp$estimate,2)," (",round(ct.Bp$conf.int[1],2),", ",round(ct.Bp$conf.int[2],2),")", sep="")


## create PLOT
ggplot(corplot) + geom_point(aes(y=lambda,x=value)) +
    facet_wrap(~dem, scales="free_x") +
    geom_errorbarh(aes(y=lambda,xmin=lcl,xmax=ucl), color='grey85') +
    geom_errorbar(aes(x=value,ymin=lcl.lam,ymax=ucl.lam), color='grey85') +
    geom_point(aes(y=lambda,x=value)) +
    #geom_smooth(aes(x=lambda,y=value), method="lm") +
    geom_hline(aes(yintercept = 1), colour="darkred", size=1, linetype = "dashed") +
    #geom_abline(slope=1,colour="grey85", linetype = "dashed", size=1) +
    geom_point(aes(y=lambda,x=value)) +
    geom_text(data=test.vals,aes(y = 1.10, x = x, label = text), hjust=0, vjust = 1, size=3) +
    
    ylab("Annual population growth rate") +
    xlab("Demographic parameter value") +
    scale_y_continuous(breaks=seq(0.8,1.1,0.1), limits=c(0.75,1.11))+
    
    theme(panel.background=element_rect(fill="white", colour="black"), 
          axis.text=element_text(size=14, color="black"),
          strip.text=element_text(size=16, color="black"),
          strip.background=element_rect(fill="white", colour="black"),
          #axis.text.x=element_text(size=12, color="black", angle=45, vjust = 1, hjust=1), 
          axis.title=element_text(size=16), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank())

ggsave("C:\\STEFFEN\\MANUSCRIPTS\\Submitted\\TRAL_IPM\\Fig4_rev.jpg", width=9, height=8)








# Step 3: Calculate the transient sensitivities for each demographic parameter, evaluated at temporal means of each parameter. 
sens_Fec <- matrix(0,dim(Fec)[1],1)
sens_Sj <- matrix(0,dim(Sj)[1],1)
sens_Sa <- matrix(0,dim(Sa)[1],1)
sens_Nbreed <- matrix(0,dim(Nbreed)[1],1)

mu_Fec <- rowMeans(Fec)
mu_Sj <- rowMeans(Sj)
mu_Sa <- rowMeans(Sa)
mu_pB <- rowMeans(prop_breed)


### CODE FROM KOONS et al. 2017 - not sure how sensitivities are calculated here
### used code from Paquet et al. 2021 - adopted for TRAL without option for sub-adult breeder (replaced prop_ad with prop_breed)
### NOT SURE HOW lambda is factored into any of these equations?

for (j in 1:dim(Fec)[1]){
  # sens_fec.kestrel.pois[i,c]<-mean_phij.kestrel.pois[i,c]*mean_prop_ad.kestrel.pois[i,c]+
  #   mean_phij.kestrel.pois[i,c]*(1-mean_prop_ad.kestrel.pois[i,c])*zj.kestrel.pois$prob.breedyoung[1,i,c]
  sens_Fec[j] <- (mu_Sj[j]*mu_pB[j])
  
  # sens_phij.kestrel.pois[i,c]<-mean_fec.kestrel.pois[i,c]*mean_prop_ad.kestrel.pois[i,c]+
  #   mean_fec.kestrel.pois[i,c]*(1-mean_prop_ad.kestrel.pois[i,c])*zj.kestrel.pois$prob.breedyoung[1,i,c]
  sens_Sj[j] <- (mu_Fec[j]*mu_pB[j])
  
  # sens_phi.kestrel.pois[i,c]<-1
  sens_Sa[j] <- 1
  
  # sens_Nad.kestrel.pois[i,c]<-(mean_fec.kestrel.pois[i,c]*mean_phij.kestrel.pois[i,c]+
  #                                mean_phi.kestrel.pois[i,c])-
  #   (mean_fec.kestrel.pois[i,c]*mean_phij.kestrel.pois[i,c]*mean_prop_ad.kestrel.pois[i,c]+
  #      mean_phi.kestrel.pois[i,c]+
  #      mean_fec.kestrel.pois[i,c]*mean_phij.kestrel.pois[i,c]*zj.kestrel.pois$prob.breedyoung[1,i,c]*(1-mean_prop_ad.kestrel.pois[i,c]))
  sens_Nbreed[j]<-(mu_Fec[j]*mu_Sj[j]+mu_Sa[j])-(mu_Fec[j]*mu_Sj[j]*mu_pB[j]+mu_Sa[j])
}

### ALTERNATIVE: calculate sensitivities in popbio, probably gives completely different number
# library(popbio)
# for (j in 1:samples){
#   seabird.vr<-list(Sj=mu_Sj[j],  Sa= mu_Sa[j], Fec= mu_Fec[j], Bp=mu_Bp[j])
#   seabird.el<-expression(
#     0,  0,  Fec,
#     (0.5*Sj), 0,  0,
#     0,  (Sj^(6)*Bp), Sa)
#   x<-vitalsens(seabird.el, seabird.vr)
#   sens_Fec[j] <- x[3,2]
#   sens_Sj[j] <- x[1,2]
#   sens_Bp[j] <- x[4,2]
#   sens_Sa[j] <- x[2,2]
# }



###transient mean contributions

cont_Fec <- matrix(0,dim(Fec)[1],1)
cont_Sj <- matrix(0,dim(Sj)[1],1)
cont_Sa <- matrix(0,dim(Sa)[1],1)
cont_Nbreed <- matrix(0,dim(Nbreed)[1],1)
cont_tot <- matrix(0,dim(Nbreed)[1],1)


for (i in 1:dim(Fec)[1]){
    dp_stoch <- cbind(Fec[i,],Sa[i,],Sj[i,],prop_breed[i,1:17])
    # Derive process variance and among demographic parameters using
    # 'shrinkage' estimates of vital rates and proportionate abundances:
    dp_varcov <- var(dp_stoch)
    sensvec <- c(sens_Fec[i], sens_Sa[i], sens_Sj[i],sens_Nbreed[i])
    # calculate demographic contributions
    contmatrix <- matrix(0,4,4)
    
    for (k in 1:4){
      for (m in 1:4){
        contmatrix[k,m] <- dp_varcov[k,m]*sensvec[k]*sensvec[m]
      }#m
    }#k
    contributions <- rowSums(contmatrix)
    
    cont_Fec[i] <- contributions[1]
    cont_Sa[i] <- contributions[2]
    cont_Sj[i] <- contributions[3]
    cont_Nbreed[i] <- contributions[4]
    cont_tot[i]<-sum(contributions[])
}


# Step 5: Calculate desired statistics (e.g. the mean and Bayesian credible  # interval) from the derived posterior distributions of the LTRE contributions. 
mean(cont_tot)
quantile(cont_tot,0.05)
quantile(cont_tot,0.95)


#mean proportion of variation in growth rate explained by variation in demographic parameters
mean(cont_Fec)/mean(cont_tot)
quantile(cont_Fec,c(0.025,0.975))/quantile(cont_tot,c(0.025,0.975))

mean(cont_Sa)/mean(cont_tot)
quantile(cont_Sa,c(0.025,0.975))/quantile(cont_tot,c(0.025,0.975))

mean(cont_Sj)/mean(cont_tot)
quantile(cont_Sj,c(0.025,0.975))/quantile(cont_tot,c(0.025,0.975))

mean(cont_Nbreed)/mean(cont_tot)
quantile(cont_Nbreed,c(0.025,0.975))/quantile(cont_tot,c(0.025,0.975))


