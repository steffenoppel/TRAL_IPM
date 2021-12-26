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

cor.test()
  
  ggplot(corplot) + geom_point(aes(x=lambda,y=value)) +
    facet_wrap(~dem, scales="free_y") +
    geom_errorbar(aes(x=lambda,ymin=lcl,ymax=ucl), color='grey85') +
    geom_errorbarh(aes(y=value,xmin=lcl.lam,xmax=ucl.lam), color='grey85') +
    geom_point(aes(x=lambda,y=value)) +
    #geom_smooth(aes(x=lambda,y=value), method="lm") +
    geom_vline(aes(xintercept = 1), colour="darkred", size=1, linetype = "dashed") +
    #geom_abline(slope=1,colour="grey85", linetype = "dashed", size=1) +
    geom_point(aes(x=lambda,y=value)) +
    
    xlab("Annual population growth rate") +
    ylab("Demographic parameter value") +
    scale_x_continuous(breaks=seq(0.75,1.1,0.05), limits=c(0.75,1.11))+
    
    theme(panel.background=element_rect(fill="white", colour="black"), 
          axis.text=element_text(size=14, color="black"),
          strip.text=element_text(size=16, color="black"),
          strip.background=element_rect(fill="white", colour="black"),
          #axis.text.x=element_text(size=12, color="black", angle=45, vjust = 1, hjust=1), 
          axis.title=element_text(size=16), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank())

  ggsave("C:\\STEFFEN\\MANUSCRIPTS\\Submitted\\TRAL_IPM\\Fig4_rev.jpg", width=14, height=12)
  







# Step 3: Calculate the transient sensitivities for each demographic parameter, evaluated at temporal means of each parameter. 
sens_Fec <- matrix(0,dim(Fec)[1],1)
sens_Sj <- matrix(0,dim(Sj)[1],1)
sens_Sa <- matrix(0,dim(Sa)[1],1)
sens_Bp <- matrix(0,dim(Bp)[1],1)
sens_Ntot <- matrix(0,dim(Ntot)[1],1)
sens_Nbreed <- matrix(0,dim(Nbreed)[1],1)

mu_Fec <- rowMeans(Fec)
mu_Sj <- rowMeans(Sj)
mu_Sa <- rowMeans(Sa)
mu_Bp <- rowMeans(Bp)
mu_Ntot <- rowMeans(Ntot)
mu_Nbreed <- rowMeans(Nbreed)


### CODE FROM KOONS et al. 2017 - not sure how sensitivities are calculated here
for (j in 1:samples){
  sens_Fec[j] <- (0.5*mu_Sj[j]*mu_Nbreed[j])
  sens_Sj[j] <- (0.5*mu_Fec[j]*mu_Nbreed[j])
  sens_Bp[j] <- (0.5*mu_Sa[j]*mu_Nbreed[j])
  sens_Sa[j] <- 1
  sens_Nbreed[j] <- (0.5*mu_F1[j]*mu_Sfj[j]+mu_Sfa[j]) -
    (0.5*mu_Sfj[j]*(mu_F1[j]*mu_nf1[j]+mu_F2[j]*mu_nf2[j])+mu_Sfa[j]) 
  sens_Ntot[j] <- (0.5*mu_F2[j]*mu_Sfj[j]+mu_Sfa[j]) -
    (0.5*mu_Sfj[j]*(mu_F1[j]*mu_nf1[j]+mu_F2[j]*mu_nf2[j])+mu_Sfa[j])
}

### calculate sensitivities in popbio
library(popbio)
for (j in 1:samples){
  seabird.vr<-list(Sj=mu_Sj[j],  Sa= mu_Sa[j], Fec= mu_Fec[j], Bp=mu_Bp[j])
  seabird.el<-expression(
    0,  0,  Fec,
    (0.5*Sj), 0,  0,
    0,  (Sj^(6)*Bp), Sa)
  x<-vitalsens(seabird.el, seabird.vr)
  sens_Fec[j] <- x[3,2]
  sens_Sj[j] <- x[1,2]
  sens_Bp[j] <- x[4,2]
  sens_Sa[j] <- x[2,2]
}

### from Paquet


###transient mean sensitivities

sens_fec.kestrel.pois <- matrix(0,samples,chains)
sens_phi.kestrel.pois <- matrix(0,samples,chains)
sens_phij.kestrel.pois <- matrix(0,samples,chains)
sens_Imm.kestrel.pois <- matrix(0,samples,chains)
sens_Nad.kestrel.pois <- matrix(0,samples,chains)


for(c in 1:chains){
  for (i in 1:samples){
    sens_fec.kestrel.pois[i,c]<-mean_phij.kestrel.pois[i,c]*mean_prop_ad.kestrel.pois[i,c]+mean_phij.kestrel.pois[i,c]*(1-mean_prop_ad.kestrel.pois[i,c])*zj.kestrel.pois$prob.breedyoung[1,i,c]
    sens_phi.kestrel.pois[i,c]<-1
    sens_phij.kestrel.pois[i,c]<-mean_fec.kestrel.pois[i,c]*mean_prop_ad.kestrel.pois[i,c]+mean_fec.kestrel.pois[i,c]*(1-mean_prop_ad.kestrel.pois[i,c])*zj.kestrel.pois$prob.breedyoung[1,i,c]
    sens_Imm.kestrel.pois[i,c]<-1
    sens_Nad.kestrel.pois[i,c]<-(mean_fec.kestrel.pois[i,c]*mean_phij.kestrel.pois[i,c]+mean_phi.kestrel.pois[i,c])-(mean_fec.kestrel.pois[i,c]*mean_phij.kestrel.pois[i,c]*mean_prop_ad.kestrel.pois[i,c]+mean_phi.kestrel.pois[i,c]+mean_fec.kestrel.pois[i,c]*mean_phij.kestrel.pois[i,c]*zj.kestrel.pois$prob.breedyoung[1,i,c]*(1-mean_prop_ad.kestrel.pois[i,c]))
  }}

###transient mean contributions

cont_fec.kestrel.pois <- matrix(0,samples,chains)
cont_phi.kestrel.pois <- matrix(0,samples,chains)
cont_phij.kestrel.pois <- matrix(0,samples,chains)
cont_Imm.kestrel.pois <- matrix(0,samples,chains)
cont_Nad.kestrel.pois<- matrix(0,samples,chains)
cont_tot.kestrel.pois <- matrix(0,samples,chains)

for(c in 1:chains){
  for (i in 1:samples){
    dp_stoch <- cbind(zj.kestrel.pois$fec[2:(nyears-1),i,c],zj.kestrel.pois$phi[2,2:(nyears-1),i,c],zj.kestrel.pois$phi[1,2:(nyears-1),i,c],Immrate.kestrel.pois[2:(nyears-1),i,c],prop_ad.kestrel.pois[2:(nyears-1),i,c])
    # Derive process variance and among demographic parameters using
    # 'shrinkage' estimates of vital rates and proportionate abundances:
    dp_varcov <- var(dp_stoch)
    sensvec <- c( sens_fec.kestrel.pois[i,c], sens_phi.kestrel.pois[i,c], sens_phij.kestrel.pois[i,c],sens_Imm.kestrel.pois[i,c],sens_Nad.kestrel.pois[i,c])
    # calculate demographic contributions
    contmatrix <- matrix(0,5,5)
    
    for (k in 1:5){
      for (m in 1:5){
        contmatrix[k,m] <- dp_varcov[k,m]*sensvec[k]*sensvec[m]
      }#m
    }#k
    contributions <- rowSums(contmatrix)
    
    cont_fec.kestrel.pois[i,c] <- contributions[1]
    cont_phi.kestrel.pois[i,c] <- contributions[2]
    cont_phij.kestrel.pois[i,c] <- contributions[3]
    cont_Imm.kestrel.pois[i,c] <- contributions[4]
    cont_Nad.kestrel.pois[i,c] <- contributions[5]
    cont_tot.kestrel.pois[i,c]<-sum(contributions[])
  }}

mean(cont_Imm.kestrel.pois)
quantile(cont_Imm.kestrel.pois,c(0.025,0.975))

#mean proportion of variation in growth rate explained by variation in immigration rate (can't be properly calculated on the posteriors)
mean(cont_Imm.kestrel.pois)/mean(cont_tot.kestrel.pois)

#Note that the LTRE contribution of "the other demographic rates, as shown in Fig. 4 can be obtained by subtracting the contribution of immigration from the total contribution
mean(cont_tot.kestrel.pois-cont_Imm.kestrel.pois)
quantile(cont_tot.kestrel.pois-cont_Imm.kestrel.pois,c(0.025,0.975))

####




# Step 4: Calculate the LTRE contributions of temporal process (co)variation # in the demographic parameters to temporal variation in the realized 
# population growth rates.
cont_F1 <- matrix(0,samples,1)
cont_F2 <- matrix(0,samples,1)
cont_Sfj <- matrix(0,samples,1)
cont_Sfa <- matrix(0,samples,1)
cont_nf1 <- matrix(0,samples,1)
cont_nf2 <- matrix(0,samples,1)
for (j in 1:samples){
  dp_stoch <- cbind(F1[j,],F2[j,],Sfj[j,],Sfa[j,],nf1[j,],nf2[j,])
  # Derive process variance and covariance among demographic parameters using
  # 'shrinkage' estimates of vital rates and proportionate abundances:
  dp_varcov <- var(dp_stoch)
  sensvec <- c(sens_F1[j],sens_F2[j],sens_Sfj[j],sens_Sfa[j], 
               sens_nf1[j],sens_nf2[j])
  # calculate demographic contributions
  contmatrix <- matrix(0,6,6)
  for (k in 1:6){
    for (m in 1:6){
      contmatrix[k,m] <- dp_varcov[k,m]*sensvec[k]*sensvec[m]
    }
  }
  contributions <- rowSums(contmatrix)
  cont_F1[j] <- contributions[1]
  cont_F2[j] <- contributions[2]
  cont_Sfj[j] <- contributions[3]
  cont_Sfa[j] <- contributions[4]
  cont_nf1[j] <- contributions[5]
  cont_nf2[j] <- contributions[6]
}

# Step 5: Calculate desired statistics (e.g. the mean and Bayesian credible  # interval) from the derived posterior distributions of the LTRE             # contributions. The following is an example for the total contribution      # from all demographic parameters.
totalcont <- cont_F1+cont_F2+cont_Sfj+cont_Sfa+cont_nf1+cont_nf2
mean(totalcont)
quantile(totalcont,0.05)
quantile(totalcont,0.95)

#The only difference for implementing Eq. 2 is that the sensitivities are evaluated at means between successive time steps, and that instead of computing (co)variation in a demographic parameter over time, one simply computes the difference between successive time steps.
#Annotated R code for implementing key steps of the transient LTRE for change in geometric mean rates of growth between two focal periods of time (Eq. 3 of the main text):
  # Step 1: Provide the symbolic matrix structure, and calculate symbolic      # derivatives of the matrix with respect to change in lower-level vital      # rates.
  matrix.elements <- expression(0.5*F1*Sfj, 0.5*F2*Sfj,
                                Sfa,        Sfa)
dF1  <- as.expression(sapply(matrix.elements,D,"F1"))
dF2  <- as.expression(sapply(matrix.elements,D,"F2"))
dSfj  <- as.expression(sapply(matrix.elements,D,"Sfj"))
dSfa  <- as.expression(sapply(matrix.elements,D,"Sfa"))

# Step 2: Compute geometric mean population growth rates for each time period # of interest, then the difference (time periods must be of equal duration). # Here we focus on the comparison of the 2006-2016 time period to that of    # 1996-2006.
samples <- 5000
loggeolam_1 <- matrix(NA,samples)
loggeolam_2 <- matrix(NA,samples)
diffgeolam <- matrix(NA,samples) 
for (j in 1:samples){
  loggeolam_1[j] <- mean(log(lam_real[j,40:49]))
  loggeolam_2[j] <- mean(log(lam_real[j,50:59]))
  diffgeolam[j] <- loggeolam_2[j] - loggeolam_1[j]
}
mean(loggeolam_1)
quantile(loggeolam_1,0.05)
quantile(loggeolam_1,0.95)
mean(loggeolam_2)
quantile(loggeolam_2,0.05)
quantile(loggeolam_2,0.95)
mean(diffgeolam)
quantile(diffgeolam,0.05)
quantile(diffgeolam,0.95)

# Step 3: Calculate population dynamics for a reference population that      # represents average initial conditions and average per time step vital rates # between the two time periods being compared.
Time <- 10
refF1 <- refF2 <- refSfj <- refSfa <- lam_realref <- matrix(NA,samples,Time)
refn <- array(NA,dim=c(2,1,samples,Time+1))
for (j in 1:samples){
  refnf1 <- (nf1[j,40] + nf1[j,50]) / 2
  refnf2 <- (nf2[j,40] + nf2[j,50]) / 2
  refn[1,1,j,1] <- refnf1
  refn[2,1,j,1] <- refnf2
  for (i in 1:Time){
    refF1[j,i] <- (F1[j,i+39] + F1[j,i+49]) / 2
    refF2[j,i] <- (F2[j,i+39] + F2[j,i+49]) / 2
    refSfj[j,i] <- (Sfj[j,i+39] + Sfj[j,i+49]) / 2
    refSfa[j,i] <- (Sfa[j,i+39] + Sfa[j,i+49]) / 2
    A <- matrix(c(
      0.5*refF1[j,i]*refSfj[j,i], 0.5*refF2[j,i]*refSfj[j,i],
      refSfa[j,i], refSfa[j,i]), nrow=2, byrow=TRUE)
    n <- A %*% refn[,,j,i]
    lam_realref[j,i] <- sum(n)
    refn[,1,j,i+1] <- n/sum(n) #store the proportionate abundances
  }
}
# temporal means of the lower-level vital rates
refF1mu <- rowMeans(refF1)
refF2mu <- rowMeans(refF2)
refSfjmu <- rowMeans(refSfj)
refSfamu <- rowMeans(refSfa)

# Step 4: Calculate differences on log scale in lower-level vital rate means # and standard deviations between time periods.
logF1mudiff <- matrix(NA,samples,1)
logF2mudiff <- matrix(NA,samples,1)
logSfjmudiff <- matrix(NA,samples,1)
logSfamudiff <- matrix(NA,samples,1)
logF1sigdiff <- matrix(NA,samples,1)
logF2sigdiff <- matrix(NA,samples,1)
logSfjsigdiff <- matrix(NA,samples,1)
logSfasigdiff <- matrix(NA,samples,1)
for (j in 1:samples){
  logF1mudiff[j] <- log(mean(F1[j,50:59])) - log(mean(F1[j,40:49]))
  logF2mudiff[j] <- log(mean(F2[j,50:59])) - log(mean(F2[j,40:49]))
  logSfjmudiff[j] <- log(mean(Sfj[j,50:59])) - log(mean(Sfj[j,40:49]))
  logSfamudiff[j] <- log(mean(Sfa[j,50:59])) - log(mean(Sfa[j,40:49]))
  logF1sigdiff[j] <- log(sqrt(var(F1[j,50:59]))) -  
    log(sqrt(var(F1[j,40:49])))
  logF2sigdiff[j] <- log(sqrt(var(F2[j,50:59]))) - 
  log(sqrt(var(F2[j,40:49])))
  logSfjsigdiff[j] <- log(sqrt(var(Sfj[j,50:59]))) - 
  log(sqrt(var(Sfj[j,40:49])))
  logSfasigdiff[j] <- log(sqrt(var(Sfa[j,50:59]))) - 
  log(sqrt(var(Sfa[j,40:49])))
}

# Step 5: Compute real-time elasticities, evaluated at the reference         # population from step 3.
# Real-time elasticities for the direct effects of change in the lower-level # vital rates.
S <- 2  # dimension of projection matrix
eAF1mu <- eAF2mu <- eASfjmu <- eASfamu <- eAF1sig <- eAF2sig <- eASfjsig <-
  eASfasig <- array(NA,dim=c(S,S,samples,Time))
tot_eAF1mu <- tot_eAF2mu <- tot_eASfjmu <- tot_eASfamu <- tot_eAF1sig <- 
  tot_eAF2sig <- tot_eASfjsig <- tot_eASfasig <- matrix(NA,samples,Time)
for (j in 1:samples){
  for (i in 1:Time){
    vr_list = list(
      F1 = refF1[j,i],
      F2 = refF2[j,i],
      Sfj = refSfj[j,i],
      Sfa = refSfa[j,i]
    )
    dA_dF1 <- matrix(sapply(dF1,eval,vr_list),ncol = 2,nrow = 2,byrow = TRUE)
    dA_dF2 <- matrix(sapply(dF2,eval,vr_list),ncol = 2,nrow = 2,byrow = TRUE)
    dA_dSfj <- matrix(sapply(dSfj,eval,vr_list),ncol = 2,nrow = 2,byrow=TRUE)
    dA_dSfa <- matrix(sapply(dSfa,eval,vr_list),ncol = 2,nrow = 2,byrow=TRUE)
    for (m in 1:S){
      for (n in 1:S){
        eAF1mu[m,n,j,i] <- refF1mu[j] * dA_dF1[m,n] * refn[n,1,j,i] / 
          lam_realref[j,i]
        tot_eAF1mu[j,i] <- sum(sum(eAF1mu[,,j,i])) 
        eAF2mu[m,n,j,i] <- refF2mu[j] * dA_dF2[m,n] * refn[n,1,j,i] / 
          lam_realref[j,i]
        tot_eAF2mu[j,i] <- sum(sum(eAF2mu[,,j,i])) 
        eASfjmu[m,n,j,i] <- refSfjmu[j] * dA_dSfj[m,n] * refn[n,1,j,i] / 
          lam_realref[j,i]
        tot_eASfjmu[j,i] <- sum(sum(eASfjmu[,,j,i])) 
        eASfamu[m,n,j,i] <- refSfamu[j] * dA_dSfa[m,n] * refn[n,1,j,i] / 
          lam_realref[j,i]
        tot_eASfamu[j,i] <- sum(sum(eASfamu[,,j,i])) 
        eAF1sig[m,n,j,i] <- (refF1[j,i] - refF1mu[j]) * dA_dF1[m,n] * 
          refn[n,1,j,i] / lam_realref[j,i]
        tot_eAF1sig[j,i] <- sum(sum(eAF1sig[,,j,i]))
        eAF2sig[m,n,j,i] <- (refF2[j,i] - refF2mu[j]) * dA_dF2[m,n] * 
          refn[n,1,j,i] / lam_realref[j,i]
        tot_eAF2sig[j,i] <- sum(sum(eAF2sig[,,j,i]))
        eASfjsig[m,n,j,i] <- (refSfj[j,i] - refSfjmu[j]) * dA_dSfj[m,n] * 
          refn[n,1,j,i] / lam_realref[j,i]
        tot_eASfjsig[j,i] <- sum(sum(eASfjsig[,,j,i])) 
        eASfasig[m,n,j,i] <- (refSfa[j,i] - refSfamu[j]) * dA_dSfa[m,n] * 
          refn[n,1,j,i] / lam_realref[j,i]
        tot_eASfasig[j,i] <- sum(sum(eASfasig[,,j,i])) 
      }
    }
  }
}
avg_eAF1mu <- rowMeans(tot_eAF1mu)
avg_eAF2mu <- rowMeans(tot_eAF2mu)
avg_eASfjmu <- rowMeans(tot_eASfjmu)
avg_eASfamu <- rowMeans(tot_eASfamu)
avg_eAF1sig <- rowMeans(tot_eAF1sig)
avg_eAF2sig <- rowMeans(tot_eAF2sig)
avg_eASfjsig <- rowMeans(tot_eASfjsig)
avg_eASfasig <- rowMeans(tot_eASfasig)

# Compute real-time elasticities for the indirect effects of past change in  # the lower-level vital rates that are channeled through perturbations to    # stage structure.
I <- diag(S)             # Identity matrix
e <- matrix(1,S,1)       # vector of 1's
# indirect elasticities
enF1mu <- enF2mu <- enSfjmu <- enSfamu <- enF1sig <- enF2sig <- enSfjsig <- 
  enSfasig <- array(0,dim=c(S,S,samples,Time))
tot_enF1mu <- tot_enF2mu <- tot_enSfjmu <- tot_enSfamu <- tot_enF1sig <-
  tot_enF2sig <- tot_enSfjsig <- tot_enSfasig <- matrix(NA,samples,Time)
for (j in 1:samples){
  # perturbation matrices
  C2F1 <- C2F2 <- C2Sfj <- C2Sfa <- C3F1 <- C3F2 <- C3Sfj <- C3Sfa <- 
    array(0,dim=c(S^2,S^2,Time)) 
  # perturbation to stage structure
  wF1mu <- wF2mu <- wSfjmu <- wSfamu <- wF1sig <- wF2sig <- wSfjsig <- 
    wSfasig <- array(0,dim=c(S,S,S,Time+1)) 
  for (i in 1:Time){
    vr_list = list(
      F1 = refF1[j,i],
      F2 = refF2[j,i],
      Sfj = refSfj[j,i],
      Sfa = refSfa[j,i]
    )
    dA_dF1 <- matrix(sapply(dF1,eval,vr_list),ncol = 2,nrow = 2,byrow = TRUE)
    dA_dF2 <- matrix(sapply(dF2,eval,vr_list),ncol = 2,nrow = 2,byrow = TRUE)
    dA_dSfj <- matrix(sapply(dSfj,eval,vr_list),ncol = 2,nrow = 2,byrow=TRUE)
    dA_dSfa <- matrix(sapply(dSfa,eval,vr_list),ncol = 2,nrow = 2,byrow=TRUE)
    mat_elements <- sapply(matrix.elements,eval,vr_list)
    A <- matrix(as.numeric(mat_elements[]),nrow=2,byrow=T)
    for (m in 1:S){
      for (n in 1:S){
        C2F1[(S+1)*m-S,(S+1)*n-S,i] <- refF1mu[j] * dA_dF1[m,n]
        C2F2[(S+1)*m-S,(S+1)*n-S,i] <- refF2mu[j] * dA_dF2[m,n]
        C2Sfj[(S+1)*m-S,(S+1)*n-S,i] <- refSfjmu[j] * dA_dSfj[m,n]
        C2Sfa[(S+1)*m-S,(S+1)*n-S,i] <- refSfamu[j] * dA_dSfa[m,n]
        C3F1[(S+1)*m-S,(S+1)*n-S,i] <- (refF1[j,i] - refF1mu[j]) * 
          dA_dF1[m,n]
        C3F2[(S+1)*m-S,(S+1)*n-S,i] <- (refF2[j,i] - refF2mu[j]) * 
          dA_dF2[m,n]
        C3Sfj[(S+1)*m-S,(S+1)*n-S,i] <- (refSfj[j,i] - refSfjmu[j]) * 
          dA_dSfj[m,n]
        C3Sfa[(S+1)*m-S,(S+1)*n-S,i] <- (refSfa[j,i] - refSfamu[j]) * 
          dA_dSfa[m,n]
        K <- I - refn[,1,j,i+1]%*%t(e)  # intermediate steps
        B <- K %*% A / lam_realref[j,i]
        gF1mu <- K %*% C2F1[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),i] %*% 
          refn[,1,j,i] / lam_realref[j,i]
        gF2mu <- K %*% C2F2[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),i] %*% 
          refn[,1,j,i] / lam_realref[j,i]
        gSfjmu <- K %*% C2Sfj[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),i] %*% 
          refn[,1,j,i] / lam_realref[j,i]
        gSfamu <- K %*% C2Sfa[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),i] %*% 
          refn[,1,j,i] / lam_realref[j,i]
        gF1sig <- K %*% C3F1[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),i] %*% 
          refn[,1,j,i] / lam_realref[j,i]
        gF2sig <- K %*% C3F2[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),i] %*% 
          refn[,1,j,i] / lam_realref[j,i]
        gSfjsig <- K %*% C3Sfj[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),i] %*% 
          refn[,1,j,i] / lam_realref[j,i]
        gSfasig <- K %*% C3Sfa[(S*m-(S-1)):(S*m),(S*n-(S-1)):(S*n),i] %*% 
          refn[,1,j,i] / lam_realref[j,i]
        wF1mu[,m,n,i+1] <- B %*% wF1mu[,m,n,i] + gF1mu
        wF2mu[,m,n,i+1] <- B %*% wF2mu[,m,n,i] + gF2mu
        wSfjmu[,m,n,i+1] <- B %*% wSfjmu[,m,n,i] + gSfjmu
        wSfamu[,m,n,i+1] <- B %*% wSfamu[,m,n,i] + gSfamu
        wF1sig[,m,n,i+1] <- B %*% wF1sig[,m,n,i] + gF1sig
        wF2sig[,m,n,i+1] <- B %*% wF2sig[,m,n,i] + gF2sig
        wSfjsig[,m,n,i+1] <- B %*% wSfjsig[,m,n,i] + gSfjsig
        wSfasig[,m,n,i+1] <- B %*% wSfasig[,m,n,i] + gSfasig
        enF1mu[m,n,j,i] <- t(e) %*% A %*% wF1mu[,m,n,i] / lam_realref[j,i] 
        tot_enF1mu[j,i] <- sum(sum(enF1mu[,,j,i]))
        enF2mu[m,n,j,i] <- t(e) %*% A %*% wF2mu[,m,n,i] / lam_realref[j,i] 
        tot_enF2mu[j,i] <- sum(sum(enF2mu[,,j,i]))
        enSfjmu[m,n,j,i] <- t(e) %*% A %*% wSfjmu[,m,n,i] / lam_realref[j,i] 
        tot_enSfjmu[j,i] <- sum(sum(enSfjmu[,,j,i]))
        enSfamu[m,n,j,i] <- t(e) %*% A %*% wSfamu[,m,n,i] / lam_realref[j,i] 
        tot_enSfamu[j,i] <- sum(sum(enSfamu[,,j,i]))
        enF1sig[m,n,j,i] <- t(e) %*% A %*% wF1sig[,m,n,i] / lam_realref[j,i] 
        tot_enF1sig[j,i] <- sum(sum(enF1sig[,,j,i]))
        enF2sig[m,n,j,i] <- t(e) %*% A %*% wF2sig[,m,n,i] / lam_realref[j,i] 
        tot_enF2sig[j,i] <- sum(sum(enF2sig[,,j,i]))
        enSfjsig[m,n,j,i] <- t(e) %*% A %*% wSfjsig[,m,n,i] / 
          lam_realref[j,i] 
        tot_enSfjsig[j,i] <- sum(sum(enSfjsig[,,j,i]))
        enSfasig[m,n,j,i] <- t(e) %*% A %*% wSfasig[,m,n,i] / 
          lam_realref[j,i] 
        tot_enSfasig[j,i] <- sum(sum(enSfasig[,,j,i]))
      }
    }
  }
}
avg_enF1mu <- rowMeans(tot_enF1mu)
avg_enF2mu <- rowMeans(tot_enF2mu)
avg_enSfjmu <- rowMeans(tot_enSfjmu)
avg_enSfamu <- rowMeans(tot_enSfamu)
avg_enF1sig <- rowMeans(tot_enF1sig)
avg_enF2sig <- rowMeans(tot_enF2sig)
avg_enSfjsig <- rowMeans(tot_enSfjsig)
avg_enSfasig <- rowMeans(tot_enSfasig)

# Step 6: Calculate vital rate contributions to the difference in geometric  # mean rates of population growth between time periods. This is a function of # logged differences in mean of vital rates, logged differences in s.d. of   # vital rates, as channeled through direct effects of perturbations to the   # moments and indirect effects of perturbations to these moments channeled   # through changes in population structure over time.
contAF1mu <- contAF2mu <- contASfjmu <- contASfamu <- contAF1sig <- contAF2sig <- contASfjsig <- contASfasig <- contnF1mu <- contnF2mu <- contnSfjmu <- contnSfamu <- contnF1sig <- contnF2sig <- contnSfjsig <- contnSfasig <- matrix(NA,samples,1)
for (j in 1:samples){
  contAF1mu[j] <- logF1mudiff[j] * avg_eAF1mu[j] 
  contAF2mu[j] <- logF2mudiff[j] * avg_eAF2mu[j]
  contASfjmu[j] <- logSfjmudiff[j] * avg_eASfjmu[j]
  contASfamu[j] <- logSfamudiff[j] * avg_eASfamu[j]
  contAF1sig[j] <- logF1sigdiff[j] * avg_eAF1sig[j]
  contAF2sig[j] <- logF2sigdiff[j] * avg_eAF2sig[j]
  contASfjsig[j] <- logSfjsigdiff[j] * avg_eASfjsig[j]
  contASfasig[j] <- logSfasigdiff[j] * avg_eASfasig[j]
  contnF1mu[j] <- logF1mudiff[j] * avg_enF1mu[j]
  contnF2mu[j] <- logF2mudiff[j] * avg_enF2mu[j]
  contnSfjmu[j] <- logSfjmudiff[j] * avg_enSfjmu[j]
  contnSfamu[j] <- logSfamudiff[j] * avg_enSfamu[j]
  contnF1sig[j] <- logF1sigdiff[j] * avg_enF1sig[j]
  contnF2sig[j] <- logF2sigdiff[j] * avg_enF2sig[j]
  contnSfjsig[j] <- logSfjsigdiff[j] * avg_enSfjsig[j]
  contnSfasig[j] <- logSfasigdiff[j] * avg_enSfasig[j]
}
# One can retrieve means and quantiles from the posterior distributions for  # each component contribution. Here we provide an example of estimating net  # contributions from each vital rate as "derived parameters", then retrieving # means and quantiles from respective posterior distributions. 
totF1 <- contAF1mu + contAF1sig + contnF1mu + contnF1sig
totF2 <- contAF2mu + contAF2sig + contnF2mu + contnF2sig
totF <- totF1 + totF2
totSfj <- contASfjmu + contASfjsig + contnSfjmu + contnSfjsig
totSfa <- contASfamu + contASfasig + contnSfamu + contnSfasig
mean(totF1)
mean(totF2)
mean(totF)
mean(totSfj)
mean(totSfa)
quantile(totF1,0.05)
quantile(totF1,0.95)
quantile(totF2,0.05)
quantile(totF2,0.95)
quantile(totF,0.05)
quantile(totF,0.95)
quantile(totSfj,0.05)
quantile(totSfj,0.95)
quantile(totSfa,0.05)
quantile(totSfa,0.95)

