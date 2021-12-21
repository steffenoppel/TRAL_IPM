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
library(grid)
library(magick)


#########################################################################
# LOAD MODEL OUTPUT FROM IPMs
#########################################################################

setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
#load("TRAL_IPM_output_2020.RData")
#load("TRAL_IPM_output_v5_Ntot_agerecruit.RData")
load("TRAL_IPM_output_FINAL.RData")
ls()



#########################################################################
# IMPLEMENT LTRE using code from Koons et al. 2017
#########################################################################
### predictions created in TRAL_IPM_FINAL.r

### need the following parameters from model
lambda.r::ann.fec
p.juv
phi.juv
p.ad
phi.ad
Ntot.breed
Ntot



## # Step 1: Using the JAGS output (App. S6; named lescipm6dadwingpiece.jags),         # compute realized population growth rates for female lesser scaup. 
library(matrixStats)
samples <- 5000
noyears <- 59
lam_real <- matrix(0,samples,noyears)
F1 <- matrix(0,samples,noyears)
F2 <- matrix(0,samples,noyears)
Sfj <- matrix(0,samples,noyears)
Sfa <- matrix(0,samples,noyears)
Smj <- matrix(0,samples,noyears)
Sma <- matrix(0,samples,noyears)
for (i in 1:noyears){
  for (j in 1:samples){
    # Account for a census timing that was different than banding timing
    F1[j,i] <- lescipm6dadwingpiece$sims.list$f1[j,i+1]
    F2[j,i] <- lescipm6dadwingpiece$sims.list$f2[j,i+1]
    Sfj[j,i] <- lescipm6dadwingpiece$sims.list$juvs[j,2,i+1]
    Sfa[j,i] <- ((lescipm6dadwingpiece$sims.list$anns[j,2,i])^(1/12))^3 * 
      ((lescipm6dadwingpiece$sims.list$anns[j,2,i+1])^(1/12))^9
    Smj[j,i] <- lescipm6dadwingpiece$sims.list$juvs[j,1,i+1]
    Sma[j,i] <- ((lescipm6dadwingpiece$sims.list$anns[j,1,i])^(1/12))^3 * 
      ((lescipm6dadwingpiece$sims.list$anns[j,1,i+1])^(1/12))^9
    lam_real[j,i] <- (lescipm6dadwingpiece$sims.list$N1f[j,i+1] + 
                        lescipm6dadwingpiece$sims.list$N2f[j,i+1]) /
      (lescipm6dadwingpiece$sims.list$N1f[j,i] +   
         lescipm6dadwingpiece$sims.list$N2f[j,i])
  }
}
tempvar_real <- rowVars(lam_real)
mean(tempvar_real)
quantile(tempvar_real,0.05)
quantile(tempvar_real,0.95)

# Step 2: Calculate stage-specific proportions of abundances for each female # age class at each time step and for each of the saved MCMC samples. 
noyears <- 60   # we just used the first 59 of these below
nf1 <- matrix(0,samples,noyears)
nf2 <- matrix(0,samples,noyears)
for (i in 1:noyears){
  for (j in 1:samples){
    nf1[j,i] <- lescipm6dadwingpiece$sims.list$N1f[j,i] /
      (lescipm6dadwingpiece$sims.list$N1f[j,i] +
         lescipm6dadwingpiece$sims.list$N2f[j,i])
    nf2[j,i] <- lescipm6dadwingpiece$sims.list$N2f[j,i] /
      (lescipm6dadwingpiece$sims.list$N1f[j,i]+
         lescipm6dadwingpiece$sims.list$N2f[j,i])
  }
}
nf1 <- nf1[,1:59]
nf2 <- nf2[,1:59]

# Step 3: Calculate the transient sensitivities for each demographic         # parameter, evaluated at temporal means of each parameter. 
sens_F1 <- matrix(0,samples,1)
sens_F2 <- matrix(0,samples,1)
sens_Sfj <- matrix(0,samples,1)
sens_Sfa <- matrix(0,samples,1)
sens_nf1 <- matrix(0,samples,1)
sens_nf2 <- matrix(0,samples,1)
mu_F1 <- rowMeans(F1)
mu_F2 <- rowMeans(F2)
mu_Sfj <- rowMeans(Sfj)
mu_Sfa <- rowMeans(Sfa)
mu_nf1 <- rowMeans(nf1)
mu_nf2 <- rowMeans(nf2)
for (j in 1:samples){
  sens_F1[j] <- (0.5*mu_Sfj[j]*mu_nf1[j])
  sens_F2[j] <- (0.5*mu_Sfj[j]*mu_nf2[j])
  sens_Sfj[j] <- (0.5*mu_F1[j]*mu_nf1[j]+0.5*mu_F2[j]*mu_nf2[j])
  sens_Sfa[j] <- 1
  sens_nf1[j] <- (0.5*mu_F1[j]*mu_Sfj[j]+mu_Sfa[j]) -
    (0.5*mu_Sfj[j]*(mu_F1[j]*mu_nf1[j]+mu_F2[j]*mu_nf2[j])+mu_Sfa[j]) 
  sens_nf2[j] <- (0.5*mu_F2[j]*mu_Sfj[j]+mu_Sfa[j]) -
    (0.5*mu_Sfj[j]*(mu_F1[j]*mu_nf1[j]+mu_F2[j]*mu_nf2[j])+mu_Sfa[j])
}

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

The only difference for implementing Eq. 2 is that the sensitivities are evaluated at means between successive time steps, and that instead of computing (co)variation in a demographic parameter over time, one simply computes the difference between successive time steps.
Annotated R code for implementing key steps of the transient LTRE for change in geometric mean rates of growth between two focal periods of time (Eq. 3 of the main text):
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

