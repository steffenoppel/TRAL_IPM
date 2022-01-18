# Appendix S2: Pseudo code and annotated R code for implementing the transient LTREs
# Example pieces of code for implementing the transient LTRE for temporal variation in realized population growth rates. 
# Pseudo code:
#   1.	Specify the matrix population model and symbolically define the vital rates that comprise each element (row, column entry) of the matrix.
# 2.	Using the symbolic representation of the matrix model and eqn 2 from the main text, symbolically represent the realized population growth rate as a function of vital rates and components of population structure (e.g. eqn S1.1). Next, compute symbolic sensitivities of realized population growth rate to changes in vital rates and components of population structure. These will be used later in step 7 to define quantitative formula. 
# 3.	Either provide empirical data for vital rate moments (mean, variance or shape parameters for calculating moments), or specify these moments for theoretical objectives. 
# 4.	Using the data from step 3, generate time series for each vital rate using parametric distributions. Alternatively, steps 3 and 4 could be replaced with generation of environmental values and a defined (or estimated) link to vital rate outcomes (e.g. with generalized linear models). Or one could simply collate a time series of estimated vital rates.
# 5.	Use the time series of vital rates from step 4 to generate matrix population models comprised of these vital rates at each time step.
# 6.	Using an estimated or chosen initial vector of structured abundance, project population dynamics and calculate realized population growth rates in a time-variant environment with eqns 1 and 2 from the main text.
# 7.	Conduct a transient LTRE analysis of variance in realized population growth rates over time using eqns 5 and 6 from the main text. This involves the following:
#   a.	calculation of realized population growth rate sensitivities to change in vital rates and population structure (evaluated at temporal mean values) using quantitative implementation of the formula developed in step 2,
# b.	calculation of temporal covariances among demographic parameters (vital rates and components of population structure), and
# c.	the multiplication of parts from a and b to attain LTRE contributions of each demographic parameter to temporal variation in realized population growth rates.
# 
# Annotated R code:
library(tidyverse)
library(data.table)
library(runjags)
library(matrixStats)
filter<-dplyr::filter
select<-dplyr::select
library(popbio)
#rm(list=ls())


### instead of simulating vital rates, we extracted real data from IPM output of TRAL
# 
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
# load("TRAL_IPM_output_FINAL_REV2021.RData")
# 
# 
# ### need the following parameters from model
# which(dimnames(TRALipm$mcmc[[1]])[[2]]=="phi.ad[26]") # phi.ad: 159-175
# which(dimnames(TRALipm$mcmc[[1]])[[2]]=="phi.juv[26]") # phi.juv: 201-217
# which(dimnames(TRALipm$mcmc[[1]])[[2]]=="ann.fec[17]") # ann.fec: 236 - 253
# which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot[1]") # Ntot: 26-43
# which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot.breed[1]") # Ntot.breed: 218-235
# which(dimnames(TRALipm$mcmc[[1]])[[2]]=="IM[18,30,3]") # IM: 2210-3829
# 
# 
# Sa <- as.matrix(TRALipm$mcmc[[1]][,c(159:175)])
# Sj <- as.matrix(TRALipm$mcmc[[1]][,c(201:217)])
# rho <- as.matrix(TRALipm$mcmc[[1]][,c(236:252)])
# all <- as.matrix(TRALipm$mcmc[[1]][,c(26:43)])
# br <- as.matrix(TRALipm$mcmc[[1]][,c(218:235)])
# IM <- as.matrix(TRALipm$mcmc[[1]][,c(2210:3829)])
# 
# for(ch in 2:nc){
#   Sa<-rbind(Sa,as.matrix(TRALipm$mcmc[[ch]][,c(159:175)]))  ## only for the years coinciding with the fecundity and pop size
#   Sj<-rbind(Sj,as.matrix(TRALipm$mcmc[[ch]][,c(201:217)]))  ## only for the years coinciding with the fecundity and pop size
#   rho<-rbind(rho,as.matrix(TRALipm$mcmc[[ch]][,c(236:252)]))  ## without the last year for which no lambda is available
#   all<-rbind(all,as.matrix(TRALipm$mcmc[[ch]][,c(26:43)]))
#   br<-rbind(br,as.matrix(TRALipm$mcmc[[ch]][,c(218:235)]))  ## only for the years coinciding with the fecundity and pop size
# }
# nb=all-br
# rm(list= ls()[!(ls() %in% c('Sa','Sj','rho','all','br','nb','IM','matrix.elements','gr','sensrho','sensSj','sensSa','sensnb','sensbr'))])
# save.image("TRAL_LTRE_input.RData")
load("TRAL_LTRE_input.RData")
input<-data.frame(Year=seq(2004:2021),Sa=c(Sa[1,],NA), Sj=c(Sj[1,],NA),Fec=c(rho[1,],NA),n.breed=br[1,],n.nonbreed=nb[1,])
fwrite(input,"LTRE_input.csv")

# Step 1: Symbolically define the structure of the matrix population model.

# matrix.elements <- expression(Sj*(1-gamma),  rho,
#                               Sj*gamma,      Sa)

### modified to fit the structure of Tristan Albatross with juveniles,non-breeders and breeders, where 'rho' is annual breeding success
matrix.elements <- expression(0,0,(rho*0.5),
                              Sj,Sa*gamma,Sa*rho*(1-gamma),
                              0,Sa*(1-gamma), Sa*gamma*(1-rho))


# Step 2: Compute symbolic sensitivities of realized population growth rate 
# to changes in vital rates and components of population structure (which
# should all be normalized) for the matrix population model, which can later
# be used to define quantitative formula in Step 8. We did not calculate the
# sensitivity to gamma because it is later held as a fixed parameter.

# gr <- expression((Sj*(1-gamma)*nj)+(rho*na)+(Sj*gamma*nj)+(Sa*na))

### modified for Tristan Albatross with nb (non-breeders) replacing nj and br (breeders) replacing na
gr <- expression((Sa*(1-gamma)*nb)+(rho*Sj*0.5*br)+(Sa*gamma*nb)+(Sa*gamma*br*(1-rho))+(Sa*(1-gamma)*br*rho))

sensrho <- D(gr,'rho')
sensSj <- D(gr,'Sj')
sensSa <- D(gr,'Sa')
sensnb <- D(gr,'nb') ## changed from nj to nb
sensbr <- D(gr,'br') ## changed from na to br


# Step 3: Define demographic data. Either provide empirical data or specify 
# data for simulated populations. For the sake of example, we simply provide
# demographic data for 1 life history in a medium-variance environment.
# 
# Sj_mean <- 0.75 # mean of juvenile survival   
# Sj_sd <- 0.176 # standard deviation of variability over time     
# Sj_a <- 0.9 # shape parameters corresponding to sd and mean for Beta dist   
# Sj_b <- 2.633
# Sa_mean <- 0.85 # mean of adult survival   
# Sa_sd <- 0.055   
# Sa_a <- 34.841   
# Sa_b <- 6.148 
# rho_mean <- 0.31 # mean offspring recruitment
# rho_sd <- 0.926
# rho_shape <- 4 # shape parameter corresponding to sd and mean for Gamma dist   
# rho_scale <- 0.463 #scale parameter corresponding to sd and mean for 
# #Gamma distribution   
# gamma <- 0.5 # demographic parameter controlling developmental delay
# 
# # Step 4: Generate time series for each vital rate. Here we generate random 
# # variability according to a stationary iid stochastic environment. 
# 
# Time <- 25 # Time horizon for each population projection
# 
# vr_sim <- data.frame(
#   Sj = rbeta(Time, Sj_a, Sj_b),
#   Sa = rbeta(Time, Sa_a, Sa_b),
#   gamma = rep(gamma, Time),
#   rho = rgamma(Time, shape = rho_shape, scale = rho_scale)
# )



########### COMPUTE ALL FOLLOWING STEPS ONCE FOR EACH MCMC SAMPLE ########

cont_rho <- as.numeric()
cont_Sj <- as.numeric()
cont_Sa <- as.numeric()
cont_nb <- as.numeric()
cont_br <- as.numeric()
cont_tot <- as.numeric()

for (sim in 1: dim(nb)[1]) {

vr_sim <- data.frame(
  Sj = Sj[sim,],
  Sa = Sa[sim,],
  gamma = rep(0.1, dim(Sa)[2]),
  rho = rho[sim,]
)



# Step 5: Calculate matrix elements, generate matrices at each time step, and
# calculate some demographic attributes for the mean matrix.
Time=dim(Sa)[2] ### length of time series for TRAL population model
# 
# mat_elements <- sapply(matrix.elements,eval,vr_sim)
# A_list <- vector("list",Time) # creates a list of matrices
# for (t in 1:Time) A_list[[t]] <- matrix(as.numeric(mat_elements[t,]),  
#                                         nrow=2,byrow=T)
# A_mean <- mean.list(A_list) # mean matrix
# A_mean_SSD <- stable.stage(A_mean) # stable stage distribution of mean matrix
# # lower-level sensitivities and elasticities associated with the mean matrix # for later use in the asymptotic LTRE
# vr_sens_elast <- vitalsens(matrix.elements, as.list(colMeans(vr_sim)))        
# 
# 
# 
# # Step 6: Projection of population dynamics in a time-variant environment and 
# # calculation of population growth rates.
# 
# n <- A_mean_SSD  # set initial stage structure to SSD of mean matrix
# 
# 
# # Some empty matrices to store summary information
# lam <- matrix(NA,Time,1)
# lam_1 <- matrix(NA,Time,1)
# n_norm <- matrix(NA,Time+1,2)
# n_norm[1,] <- n
# 
# for (t in 1:Time){
#   n <- A_list[[t]] %*% n
#   lam[t] <- sum(n) # realized population growth rate for the time step
#   n <- n/sum(n) # normalize the population structure
#   n_norm[t+1,] <- n # store the normalized population structure
#   # asymptotic growth rate of matrix at time t
#   lam_1[t] <- max(Re(eigen(A_list[[t]])$values)) 
# }
# n_norm <- n_norm[1:Time,]

### modified for TRAL Population model by creating n_norm from IPM
# first column is proportion of non-breeders, second column is proportion of breeders

n_norm <- matrix(NA,Time+1,2)
n_norm[,1]<-(nb[sim,])/all[sim,]
n_norm[,2]<-(br[sim,])/all[sim,]



# Step 7: Conduct the transient LTRE analysis.

# Compute sensitivities of realized growth rate evaluated at temporal means 
# of each vital rate and normalized component of population structure.  

sens_rho <- mean(n_norm[,2]) * mean(vr_sim$Sj) * 0.5 - mean(vr_sim$Sa)*mean(n_norm[,2])*(0.1) + mean(vr_sim$Sa)*mean(n_norm[,2])*(1-0.1) ## should be sensrho
sens_Sj <- mean(n_norm[,2]) * 0.5 * mean(vr_sim$rho) ## should be sensSj
sens_Sa <- mean(n_norm[,1]) * (1-0.1) + mean(n_norm[,1]) * (0.1) +
  mean(n_norm[,2])*(0.1)*(1-mean(vr_sim$rho)) + mean(n_norm[,2])*(1-0.1)*(mean(vr_sim$rho))  ## should be sensSa with gamma=0.1
sens_nb <- mean(vr_sim$Sa) * (1-0.1) + mean(vr_sim$Sa) * (0.1)  ## should be sensnb
sens_br <- (mean(vr_sim$rho) * mean(vr_sim$Sj) * 0.5) + mean(vr_sim$Sa)*(0.1)*(1-mean(vr_sim$rho)) + mean(vr_sim$Sa)*(1-0.1)*(mean(vr_sim$rho))  ## should be sensbr
sensvec <- c(sens_rho,sens_Sj,sens_Sa,sens_nb,sens_br)


#### these yield more sensible results (from Paquet 2021) but I don't understand why these equations are formulated like that:
# sens_rho <- mean(n_norm[,2]) * mean(vr_sim$Sj)
# sens_Sj <- mean(n_norm[,2]) * mean(vr_sim$rho)
# sens_Sa <- 1
# sens_br<-(mean(vr_sim$rho)*mean(vr_sim$Sj)+mean(vr_sim$Sa))-(mean(vr_sim$rho)*mean(vr_sim$Sj)*mean(n_norm[,2])+mean(vr_sim$Sa))


# Compute temporal variance-covariance matrix for relevant vital rates and   # components of population structure.
dp_stoch <- cbind(vr_sim$rho,vr_sim$Sj,vr_sim$Sa,n_norm[,1],n_norm[,2])
dp_varcov <- var(dp_stoch)

# Calculate contributions of both vital rates and components of population
# structure to the variance in realized finite population growth rates over
# time.
contmatrix <- matrix(0,5,5)
for (k in 1:5){
  for(m in 1:5){
    contmatrix[k,m] <- dp_varcov[k,m]*sensvec[k]*sensvec[m]
  }
}
contributions <- rowSums(contmatrix)
cont_rho[sim] <- contributions[1]
cont_Sj[sim] <- contributions[2]
cont_Sa[sim] <- contributions[3]
cont_nb[sim] <- contributions[4]
cont_br[sim] <- contributions[5]
cont_tot[sim] <- cont_rho[sim]+cont_Sj[sim]+cont_Sa[sim]+cont_nb[sim]+cont_br[sim]

} ## end loop over all MCMC simulations



# Step 8: Calculate the mean and Bayesian credible interval from the derived posterior distributions of the LTRE contributions. 
mean(cont_tot)
quantile(cont_tot,c(0.025,0.975))


#mean proportion of variation in growth rate explained by variation in demographic parameters
mean(cont_rho)/mean(cont_tot)
quantile(cont_rho,c(0.025,0.975))/quantile(cont_tot,c(0.025,0.975))

mean(cont_Sa)/mean(cont_tot)
quantile(cont_Sa,c(0.025,0.975))/quantile(cont_tot,c(0.025,0.975))

mean(cont_Sj)/mean(cont_tot)
quantile(cont_Sj,c(0.025,0.975))/quantile(cont_tot,c(0.025,0.975))

mean(cont_br)/mean(cont_tot)
quantile(cont_br,c(0.025,0.975))/quantile(cont_tot,c(0.025,0.975))