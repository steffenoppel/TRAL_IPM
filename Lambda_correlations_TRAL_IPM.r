#=========================================================================================================
#=========================================================================================================

#  Correlations between demographic rates and population growth

#=========================================================================================================
#=========================================================================================================
## contributed by Cat Horswill 9 Feb 2022
## modified by Steffen Oppel 9 Feb 2022
## requested by reviewer 

library(tidyverse)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select


setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
load("TRAL_IPM_output_REV2022_FINAL.RData")
#load("TRAL_LTRE_input.RData")
#load("F:/5. ZSL/Collaborations/Steffen Oppel/LTRE_input_mcmc.RData")



############## PREPARE INPUT DATA FOR CORRELATION FROM ALL MCMC SAMPLES #########################

str(TRALipm$mcmc)
retain<-parameters[c(13,4,11,12,14,7)]

### need the following parameters from model
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="lambda[1]") # lambda: 26-42
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="phi.ad[26]") # phi.ad: 177-193
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="phi.juv[26]") # phi.juv: 220-236
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="ann.fec[17]") # ann.fec: 256 - 273
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="breed.prop[17]") # breed.prop: 4-20
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot.breed[17]") # Ntot.breed: 238-254
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="IM[14,2,2]") # Ntot.breed: 238-254

as.matrix(TRALipm$mcmc[[1]])[,885:894]


retain
retaincols<-c(4:20, # breed.prop
              177:193, #phi.ad
              220:236, # phi.juv
              256:272, # ann.fec
              26:42, #lambda
              238:254) #Ntot.breed

### EXTRACT ALL FROM THE MODEL OUTPUT AND SAVE IN DIFFERENT LIST

LTRE_input_mcmc<-as.matrix(TRALipm$mcmc[[1]])[,retaincols]
str(LTRE_input_mcmc)

for(ch in 2:nc){
  LTRE_input_mcmc<-rbind(LTRE_input_mcmc,as.matrix(TRALipm$mcmc[[ch]])[,retaincols])
}

str(LTRE_input_mcmc)



############## PREPARE INPUT DATA FOR CORRELATION FROM ALL MCMC SAMPLES #########################

dimnames(LTRE_input_mcmc)[[2]]
phi.ad<-LTRE_input_mcmc[,18:34]
phi.juv<-LTRE_input_mcmc[,35:51]
ann.fec<-LTRE_input_mcmc[,52:68]
breed.prop<-LTRE_input_mcmc[,1:17]
lambda.IPM<-LTRE_input_mcmc[,69:85]
Nbreed<-LTRE_input_mcmc[,86:102]

#estimate lambda from the IPM output of Ntot
# n.years<-18
# ntot<-LTRE_input_mcmc[,76:93]
# getlam <- function(x) x[,2:n.years] / x[,1:(n.years-1)]
# lambda.IPM<-getlam(ntot)

n.draws <- nrow(phi.ad) # Determine number of MCMC draws
corr <- matrix(NA, ncol=5, nrow=n.draws) # Create object to hold results

#To accommodate for the uncertainty associated with the parameter estimates the correlation coefficients between a 
#demographic rate and the population growth rate must be calculated for each Markov chain Monte Carlo (MCMC)-
#based draw of the posterior distributions of both to obtain the posterior distribution of the
#correlation coefficient (e.g., Baillie et al., 2009; Schaub et al., 2013; Tempel et al., 2014; Saunders et al., 2018).

for (s in 1:n.draws){ # Loop over all MCMC draws and get correlations
  corr[s,1] <- cor(phi.ad[s,], lambda.IPM[s,])
  corr[s,2] <- cor(phi.juv[s,], lambda.IPM[s,])
  corr[s,3] <- cor(ann.fec[s,], lambda.IPM[s,])
  corr[s,4] <- cor(breed.prop[s,], lambda.IPM[s,])
  corr[s,5] <- cor(Nbreed[s,], lambda.IPM[s,])
}

###############################################################################
###############################################################################
# Calculate posterior summaries for the correlation coefficients
# Posterior means
apply(corr, 2, mean)
# 95% credible intervals
cri <- function(x) quantile(x, c(0.025, 0.975))
apply(corr, 2, cri)
###############################################################################
###############################################################################


#########
#PLOT
library(plotrix)


aphi<-apply(phi.ad, 2, mean)
cri.ap<-apply(phi.ad, 2, cri)

jphi<-apply(phi.juv, 2, mean)
cri.jp<-apply(phi.juv, 2, cri)

fec<-apply(ann.fec, 2, mean)
cri.f<-apply(ann.fec, 2, cri)

prop<-apply(breed.prop, 2, mean)
cri.p<-apply(breed.prop, 2, cri)

lam<-apply(lambda.IPM, 2, mean)
cri.lam<-apply(lambda.IPM, 2, cri)



layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),
       heights=c(1,1.1), widths=c(1.1,1))

par(mar=c(3,3.3,0,0), mgp=c(1.9,0.2,0))
plot(aphi, lam, type="p", ylab="Population growth rate", yaxt="n", xaxt="n", xlab="Adult survival", 
     ylim=c(0.85,1.075),xlim=c(0.82,1), cex.lab=1.1, pch=16)
plotCI(x=aphi, y=lam, liw=aphi-cri.ap[1,], uiw=cri.ap[2,]-aphi, err="x", add=TRUE, cex=1, col="dark grey", sfrac=0)
plotCI(x=aphi, y=lam, liw=lam-cri.lam[1,], uiw=cri.lam[2,]-lam, err="y", add=TRUE, cex=1, col="dark grey", sfrac=0)
points(aphi, lam,pch=16)
axis(1,at=c(0.85,0.9,0.95, 1.0),tck=0.04, col="dark grey")
axis(2,at=c(0.9,1.0,1.1),tck=0.04, col="dark grey", las=2)
axis(2,at=c(0.85,0.95,1.05),labels=c("","",""),tck=0.02, col="dark grey", las=2)

box(col="dark grey", lwd=1)
text(0.83,1.06,"A", cex=1.2)


par(mar=c(3,1,0,0))
plot(jphi, lam, type="p", ylab="", yaxt="n", xaxt="n", xlab="Juvenile survival", 
     ylim=c(0.85,1.075), xlim=c(0.5,0.95), cex.lab=1.1, pch=16)
plotCI(x=jphi, y=lam, liw=jphi-cri.jp[1,], uiw=cri.jp[2,]-jphi, err="x", add=TRUE, cex=1, col="dark grey", sfrac=0)
plotCI(x=jphi, y=lam, liw=lam-cri.lam[1,], uiw=cri.lam[2,]-lam, err="y", add=TRUE, cex=1, col="dark grey", sfrac=0)
points(jphi, lam,pch=16)
axis(1,at=c(0.5,0.7,0.9),tck=0.04, col="dark grey")
axis(2,at=c(0.9,1.0,1.1),labels=c("","",""), tck=0.04, col="dark grey", las=2)
axis(2,at=c(0.85,0.95,1.05),labels=c("","",""),tck=0.02, col="dark grey", las=2)

box(col="dark grey", lwd=1)
text(0.52,1.06,"B", cex=1.2)


par(mar=c(3,3.3,0,0), mgp=c(1.9,0.2,0))
plot(fec, lam, type="p", ylab="Population growth rate", yaxt="n", xaxt="n", xlab="Fecundity", 
     ylim=c(0.85,1.075),xlim=c(0.1,0.6), cex.lab=1.1, pch=16)
plotCI(x=fec, y=lam, liw=fec-cri.f[1,], uiw=cri.f[2,]-fec, err="x", add=TRUE, cex=1, col="dark grey", sfrac=0)
plotCI(x=fec, y=lam, liw=lam-cri.lam[1,], uiw=cri.lam[2,]-lam, err="y", add=TRUE, cex=1, col="dark grey", sfrac=0)
points(fec, lam,pch=16)
axis(1,at=c(0.2,0.4, 0.6),tck=0.04, col="dark grey")
axis(2,at=c(0.9,1.0,1.1),tck=0.04, col="dark grey", las=2)
axis(2,at=c(0.85,0.95,1.05),labels=c("","",""),tck=0.02, col="dark grey", las=2)

box(col="dark grey", lwd=1)
text(0.12,1.06,"C", cex=1.2)


par(mar=c(3,1,0,0))
plot(prop, lam, type="p", ylab="", yaxt="n", xaxt="n", xlab="Breeding propensity", 
     ylim=c(0.85,1.075), xlim=c(0.3,0.75), cex.lab=1.1, pch=16)
plotCI(x=prop, y=lam, liw=prop-cri.p[1,], uiw=cri.p[2,]-prop, err="x", add=TRUE, cex=1, col="dark grey", sfrac=0)
plotCI(x=prop, y=lam, liw=lam-cri.lam[1,], uiw=cri.lam[2,]-lam, err="y", add=TRUE, cex=1, col="dark grey", sfrac=0)
points(prop, lam,pch=16)
axis(1,at=c(0.3,0.5,0.7),tck=0.04, col="dark grey")
axis(2,at=c(0.9,1.0,1.1),labels=c("","",""), tck=0.04, col="dark grey", las=2)
axis(2,at=c(0.85,0.95,1.05),labels=c("","",""),tck=0.02, col="dark grey", las=2)

box(col="dark grey", lwd=1)
text(0.32,1.06,"D", cex=1.2)
