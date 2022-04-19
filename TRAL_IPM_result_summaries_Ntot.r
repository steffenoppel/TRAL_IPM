##########################################################################
#
# TRISTAN ALBATROSS INTEGRATED POPULATION MODEL 2001-2018
#
##########################################################################
# based on output created in TRAL_IPM_v3.r
# includes JAGS output from 4 scenarios of AYNA population trajectory

# changed on 3 April 2020 to adjust for reduced parameters monitored
# changed on 22 April to incorporate 3 scenarios

# updated 15 January 2021 to include new output from m-array
# split into new file on 21 January 2021 to include new output for total population size

# updated 9 August 2021 to include final output in different form

# updated 4 Oct 2021 to include Figure S2 (showing no difference in effort formulation)

# updated 21 Dec 2021 to include revision suggestion by referees

### TO DO LIST
# (1) the probability that the given management scenario would result in a larger population size in 2050 relative to the size predicted with no change
# (2) cumulative probability of annual population extinction (or quasi-extinction if a threshold other than zero is of interest) under each scenario
# (3) the probability that future lambda is >1
# (4) the probability that pop in 2050 is greater than pop in 2004
# Plot of correlation between lambda and demographic estimates (phi.ad, phi.juv, p.ad, ann.fec)

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
#load("TRAL_IPM_output_FINAL_REV2021.RData")
load("TRAL_IPM_output_REV2022_FINAL.RData")
imgTRAL<-image_read("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\PR_Comms\\Icons\\alby 4.jpg") %>% image_transparent("white", fuzz=5)
TRALicon <- rasterGrob(imgTRAL, interpolate=TRUE)



#########################################################################
# PRODUCE OUTPUT TABLES THAT COMBINE ALL 3 SCENARIOS
#########################################################################

## based on the parameters defined in the published code

summary_tralipm <- summary(TRALipm, vars=parameters[c(1:3,6,7,9,10,11,12,13,14)])   ## remove sigmas and others we don't need
predictions <- data.frame(summary_tralipm,
                          parameter = row.names(summary_tralipm))
predictions$parameter



## write output into file ##
export<-predictions %>% filter(!grepl("lambda",parameter)) %>%
  #filter(!grepl("Ntot.breed",parameter)) %>%
  #filter(!grepl("agebeta",parameter)) %>%
  mutate(Year=c(
    rep(NA,3),         ## for mean phi and fec
    #seq(2004,2021,1),   ## for breed.prop
    rep(NA,4),         ## for growth rates
    #seq(2004,2020,1), ##  for lambda
    seq(2004,2021,1),   ## for N.tot
    rep(seq(2022,2051,1),each=3), ##  for Ntot.f with 3 scenarios
    seq(2004,2021,1),   ## for Ntot.breed
    rep(seq(1979,2021,1), 2), ##  for phi.ad and phi.juv
    seq(2004,2021,1)   ## for ann.fec
    #rep(NA,4),         ## for mean p
    #seq(1979,2021,1) ##  for p.ad
  )) %>%     ## for deviance and agebeta
  mutate(demographic=parameter) %>%
  mutate(demographic=ifelse(grepl("fec",parameter,perl=T,ignore.case = T)==T,"fecundity",demographic))%>%
  mutate(demographic=ifelse(grepl("phi",parameter,perl=T,ignore.case = T)==T,"survival",demographic))%>%
  mutate(demographic=ifelse(grepl("Ntot",parameter,perl=T,ignore.case = T)==T,"pop.size",demographic)) %>%
  mutate(demographic=ifelse(grepl("growth",parameter,perl=T,ignore.case = T)==T,"growth.rate",demographic)) %>%
  mutate(demographic=ifelse(grepl("agebeta",parameter,perl=T,ignore.case = T)==T,"agebeta",demographic)) %>%
  rename(Rhat=psrf) %>%
  arrange(demographic,Year)
tail(export)

## assess convergence and sample size
summary(export$Rhat)
summary(export$SSeff)

#write.table(export,"TRAL_Gough_IPM_estimates_2022.csv", sep=",", row.names=F)




#########################################################################
# SUMMARIES FOR TEXT
#########################################################################

### change in breeding pop size
jags.data$R[jags.data$R<1000]<-NA
summary(lm(jags.data$R~seq(2004,2021,1)))
range(jags.data$R, na.rm=T)

# bsout<-export %>% filter(grepl("Ntot.breed",parameter)) %>% 
#   select(Year,Median,Lower95,Upper95)
# summary(lm(Median~Year,data=bsout))
# range(bsout$Median, na.rm=T)


##change in breeding success
jags.data$J[is.na(jags.data$R)]<-NA
jags.data$success<-jags.data$J/jags.data$R
#summary(lm(success~Year,data=PROD.DAT))
bsout<-export %>% filter(demographic=="fecundity") %>% filter(!is.na(Year)) %>%
  select(Year,Median,Lower95,Upper95)
summary(lm(Median~Year,data=bsout))
range(bsout$Median, na.rm=T)

### FINAL POPULATION SIZES ###
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot.f[1,30]")


### POPULATION SIZES IN RESULTS SECTION ###
(export %>% filter(parameter %in% c("Ntot[1]","Ntot[18]")) %>% select(Median, Lower95,Upper95)) *2
(export %>% filter(parameter %in% c("Ntot.f[2,30]")) %>% select(Median, Lower95,Upper95)) *2

export %>% filter(demographic=="fecundity") %>% select(Year, Median, Lower95,Upper95) %>% arrange(Median)



#########################################################################
# PRODUCE TABLE 1 THAT SUMMARISES DEMOGRAPHIC RATES
#########################################################################

TABLE1<-export %>% 
  filter(!grepl("Ntot",parameter)) %>%
  filter(parameter %in% c("mean.phi.ad",
                          "mean.phi.juv",
                          "mean.fec",
                          "pop.growth.rate",
                          "fut.growth.rate[1]",
                          "fut.growth.rate[2]",
                          "fut.growth.rate[3]"
                          #"mean.propensity",
                          #"mean.recruit",
                           )) %>%
  mutate(order=c(3,4,5,6,7,1,2)) %>%
  arrange(order) %>%
  mutate(CredInt=paste(round(Lower95,3), "-", round(Upper95,3))) %>%
  select(parameter,Median,CredInt,Lower95,Upper95)

#write.table(TABLE1,"TRAL_demographic_estimates_2022.csv", sep=",", row.names=F)





#########################################################################
# PRODUCE OUTPUT GRAPH FIGURE 2 THAT SHOWS ESTIMATES FOR POPULATION TREND
#########################################################################
## INCLUDED DIFFERENT SCENARIOS ON 22 APRIL 2020
## scenario 1: projection with no changes in demography
## scenario 2: successful mouse eradication in 2021 - fecundity doubles
## scenario 3: increasing mouse impacts on adult survival (adult survival decreases by 10%)



## PREPARE PLOTTING DATAFRAME
plot1_df <- export %>%
  rename(lcl=Lower95,ucl=Upper95) %>%
  filter(grepl("Ntot",parameter,perl=T,ignore.case = T)) %>%
  filter(!grepl("Ntot.breed",parameter)) %>%
  arrange(Year) %>%
  mutate(Scenario="no future change") %>%
  mutate(Scenario=if_else(grepl("f\\[2",parameter,perl=T,ignore.case = T), "after successful mouse eradication",if_else(grepl("f\\[3",parameter,perl=T,ignore.case = T),"no mouse eradication and worsening impacts",Scenario))) %>%
  mutate(ucl=if_else(ucl>15000,15000,ucl)) %>%
  filter(!(Median<500 & Year<2020)) #%>%
#mutate(Median=ifelse(Year>2021 & Scenario=="status quo",NA,Median)) %>%


## CREATE PLOT FOR POP TREND AND SAVE AS PDF

TRAL.pop<-data.frame(Year=seq(2004,2021,1),tot=jags.data$R,line="observed trend")
ggplot(plot1_df) + 
  geom_line(aes(y=Median*2, x=Year, colour=Scenario), size=1)+   #
  geom_ribbon(aes(x=Year, ymin=lcl*2,ymax=ucl*2, fill=Scenario),alpha=0.3)+ #
  #scale_color_manual(values=c('#4393c3','#d6604d','#b2182b')) +
  #scale_fill_manual(values=c('#4393c3','#d6604d','#b2182b')) +
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  
  ## add the breeding pair count data
  geom_point(data=TRAL.pop[TRAL.pop$tot>500 & TRAL.pop$tot<2395,],aes(y=tot*2, x=Year),col="black", size=2.5)+
  geom_smooth(data=TRAL.pop[TRAL.pop$tot>500 & TRAL.pop$tot<2395,],aes(y=tot*2, x=Year, lty=line),method="lm",se=T,col="grey12", size=1)+

  #ylab() +
  #xlab("Year") +
  #geom_vline(aes(xintercept = 2022), colour="gray15", linetype = "dashed", size=1) +
  geom_segment(aes(x = 2022, y = 0, xend = 2022, yend = 14000), colour="gray15", linetype = "dashed", size=1)+
  scale_y_continuous(breaks=seq(0,18000,2000), limits=c(0,20000),expand = c(0, 0))+
  scale_x_continuous(breaks=seq(2005,2050,5), limits=c(2004,2050))+
  #scale_linetype_manual(name="Breeding population",label="observed trend") +
  labs(x="Year", y="Tristan Albatross Population Size (Individuals)",
       col="Total population scenario",
       fill="Total population scenario",
       linetype="Breeding population") +
  
  ### add the bird icons
  annotation_custom(TRALicon, xmin=2045, xmax=2050, ymin=15000, ymax=21000) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        legend.text=element_text(size=14),
        legend.title = element_text(size=16),
        legend.position=c(0.21,0.82),
        panel.grid.major = element_line(size=.1, color="grey94"),
        #panel.grid.major.y = element_line(size=.1, color="grey37"), 
        #panel.grid.major.x = element_blank(), 
        #panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))

#ggsave("TRAL_IPM_pop_trend_Gough_2004_2050_Ntot.jpg", width=14, height=8)
ggsave("C:\\STEFFEN\\MANUSCRIPTS\\Submitted\\TRAL_IPM\\Fig2_revised.jpg", width=14, height=8)


## ALTERNATIVE PLOT WITH TWO SEPARATE AXES
# 
# ggplot(plot1_df) + 
#   geom_line(aes(y=Median*2, x=Year, colour=factor(Scenario)), size=1) +   #
#   geom_ribbon(aes(x=Year, ymin=lcl*2,ymax=ucl*2, fill=Scenario),alpha=0.3)+ #
#   geom_point(data=TRAL.pop[TRAL.pop$tot>500 & TRAL.pop$tot<2395,],aes(y=tot*5, x=Year),col="firebrick", size=2.5)+
#   geom_point(data=TRAL.pop[TRAL.pop$Year %in% c(2001,2011),],aes(y=tot*5, x=Year),col="salmon", size=2)+
#   ylab("\nGlobal Tristan Albatross Population (Individuals)") +
#   xlab("Year") +
#   labs(color = "Scenario", fill = "Scenario")+
#   scale_y_continuous(breaks=seq(0,27000,2000), limits=c(0,27000),
#                      sec.axis = sec_axis(~ . / 5, 
#                                          name = "Number of Breeding Pairs on Gough\n"))+
#   scale_x_continuous(breaks=seq(2001,2050,5))+
#   theme_bw()+
#   theme( 
#     axis.text=element_text(size=16, color="black"), 
#     axis.title=element_text(size=18),
#     legend.text=element_text(size=14),
#     legend.title = element_text(size=16),
#     legend.position=c(0.1,0.9),
#     #panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank()#, 
#     #panel.border = element_blank()
#   )






######################################################################################
# FIGURE 3 SHOWING PROBABILITY OF LAMBDA <1
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

## PROBABILITY OF FUTURE GROWTH
fut.lam.samp %>% gather(key="scenario",value="lam") %>% mutate(pos=ifelse(lam>1,1,0)) %>%
  group_by(scenario) %>%
  summarise(prob=sum(pos)/dim(fut.lam.samp)[1])


## quantify probability of lambda>1
fut.lam.samp %>% rename(nochange=`fut.growth.rate[1]`,erad=`fut.growth.rate[2]`,worse=`fut.growth.rate[3]`) %>%
  mutate(comp12=ifelse(nochange<1,0,1),comp23=ifelse(erad<1,0,1),comp13=ifelse(worse<1,0,1)) %>%
  summarise(prob12=mean(comp12),prob23=mean(comp23),prob13=mean(comp13))

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
  guides(fill = guide_legend(reverse=TRUE)) +
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  scale_x_continuous(breaks=seq(0.9,1.04,0.01), limits=c(0.9,1.04))+
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14),
        legend.position="top", ##c(0.14,0.88),
        panel.grid.major = element_line(size=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))


ggsave("C:\\STEFFEN\\MANUSCRIPTS\\Submitted\\TRAL_IPM\\Fig3_revised.jpg", width=14, height=8)




#########################################################################
# SUMMARY OF POPULATION SIZES
#########################################################################
plot1_df %>% filter(Year==2051) %>% select(Scenario, Median, lcl,ucl) %>% mutate(Median=Median*2, lcl=lcl*2,ucl=ucl*2)
plot1_df %>% filter(Year==2004) %>% select(Scenario, Median, lcl,ucl) %>% mutate(Median=Median*2, lcl=lcl*2,ucl=ucl*2)
plot1_df %>% filter(Year==2021) %>% select(Scenario, Median, lcl,ucl) %>% mutate(Median=Median*2, lcl=lcl*2,ucl=ucl*2)
#fwrite(plot1_df,"TRAL_pop_projections_FINAL.csv")

### FINAL POPULATION SIZES ###
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot.f[2,30]")
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot[1]")
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot[18]")

## collate all samples
# fut.pop.samp<-data.frame()
# for(ch in 1:nc){
#   fut.pop.samp<-bind_rows(fut.pop.samp,as.data.frame((TRALipm$mcmc[[ch]])[,61]))
# }
# median(fut.pop.samp$var1)*2


#########################################################################
# CALCULATE BENEFIT OF ERADICATION
#########################################################################

### MAXIMUM BENEFIT ###
plot1_df %>% filter(Year==2050) %>% 
  mutate(benefit=max(Median)/min(Median),
         benefit.lcl=max(lcl)/min(lcl),
         benefit.ucl=max(ucl)/min(ucl))

### MINIMUM BENEFIT ###
plot1_df %>% filter(Year==2050) %>% 
  mutate(benefit=max(Median)/median(Median),
         benefit.lcl=max(lcl)/median(lcl),
         benefit.ucl=max(ucl)/median(ucl))

### BENEFIT EVEN WITH FAILURE ###
plot1_df %>% filter(Year==2050) %>% 
  mutate(benefit=median(Median)/min(Median),
         benefit.lcl=median(lcl)/min(lcl),
         benefit.ucl=median(ucl)/min(ucl))






#########################################################################
# FIG. S1 - PRODUCE OUTPUT GRAPH THAT SHOWS ESTIMATES FOR PRODUCTIVITY
#########################################################################
bsout<-bsout %>% filter(Year<2022)

## CREATE PLOT FOR POP TREND AND SAVE AS PDF
ggplot(bsout) + 
  geom_point(aes(y=Median, x=Year), size=2, colour="black")+   #
  geom_errorbar(aes(ymin=Lower95, ymax=Upper95, x=Year), width=0.2)+   #
  geom_smooth(aes(y=Median, x=Year),method="lm",se=T,col="grey12", size=1)+
  ylab("Breeding success of Tristan Albatross") +
  xlab("Year") +
  geom_hline(aes(yintercept = 0.63), colour="gray15", linetype = "dashed", size=1) +
  scale_y_continuous(breaks=seq(0,0.8,0.2), limits=c(0,0.8))+
  scale_x_continuous(breaks=seq(2004,2021,2), limits=c(2003.7,2021.2))+
  
  ### add the bird icons
  #annotation_custom(TRALicon, xmin=2045, xmax=2050, ymin=16000, ymax=20000) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        legend.text=element_text(size=14),
        legend.title = element_text(size=16),
        legend.position=c(0.26,0.9),
        #panel.grid.major = element_blank(), 
        #panel.border = element_blank(),
        panel.grid.minor = element_blank())

ggsave("C:\\STEFFEN\\MANUSCRIPTS\\Submitted\\TRAL_IPM\\FigS1_rev.jpg", width=14, height=8)





######################################################################################
# FIG. S2: CALCULATE PROBABILITY OF EXTINCTION AND PROBABILITY OF HIGHER POP IN MGMT SCENARIOS
######################################################################################

# (1) the probability that the given management scenario would result in a larger population size in 2050 relative to the size predicted with no change
# (2) cumulative probability of annual population extinction (or quasi-extinction if a threshold other than zero is of interest) under each scenario
# (3) the probability that future lambda is >1
# (4) the probability that pop in 2050 is greater than pop in 2004

## get all mcmc samples for Ntot.f[1,30] and Ntot.f[2,30] and Ntot.f[3,30]
str(TRALipm)
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot.f[1,30]")
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot.f[2,30]")
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot.f[3,30]")
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot[1]")

## collate all samples
Ntot.f.samp<-data.frame()
for(ch in 1:nc){
  Ntot.f.samp<-bind_rows(Ntot.f.samp,as.data.frame((TRALipm$mcmc[[ch]])[,c(9,114:116)]))
}
head(Ntot.f.samp)
dim(Ntot.f.samp)


## quantify probability of difference
## always 0 because of deterministic nature of projection
Ntot.f.samp %>% rename(past=`Ntot[1]`,nochange=`Ntot.f[1,30]`,erad=`Ntot.f[2,30]`,worse=`Ntot.f[3,30]`) %>%
  mutate(comp12=ifelse(nochange>erad,0,1),comp23=ifelse(worse>erad,0,1),comp13=ifelse(nochange<worse,0,1)) %>%
  mutate(past1=ifelse(nochange<past,0,1),past3=ifelse(worse<past,0,1),past2=ifelse(erad<past,0,1)) %>%
  mutate(ext1=ifelse(nochange<500,1,0),ext2=ifelse(erad<500,1,0),ext3=ifelse(worse<500,1,0)) %>%
  summarise(prob12=mean(comp12),prob23=mean(comp23),prob13=mean(comp13),
            p.past1=mean(past1),p.past2=mean(past2),p.past3=mean(past3),
            p.ext1=mean(ext1),p.ext2=mean(ext2),p.ext3=mean(ext3))

## randomise comparison
Ntot.f.samp %>% rename(past=`Ntot[1]`,nochange=`Ntot.f[1,30]`,erad=`Ntot.f[2,30]`,worse=`Ntot.f[3,30]`) %>%
  mutate(comp12=ifelse(sample(nochange)>sample(erad),0,1),comp23=ifelse(sample(worse)>sample(erad),0,1),comp13=ifelse(sample(nochange)<sample(worse),0,1)) %>%
  mutate(past1=ifelse(nochange<past,0,1),past3=ifelse(worse<past,0,1),past2=ifelse(erad<past,0,1)) %>%
  mutate(ext1=ifelse(nochange<500,1,0),ext2=ifelse(erad<500,1,0),ext3=ifelse(worse<500,1,0)) %>%
  summarise(prob12=mean(comp12),prob23=mean(comp23),prob13=mean(comp13),
            p.past1=mean(past1),p.past2=mean(past2),p.past3=mean(past3),
            p.ext1=mean(ext1),p.ext2=mean(ext2),p.ext3=mean(ext3))

## plot histograms for future pop after 30 years
Ntot.f.samp %>% select(-`Ntot[1]`) %>% rename(nochange=`Ntot.f[1,30]`,erad=`Ntot.f[2,30]`,worse=`Ntot.f[3,30]`) %>%
  gather(key="Scenario",value="N") %>%
  mutate(N = N*2) %>%  ## multiply by 2 to add in males
  mutate(Scenario=if_else(Scenario=="erad", "after successful mouse eradication",if_else(Scenario=="worse","unsuccessful mouse eradication and worse impacts","no future change"))) %>%
  
  ggplot(aes(x = N, fill = Scenario)) +                       # Draw overlaying histogram
  geom_histogram(position = "identity", alpha = 0.2, bins = 80, aes(y = ..density..), color="black") +
  geom_density(alpha=0.5) +
  geom_segment(aes(x = 4981*2, y = 0, xend = 4981*2, yend = 0.001), colour="gray15", linetype = "dashed", size=1)+
  
  scale_y_continuous(breaks=seq(0,0.00125,0.00025), limits=c(0,0.0015),labels=seq(0,0.125,0.025), expand = c(0, 0))+
  labs(x="Tristan Albatross population size in 2050 (Individuals)", y="Probability density",
     fill="Scenario") +
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +

  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        legend.text=element_text(size=14),
        legend.title = element_text(size=16),
        legend.position=c(0.7,0.82),
        panel.grid.major = element_line(size=.1, color="grey94"),
        #panel.grid.major.y = element_line(size=.1, color="grey37"), 
        #panel.grid.major.x = element_blank(), 
        #panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))


ggsave("C:\\STEFFEN\\MANUSCRIPTS\\Submitted\\TRAL_IPM\\FigS2_revised.jpg", width=14, height=8)





#############################################################################################
# PRODUCE PLOT OF ANNUAL LAMBDA AGAINST BREEDING SUCCESS, SURVIVAL, AND RESIGHT PROBABILITY
#############################################################################################

str(TRALipm$mcmc)
#retain<-parameters[c(8,15,4,11,12,14,9,13,18)]


############ PLOT RAW CORRELATIONS BETWEEN DEMOGRAPHIC PARAMETERS #####
### need the following parameters from model
## # Step 1: Using the JAGS output (named TRALipm$mcmc),         # compute realized population growth rates for Tristan Albatross 
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="lambda[1]") # lambda: 26-42
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="phi.ad[43]") # phi.ad: 177-194
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="phi.juv[43]") # phi.juv: 220-237
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="ann.fec[1]") # ann.fec: 256 - 273
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="breed.prop[18]") # breed.prop: 4-21
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot.breed[18]") # Ntot.breed: 238-255
library(matrixStats)
lam <- as.matrix(TRALipm$mcmc[[1]][,c(26:42)])
Sa <- as.matrix(TRALipm$mcmc[[1]][,c(177:194)])
Sj <- as.matrix(TRALipm$mcmc[[1]][,c(220:237)])
Fec <- as.matrix(TRALipm$mcmc[[1]][,c(256:273)])
Bp <- as.matrix(TRALipm$mcmc[[1]][,c(238:255)])
breed.prop <- as.matrix(TRALipm$mcmc[[1]][,c(4:21)])


# Account for different start of time series for survival and fecundity data
for(ch in 2:nc){
  lam<-rbind(lam,as.matrix(TRALipm$mcmc[[ch]][,c(26:42)]))
  Sa<-rbind(Sa,as.matrix(TRALipm$mcmc[[ch]][,c(177:194)]))  ## only for the years coinciding with the fecundity and pop size
  Sj<-rbind(Sj,as.matrix(TRALipm$mcmc[[ch]][,c(220:237)]))  ## only for the years coinciding with the fecundity and pop size
  Fec<-rbind(Fec,as.matrix(TRALipm$mcmc[[ch]][,c(256:273)]))  ## without the last year for which no lambda is available
  Bp<-rbind(Bp,as.matrix(TRALipm$mcmc[[ch]][,c(238:255)]))  ## only for the years coinciding with the fecundity and pop size
  breed.prop<-rbind(breed.prop,as.matrix(TRALipm$mcmc[[ch]][,c(4:21)]))  ## only for the years coinciding with the fecundity and pop size
}


dim(lam)
corplot.l<-data.frame(lambda=colMeans(lam),lcl.lam=apply(lam,2,quantile,probs = 0.025, na.rm = TRUE),ucl.lam=apply(lam,2,quantile,probs = 0.975, na.rm = TRUE),Year=c(2004:2020))
corplot.f<-data.frame(parm=dimnames(Fec)[[2]],value=colMeans(Fec),Year=c(2004:2021),dem="Breeding success",lcl=apply(Fec,2,quantile,probs = 0.025, na.rm = TRUE),ucl=apply(Fec,2,quantile,probs = 0.975, na.rm = TRUE))
corplot.Sa<-data.frame(parm=dimnames(Sa)[[2]],value=colMeans(Sa),Year=c(2004:2021),dem="Adult survival",lcl=apply(Sa,2,quantile,probs = 0.025, na.rm = TRUE),ucl=apply(Sa,2,quantile,probs = 0.975, na.rm = TRUE))
corplot.Sj<-data.frame(parm=dimnames(Sj)[[2]],value=colMeans(Sj),Year=c(2004:2021),dem="Juvenile survival",lcl=apply(Sj,2,quantile,probs = 0.025, na.rm = TRUE),ucl=apply(Sj,2,quantile,probs = 0.975, na.rm = TRUE))
corplot.Bp<-data.frame(parm=dimnames(Bp)[[2]],value=colMeans(Bp),Year=c(2004:2021),dem="N breeding pairs",lcl=apply(Bp,2,quantile,probs = 0.025, na.rm = TRUE),ucl=apply(Bp,2,quantile,probs = 0.975, na.rm = TRUE))
corplot.breed.prop<-data.frame(parm=dimnames(breed.prop)[[2]],value=colMeans(breed.prop),Year=c(2004:2021),dem="Breeding propensity",lcl=apply(breed.prop,2,quantile,probs = 0.025, na.rm = TRUE),ucl=apply(breed.prop,2,quantile,probs = 0.975, na.rm = TRUE))


corplot<-bind_rows(corplot.f,corplot.Sa,corplot.Sj,corplot.Bp,corplot.breed.prop) %>% filter(Year<2021) %>%
  left_join(corplot.l, by="Year")

### Pearson correlation tests
test.vals<-corplot %>% filter(dem!="Breeding propensity") %>% group_by(dem) %>% summarise(x=min(lcl)) %>% mutate(text=NA)
ct.f<-cor.test(corplot.l$lambda,corplot.f$value[1:17])
test.vals$text[2]<-paste("r = ",round(ct.f$estimate,2)," (",round(ct.f$conf.int[1],2),", ",round(ct.f$conf.int[2],2),")", sep="")

ct.Sa<-cor.test(corplot.l$lambda,corplot.Sa$value[1:17])
test.vals$text[1]<-paste("r = ",round(ct.Sa$estimate,2)," (",round(ct.Sa$conf.int[1],2),", ",round(ct.Sa$conf.int[2],2),")", sep="")

ct.Sj<-cor.test(corplot.l$lambda,corplot.Sj$value[1:17])
test.vals$text[3]<-paste("r = ",round(ct.Sj$estimate,2)," (",round(ct.Sj$conf.int[1],2),", ",round(ct.Sj$conf.int[2],2),")", sep="")

ct.Bp<-cor.test(corplot.l$lambda,corplot.Bp$value[corplot.Bp$Year<2021])
test.vals$text[4]<-paste("r = ",round(ct.Bp$estimate,2)," (",round(ct.Bp$conf.int[1],2),", ",round(ct.Bp$conf.int[2],2),")", sep="")

ct.breed.prop<-cor.test(corplot.l$lambda,corplot.breed.prop$value[corplot.breed.prop$Year<2021])
#test.vals$text[5]<-paste("r = ",round(ct.breed.prop$estimate,2)," (",round(ct.breed.prop$conf.int[1],2),", ",round(ct.breed.prop$conf.int[2],2),")", sep="")


## create PLOT
corplot %>% filter(dem!="Breeding propensity") %>%
ggplot() + geom_point(aes(y=lambda,x=value)) +
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

#ggsave("C:\\STEFFEN\\MANUSCRIPTS\\Submitted\\TRAL_IPM\\FigS3_rev.jpg", width=14, height=8)









#########################################################################
# PRODUCE SUPPLEMENTARY FIGURE COMPARING CONT EFFORT AND TWO-INTERCEPT MODEL
#########################################################################

### load the model output from constant effort predictions
load("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\oldTRAL_IPM\\TRAL_IPM_output_FINAL_conteffort.RData")


## write output into file ##
exportconteff<-predictions %>%
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

TABLE1coneff<-exportconteff %>% 
  filter(!grepl("Ntot",parameter)) %>%
  filter(parameter %in% c("fut.growth.rate[1]",
                          "fut.growth.rate[2]",
                          "fut.growth.rate[3]",
                          "mean.fec",
                          "pop.growth.rate",
                          "mean.phi.ad",
                          "mean.phi.juv" )) %>%
  mutate(Model="continuous observation effort") 


FIGS2DATA<-TABLE1 %>%
  mutate(Model="categorical observation effort") %>%
  rename(lcl=Lower95, ucl=Upper95) %>%
  bind_rows(TABLE1coneff) %>%
  select(Model, parameter,Median, lcl,ucl) %>%
  mutate(plotorder=rep(c(3,4,5,6,7,1,2),2)) %>%
  mutate(plotorder=if_else(Model=="categorical observation effort",plotorder-0.2, plotorder+0.2)) %>%
  arrange(plotorder) %>%
  mutate(parameter=rep(c("adult survival",
                         "juvenile survival",
                         "breeding success",
                         "pop.growth (past)",
                         "pop.growth (future) - no change",
                         "pop.growth (future) - eradication",
                         "pop.growth (future) - worse mice"),each=2))

## CREATE PLOT FOR COMPARISON OF PARAMETER ESTIMATES
ggplot(FIGS2DATA) + 
  geom_point(aes(y=Median, x=plotorder, colour=Model), size=1)+   #
  geom_errorbar(aes(ymin=lcl, ymax=ucl, x=plotorder, colour=Model), width=0.1)+   #
  ylab("\nParameter estimate\n") +
  xlab("Integrated population model parameter") +
  scale_y_continuous(breaks=seq(0.2,1.1,0.1), limits=c(0.2,1.1))+
  scale_x_continuous(breaks=seq(1,7,1), limits=c(0,8), labels=FIGS2DATA$parameter[seq(1,13,2)])+
  
  ### add the bird icons
  #annotation_custom(TRALicon, xmin=2045, xmax=2050, ymin=16000, ymax=20000) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=18, color="black"), 
        axis.text.x=element_text(size=14, color="black", angle=45,hjust = 1),        
        axis.title=element_text(size=20),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),
        legend.key=element_blank(),
        #legend.box.background =  = element_blank(),
        legend.position=c(0.80,0.20),
        #panel.grid.major = element_blank(), 
        #panel.border = element_blank(),
        panel.grid.minor = element_blank())

ggsave("C:\\STEFFEN\\MANUSCRIPTS\\Submitted\\TRAL_IPM\\FigS2.jpg", width=14, height=8)

