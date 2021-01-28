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
# adjusted to final monitored parameters for presentation in a manuscript - survival over time


library(tidyverse)
library(jagsUI)
library(data.table)
#library(nimble)
filter<-dplyr::filter
select<-dplyr::select


#########################################################################
# LOAD MODEL OUTPUT FROM IPMs
#########################################################################

setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
#load("TRAL_IPM_output_2020.RData")
load("TRAL_IPM_output_v5_Ntot.RData")




#########################################################################
# PRODUCE OUTPUT TABLES THAT COMBINE ALL 3 SCENARIOS
#########################################################################


## write output into file ##
export<-as.data.frame(TRALipm$summary) %>% select(c(1,5,2,3,7,8)) %>%
  setNames(c('Mean', 'Median','SD','lcl', 'ucl','Rhat')) %>%
  mutate(parameter=row.names(TRALipm$summary)) %>%
  mutate(Year=c(#seq(2004,2020,1),   ## for ann.fec
                rep(NA,9),         ## for mean phi, p, and growth rates
                seq(2004,2020,1),   ## for N.tot
                rep(seq(2021,2050,1),each=3), ##  for Ntot.f with 3 scenarios
                rep(seq(1978,2019,1),2), ##  for phi.ad and phi.juv
                rep(NA,1))) %>%     ## for deviance
  mutate(demographic=parameter) %>%
  mutate(demographic=ifelse(grepl("Ntot",parameter,perl=T,ignore.case = T)==T,"pop.size",demographic)) %>%
  mutate(demographic=ifelse(grepl("phi.juv",parameter,perl=T,ignore.case = T)==T,"surv.juv",demographic)) %>%
  mutate(demographic=ifelse(grepl("phi.ad",parameter,perl=T,ignore.case = T)==T,"surv.ad",demographic)) %>%
  arrange(demographic,Year)
tail(export)
hist(export$Rhat)

write.table(export,"TRAL_Gough_IPM_estimates_2020_v5.csv", sep=",", row.names=F)




#########################################################################
# PRODUCE TABLE 1 THAT SUMMARISES DEMOGRAPHIC RATES
#########################################################################

TABLE1<-export %>% 
    filter(is.na(Year)==T) %>%
    filter(!parameter=="deviance") %>%
    select(-Year,-demographic)

write.table(TABLE1,"TRAL_demographic_estimates_2020.csv", sep=",", row.names=F)







#########################################################################
# FIGURE 1: ADULT AND JUVENILE SURVIVAL PROBABILITY
#########################################################################
# LOAD ICON
library(grid)
library(magick)
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\PR_Comms\\Icons")
imgTRAL<-image_read("alby 4.jpg") %>% image_transparent("white", fuzz=5)
TRALicon <- rasterGrob(imgTRAL, interpolate=TRUE)



### subset the data and create plot ##

export %>% filter(grepl("phi",parameter,perl=T,ignore.case = T)==T) %>%
  filter(!parameter=="mean.phi.ad") %>%
  filter(!parameter=="mean.phi.juv") %>%
  arrange(demographic,Year) %>%
  mutate(Age=rep(c("adult","juvenile"), each=42)) %>%

  ggplot(aes(y=Median, x=Year, col=Age)) + geom_point(size=2.5, position=position_dodge(width=0.35))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl, col=Age), width=.1, position=position_dodge(width=0.35))+
  geom_hline(data=TABLE1[TABLE1$parameter=="mean.phi.juv",],aes(yintercept=Median), color='#00BFC4', size=2)+
  geom_hline(data=TABLE1[TABLE1$parameter=="mean.phi.ad",],aes(yintercept=Median), color='#F8766D', size=2)+
  ylab("Apparent annual survival probability") +
  scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0,1))+
  scale_x_continuous(breaks=seq(1980,2020,5))+
  
  ### add the bird icons
  annotation_custom(TRALicon, xmin=1980, xmax=1985, ymin=0, ymax=0.2) +
  
  theme(panel.background=element_rect(fill="white", colour="black"),
        axis.text=element_text(size=18, color="black"),
        axis.title=element_text(size=20),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=16, color="black"),
        legend.background = element_rect(),
        legend.key = element_blank(),
        legend.position=c(0.9,0.1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
ggsave("FIG_1_TRAL_annual_survival_probability.pdf", width=8, height=6)











#########################################################################
# FIGURE 2 - POPULATION TRAJECTORY
#########################################################################

## SUBSET DATA AND CREATE PLOT
export %>%
  filter(grepl("Ntot",parameter,perl=T,ignore.case = T)) %>%
  arrange(Year) %>%
  mutate(Scenario="status quo") %>%
  mutate(Scenario=if_else(grepl("f\\[2",parameter,perl=T,ignore.case = T), "no mice",if_else(grepl("f\\[3",parameter,perl=T,ignore.case = T),"worse mice",Scenario))) %>%
  mutate(ucl=if_else(ucl>15000,15000,ucl)) %>%
  filter(!(Median<500 & Year<2020)) %>%

  ggplot() + 
  geom_line(aes(y=Median*2, x=Year, colour=Scenario), size=1)+   #
  geom_ribbon(aes(x=Year, ymin=lcl*2,ymax=ucl*2, fill=Scenario),alpha=0.3)+ #
  geom_point(data=TRAL.pop[TRAL.pop$tot>500 & TRAL.pop$tot<2395,],aes(y=tot*5, x=Year),col="firebrick", size=2.5)+
  geom_point(data=TRAL.pop[TRAL.pop$Year %in% c(2001,2011),],aes(y=tot*5, x=Year),col="salmon", size=2)+
  
  ### add the bird icons
  annotation_custom(TRALicon, xmin=2045, xmax=2050, ymin=18000, ymax=20000) +
  

  ylab("Global Tristan Albatross population (individuals)") +
  scale_y_continuous(breaks=seq(0,20000,2000), limits=c(0,20000),
                     sec.axis = sec_axis(~ . / 5, name = "breeding pairs (on Gough)"))+
  scale_x_continuous(breaks=seq(2001,2050,5))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=16, color="black"),
        legend.background = element_rect(),
        legend.key = element_blank(),
        legend.position=c(0.1,0.9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

#ggsave("TRAL_IPM_pop_trend_Gough_2001_2050_3scenarios.pdf", width=14, height=8)
ggsave("FIG2_TRAL_IPM_pop_trend_Gough_2004_2050_Ntot.jpg", width=14, height=8)


