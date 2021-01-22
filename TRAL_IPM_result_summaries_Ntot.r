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
  #mutate(parameter=ifelse(grepl("ann.surv\\[1,",parameter,perl=T,ignore.case = T)==T,"juv.survival",parameter)) %>%
  #mutate(parameter=ifelse(grepl("ann.surv\\[2,",parameter,perl=T,ignore.case = T)==T,"adult.survival",parameter)) %>%
  # mutate(Year=c(seq(2004,2020,1),rep(seq(2021,2050,1),each=3),   ## for Ntot - with 3 scenarios
  #               #seq(2001,2019,1),   ## for lambda
  #               #seq(2020,2049,1),   ## fut.lambda
  #               seq(2004,2020,1),   ## ann.fec
  #               rep(seq(1979.5,2020.5,1),each=2), ##  juv and adult survival
  #               rep(NA,10))) %>%
  mutate(Year=c(seq(2004,2020,1),   ## for ann.fec
                rep(NA,9),         ## for mean phi, p, and growth rates
                seq(2004,2020,1),   ## for N.tot
                rep(seq(2021,2050,1),each=3), ##  for Ntot.f with 3 scenarios
                seq(2004,2019,1), ##  for lambda 
                seq(1978,2020,1), ##  for p.ad 
                rep(NA,1))) %>%     ## for deviance
  mutate(demographic=parameter) %>%
  mutate(demographic=ifelse(grepl("ann.fec",parameter,perl=T,ignore.case = T)==T,"fecundity",demographic))%>%
  mutate(demographic=ifelse(grepl("p.ad",parameter,perl=T,ignore.case = T)==T,"breed.propensity",demographic))%>%
  mutate(demographic=ifelse(grepl("Ntot",parameter,perl=T,ignore.case = T)==T,"pop.size",demographic)) %>%
  mutate(demographic=ifelse(grepl("lambda",parameter,perl=T,ignore.case = T)==T,"growth.rate",demographic)) %>%
  arrange(demographic,Year)
tail(export)
hist(export$Rhat)

write.table(export,"TRAL_Gough_IPM_estimates_2020_v5.csv", sep=",", row.names=F)




#########################################################################
# PRODUCE TABLE 1 THAT SUMMARISES DEMOGRAPHIC RATES
#########################################################################

TABLE1<-export %>% mutate(parameter=ifelse(grepl("ann.fec",parameter,perl=T,ignore.case = T)==T,"fecundity",parameter))%>%
    mutate(parameter=ifelse(grepl("mean.p.juv",parameter,perl=T,ignore.case = T)==T,"recruitment",parameter))%>%
    mutate(parameter=ifelse(grepl("mean.p.ad",parameter,perl=T,ignore.case = T)==T,"breed.propensity",parameter))%>%
    mutate(parameter=ifelse(grepl("Ntot",parameter,perl=T,ignore.case = T)==T,"pop.size",parameter)) %>%
  group_by(parameter) %>%
  summarise(median=mean(Median),lcl=mean(lcl),ucl=mean(ucl)) %>%
  filter(!parameter=="deviance") %>%
  filter(!parameter=="pop.size") %>%
  filter(!parameter=="mean.fec") %>%
  filter(!grepl("lambda",parameter,perl=T,ignore.case = T)==T)

write.table(TABLE1,"TRAL_demographic_estimates_2020.csv", sep=",", row.names=F)




#########################################################################
# PLOT ADULT RETURN AND RESIGHT PROBABILITY
#########################################################################

export %>% filter(grepl("p.ad",parameter,perl=T,ignore.case = T)==T) %>%
  filter(!parameter=="mean.p.ad") %>%
  arrange(Year) %>%
  
  ggplot(aes(x=Year,y=Median)) + geom_point(size=2) + geom_smooth(method='lm')  + theme_bw()




#########################################################################
# PRODUCE OUTPUT GRAPH THAT SHOWS ESTIMATES FOR POPULATION TREND
#########################################################################
## INCLUDED DIFFERENT SCENARIOS ON 22 APRIL 2020
## scenario 1: projection with no changes in demography
## scenario 2: successful mouse eradication in 2021 - fecundity doubles
## scenario 3: increasing mouse impacts on adult survival (adult survival decreases by 10%)

## CHANGED ON 4 MAY 2020 to address Andrew Callender's concerns


## PLOT THE BASELINE SCENARIO
export %>%
  filter(grepl("Ntot",parameter,perl=T,ignore.case = T)) %>%
  arrange(Year) %>%
  mutate(Scenario="status quo") %>%
  mutate(Scenario=if_else(grepl("f\\[2",parameter,perl=T,ignore.case = T), "no mice",if_else(grepl("f\\[3",parameter,perl=T,ignore.case = T),"worse mice",Scenario))) %>%
  mutate(ucl=if_else(ucl>15000,15000,ucl)) %>%
  filter(!(Median<500 & Year<2020)) %>%
  #mutate(Median=ifelse(Year>2021 & Scenario=="status quo",NA,Median)) %>%
  
## CREATE PLOT FOR POP TREND AND SAVE AS PDF

  ggplot() + 
  geom_line(aes(y=Median*2, x=Year, colour=Scenario), size=1)+   #
  geom_ribbon(aes(x=Year, ymin=lcl*2,ymax=ucl*2, fill=Scenario),alpha=0.3)+ #
  geom_point(data=TRAL.pop[TRAL.pop$tot>500 & TRAL.pop$tot<2395,],aes(y=tot*5, x=Year),col="firebrick", size=2.5)+
  geom_point(data=TRAL.pop[TRAL.pop$Year %in% c(2001,2011),],aes(y=tot*5, x=Year),col="salmon", size=2)+

  ylab("Global Tristan Albatross population (individuals)") +
  scale_y_continuous(breaks=seq(0,20000,2000), limits=c(0,20000),
                     sec.axis = sec_axis(~ . / 5, name = "breeding pairs (on Gough)"))+
  scale_x_continuous(breaks=seq(2001,2050,5))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        legend.position=c(0.1,0.9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

#ggsave("TRAL_IPM_pop_trend_Gough_2001_2050_3scenarios.pdf", width=14, height=8)
ggsave("TRAL_IPM_pop_trend_Gough_2004_2050_Ntot.jpg", width=14, height=8)


