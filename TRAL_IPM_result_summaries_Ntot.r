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
imgTRAL<-image_read("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\PR_Comms\\Icons\\alby 4.jpg") %>% image_transparent("white", fuzz=5)
TRALicon <- rasterGrob(imgTRAL, interpolate=TRUE)



#########################################################################
# PRODUCE OUTPUT TABLES THAT COMBINE ALL 3 SCENARIOS
#########################################################################
### predictions created in TRAL_IPM_FINAL.r

## write output into file ##
export<-predictions %>%
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
  rename(Rhat=psrf) %>%
  arrange(demographic,Year)
tail(export)
hist(export$Rhat)
hist(export$SSeff)

#write.table(export,"TRAL_Gough_IPM_estimates_2021_FINAL.csv", sep=",", row.names=F)




#########################################################################
# SUMMARIES FOR TEXT
#########################################################################
## NEED TO DO: base these regressions on IPM estimates

### change in breeding pop size
PROD.DAT$R[PROD.DAT$R<1000]<-NA
summary(lm(R~Year,data=PROD.DAT))
range(PROD.DAT$R, na.rm=T)

##change in breeding success
PROD.DAT$J[is.na(PROD.DAT$R)]<-NA
PROD.DAT$success<-PROD.DAT$J/PROD.DAT$R
#summary(lm(success~Year,data=PROD.DAT))
bsout<-export %>% filter(demographic=="fecundity") %>% filter(!is.na(Year)) %>%
  select(Year,Median,Lower95,Upper95)
summary(lm(Median~Year,data=bsout))
range(bsout$Median, na.rm=T)


#########################################################################
# PRODUCE TABLE 1 THAT SUMMARISES DEMOGRAPHIC RATES
#########################################################################

TABLE1<-export %>% 
  filter(!grepl("Ntot",parameter)) %>%
  filter(parameter %in% c("fut.growth.rate[1]",
                          "fut.growth.rate[2]",
                          "fut.growth.rate[3]",
                          "mean.fec",
                          #"mean.propensity",
                          #"mean.recruit",
                          "pop.growth.rate",
                          "mean.phi.ad",
                          "mean.phi.juv" )) 

#write.table(TABLE1,"TRAL_demographic_estimates_2021.csv", sep=",", row.names=F)





#########################################################################
# PRODUCE OUTPUT GRAPH THAT SHOWS ESTIMATES FOR POPULATION TREND
#########################################################################
## INCLUDED DIFFERENT SCENARIOS ON 22 APRIL 2020
## scenario 1: projection with no changes in demography
## scenario 2: successful mouse eradication in 2021 - fecundity doubles
## scenario 3: increasing mouse impacts on adult survival (adult survival decreases by 10%)



## PREPARE PLOTTING DATAFRAME
plot1_df <- export %>%
  rename(lcl=Lower95,ucl=Upper95) %>%
  filter(grepl("Ntot",parameter,perl=T,ignore.case = T)) %>%
  arrange(Year) %>%
  mutate(Scenario="past, and no future change") %>%
  mutate(Scenario=if_else(grepl("f\\[2",parameter,perl=T,ignore.case = T), "after successful mouse eradication",if_else(grepl("f\\[3",parameter,perl=T,ignore.case = T),"unsuccessful mouse eradication and worsening impacts",Scenario))) %>%
  mutate(ucl=if_else(ucl>15000,15000,ucl)) %>%
  filter(!(Median<500 & Year<2020)) #%>%
#mutate(Median=ifelse(Year>2021 & Scenario=="status quo",NA,Median)) %>%


## CREATE PLOT FOR POP TREND AND SAVE AS PDF
TRAL.pop$line="observed trend"
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
  scale_y_continuous(breaks=seq(0,18000,2000), limits=c(0,20000),expand = c(0, 0))+
  scale_x_continuous(breaks=seq(2005,2050,5), limits=c(2004,2050))+
  #scale_linetype_manual(name="Breeding population",label="observed trend") +
  labs(x="Year", y="\nTristan Albatross Population Size (Individuals)\n",
       col="Total population scenario",
       fill="Total population scenario",
       linetype="Breeding population") +
  
  ### add the bird icons
  annotation_custom(TRALicon, xmin=2045, xmax=2050, ymin=16000, ymax=20000) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        legend.text=element_text(size=14),
        legend.title = element_text(size=16),
        legend.position=c(0.26,0.82),
        panel.grid.major = element_line(size=.1, color="grey94"),
        #panel.grid.major.y = element_line(size=.1, color="grey37"), 
        #panel.grid.major.x = element_blank(), 
        #panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))

ggsave("TRAL_IPM_pop_trend_Gough_2004_2050_Ntot.jpg", width=14, height=8)
ggsave("C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\TRAL_IPM\\Fig1.jpg", width=14, height=8)




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



#########################################################################
# PRODUCE OUTPUT GRAPH THAT SHOWS ESTIMATES FOR PRODUCTIVITY
#########################################################################
bsout<-bsout %>% filter(Year<2021)

## CREATE PLOT FOR POP TREND AND SAVE AS PDF
ggplot(bsout) + 
  geom_point(aes(y=Median, x=Year), size=2, colour="firebrick")+   #
  geom_errorbar(aes(ymin=Lower95, ymax=Upper95, x=Year), width=0.2)+   #
  geom_smooth(aes(y=Median, x=Year),method="lm",se=T,col="grey12", size=1)+
  ylab("\nBreeding success of Tristan Albatross\n") +
  xlab("Year") +
  scale_y_continuous(breaks=seq(0,0.8,0.2), limits=c(0,0.8))+
  scale_x_continuous(breaks=seq(2004,2020,2), limits=c(2004,2020))+
  
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

ggsave("C:\\STEFFEN\\MANUSCRIPTS\\in_prep\\TRAL_IPM\\FigS2.jpg", width=14, height=8)

