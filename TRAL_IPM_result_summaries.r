##########################################################################
#
# TRISTAN ALBATROSS INTEGRATED POPULATION MODEL 2001-2018
#
##########################################################################
# based on output created in TRAL_IPM_v3.r
# includes JAGS output from 4 scenarios of AYNA population trajectory

# changed on 3 April 2020 to adjust for reduced parameters monitored
# changed on 22 April to incorporate 3 scenarios


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
load("TRAL_IPM_output_v4.RData")




#########################################################################
# PRODUCE OUTPUT TABLES THAT COMBINE ALL 4 SCENARIOS
#########################################################################


## write output into file ##
export<-as.data.frame(TRALipm$summary) %>% select(c(1,5,2,3,7,8)) %>%
  setNames(c('Mean', 'Median','SD','lcl', 'ucl','Rhat')) %>%
  mutate(parameter=row.names(TRALipm$summary)) %>%
  mutate(parameter=ifelse(grepl("ann.surv\\[1,",parameter,perl=T,ignore.case = T)==T,"juv.survival",parameter)) %>%
  mutate(parameter=ifelse(grepl("ann.surv\\[2,",parameter,perl=T,ignore.case = T)==T,"adult.survival",parameter)) %>%
  mutate(Year=c(seq(2001,2020,1),rep(seq(2021,2050,1),each=3),   ## for Ntot - with 3 scenarios
                #seq(2001,2019,1),   ## for lambda
                #seq(2020,2049,1),   ## fut.lambda
                seq(2001,2020,1),   ## ann.fec
                rep(seq(2000.5,2018.5,1),each=2), ##  juv and adult survival
                rep(NA,10))) %>%
  filter(!parameter=="juv.survival") %>%
  arrange(parameter,Year)
tail(export)
hist(export$Rhat)

write.table(export,"TRAL_Gough_IPM_estimates_2020_v5.csv", sep=",", row.names=F)



## write output into file ##
# exportmouse<-as.data.frame(TRALmouse$summary) %>% select(c(1,5,2,3,7,8)) %>%
#   setNames(c('Mean', 'Median','SD','lcl', 'ucl','Rhat')) %>%
#   mutate(parameter=row.names(TRALipm$summary)) %>%
#   mutate(parameter=ifelse(grepl("ann.surv\\[1,",parameter,perl=T,ignore.case = T)==T,"juv.survival",parameter)) %>%
#   mutate(parameter=ifelse(grepl("ann.surv\\[2,",parameter,perl=T,ignore.case = T)==T,"adult.survival",parameter)) %>%
#   mutate(Year=c(seq(2001,2049,1),   ## for Ntot
#                 seq(2001,2018,1),   ## for lambda
#                 seq(2019,2048,1),   ## fut.lambda
#                 seq(2001,2019,1),   ## ann.fec
#                 rep(seq(2000.5,2017.5,1),each=2), ##  juv and adult survival
#                 rep(NA,8))) 
# hist(exportmouse$Rhat)



#########################################################################
# PRODUCE TABLE 1 THAT SUMMARISES DEMOGRAPHIC RATES
#########################################################################

TABLE1<-export %>% mutate(parameter=ifelse(grepl("ann.fec",parameter,perl=T,ignore.case = T)==T,"fecundity",parameter))%>%
    mutate(parameter=ifelse(grepl("imm.rec",parameter,perl=T,ignore.case = T)==T,"recruitment",parameter))%>%
    mutate(parameter=ifelse(grepl("skip.",parameter,perl=T,ignore.case = T)==T,"breed.propensity",parameter))%>%
    mutate(parameter=ifelse(grepl("Ntot.breed",parameter,perl=T,ignore.case = T)==T,"pop.size",parameter)) %>%
  mutate(parameter=ifelse(grepl("lambda",parameter,perl=T,ignore.case = T)==T,"population.growth.rate",parameter)) %>%
  group_by(parameter) %>%
  summarise(median=mean(Median),lcl=mean(lcl),ucl=mean(ucl)) %>%
  filter(!parameter=="deviance") %>%
  filter(!parameter=="pop.size") %>%
  filter(!parameter=="beta[2]") %>%
  mutate(parameter=ifelse(parameter=="beta[1]","juvenile.survival",parameter))

write.table(TABLE1,"TRAL_demographic_estimates_2020.csv", sep=",", row.names=F)






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
  filter(grepl("Ntot.breed",parameter,perl=T,ignore.case = T)) %>%
  arrange(Year) %>%
  mutate(Scenario="status quo") %>%
  mutate(Scenario=if_else(grepl("f\\[2",parameter,perl=T,ignore.case = T), "no mice",if_else(grepl("f\\[3",parameter,perl=T,ignore.case = T),"worse mice",Scenario))) %>%
  mutate(ucl=if_else(ucl>10000,9999,ucl)) %>%
  filter(!(Median<500 & Year<2020)) %>%
  #mutate(Median=ifelse(Year>2021 & Scenario=="status quo",NA,Median)) %>%
  
## CREATE PLOT FOR POP TREND AND SAVE AS PDF

  ggplot() + 
  geom_line(aes(y=Median, x=Year, colour=Scenario), size=1)+   #
  geom_ribbon(aes(x=Year, ymin=lcl,ymax=ucl, fill=Scenario),alpha=0.3)+ #
  geom_hline(aes(yintercept=max(TRAL.pop$tot)),col="darkgray",linetype=2,size=1.5)+
  geom_hline(aes(yintercept=min(TRAL.pop$tot[TRAL.pop$tot>500])),col="darkgray",linetype=2,size=1.5)+
  geom_point(data=TRAL.pop[TRAL.pop$tot>500 & TRAL.pop$tot<2395,],aes(y=tot, x=Year),col="firebrick", size=2.5)+
  geom_point(data=TRAL.pop[TRAL.pop$Year %in% c(2001,2011),],aes(y=tot, x=Year),col="salmon", size=2)+

  ylab("Breeding TRAL pairs on Gough") +
  scale_y_continuous(breaks=seq(0,3000,500), limits=c(0,3000))+
  scale_x_continuous(breaks=seq(2001,2050,5))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

#ggsave("TRAL_IPM_pop_trend_Gough_2001_2050_3scenarios.pdf", width=14, height=8)
ggsave("TRAL_IPM_pop_trend_Gough_2001_2050_3scenarios_v5.jpg", width=14, height=8)




### PLOT SIMPLE LINEAR REGRESSION FOR ANDREW

summary(lm(tot~Year, data=TRAL.pop[TRAL.pop$tot>500 & TRAL.pop$tot<2395,]))
summary(lm(tot~Year, data=TRAL.pop[TRAL.pop$tot>500,]))
  
ggplot(TRAL.pop[TRAL.pop$tot>500,]) + 
  geom_smooth(aes(y=tot, x=Year),method="lm") +
  geom_point(aes(y=tot, x=Year),col="firebrick", size=2.5)+
  #geom_point(data=TRAL.pop[TRAL.pop$Year %in% c(2001,2011),],aes(y=tot, x=Year),col="salmon", size=2)+
  ylab("Breeding TRAL pairs on Gough") +
  scale_y_continuous(breaks=seq(0,3000,500), limits=c(0,3000))+
  scale_x_continuous(breaks=seq(2000,2020,2))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

ggsave("TRAL_2001_2020_linear_trend.jpg", width=14, height=8)









#########################################################################
# PRODUCE OUTPUT GRAPH FOR SURVIVAL AND FECUNDITY
#########################################################################


## CREATE PLOT FOR SURVIVAL AND SAVE AS PDF
pdf("TRAL_IPM_survival_Gough_2001_2019.pdf", width=11, height=8)
export %>% filter(grepl("survival",parameter,perl=T,ignore.case = T)) %>%
  mutate(Year=ifelse(parameter=="adult.survival",Year-0.15,Year+0.15)) %>%  ### to avoid overplotting

  ggplot(aes(y=Median, x=Year, colour=parameter)) + geom_point(size=2.5)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  geom_hline(yintercept=export$Median[export$parameter=="beta[1]"], colour="#619CFF") +
  geom_hline(yintercept=export$Median[export$parameter=="beta[2]"], colour="#619CFF") +
  ylab("Apparent annual survival probability") +
  scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0,1))+
  scale_x_continuous(breaks=seq(2001,2019,2))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()





## CREATE PLOT FOR FECUNDITY AND SAVE AS PDF
pdf("TRAL_IPM_fecundity_Gough_2001_2019.pdf", width=11, height=8)
export %>% filter(grepl("ann.fec",parameter,perl=T,ignore.case = T)) %>%
  
  ggplot(aes(y=Median, x=Year)) + geom_point(size=2.5)+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  geom_hline(yintercept=export$Median[export$parameter=="mean.fec"], colour="#619CFF") +
  ylab("Annual fecundity") +
  scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0,1))+
  scale_x_continuous(breaks=seq(2001,2020,2))+
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CORRELATION OF POP GROWTH RATE WITH DEMOGRAPHIC PARAMETERS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## no longer sensible after revisions in April 2020, due to constant pop growth rate

# ##################  ####################################
# #pdf("TRAL_lambda_correlations.pdf", width=9, height=9)
# l.fitted<-l.lower<-l.upper<-ad.fitted<-ad.lower<-ad.upper<-ju.fitted<-ju.lower<-ju.upper<-pr.fitted<-pr.lower<-pr.upper<-bp.fitted<-bp.lower<-bp.upper<-rec.fitted<-rec.lower<-rec.upper<-numeric()
# year<-c(2001:2019)
# n.years<-length(year)
# for (i in 1:(n.years-1)){
#   # l.fitted[i]<-quantile(TRALipm$sims.list$pop.growth.rate[,i], 0.5)
#   # l.lower[i]<-quantile(TRALipm$sims.list$pop.growth.rate[,i], 0.025)
#   # l.upper[i]<-quantile(TRALipm$sims.list$pop.growth.rate[,i], 0.975)
#   
#   pr.fitted[i]<-quantile(TRALipm$sims.list$ann.fec[,i], 0.5)
#   pr.lower[i]<-quantile(TRALipm$sims.list$ann.fec[,i], 0.025)
#   pr.upper[i]<-quantile(TRALipm$sims.list$ann.fec[,i], 0.975)
#   
#   # bp.fitted[i]<-quantile(TRALipm$sims.list$skip.prob[,i], 0.5)
#   # bp.lower[i]<-quantile(TRALipm$sims.list$skip.prob[,i], 0.025)
#   # bp.upper[i]<-quantile(TRALipm$sims.list$skip.prob[,i], 0.975)
#   
#   # rec.fitted[i]<-quantile(TRALipm$sims.list$imm.rec[,i], 0.5)
#   # rec.lower[i]<-quantile(TRALipm$sims.list$imm.rec[,i], 0.025)
#   # rec.upper[i]<-quantile(TRALipm$sims.list$imm.rec[,i], 0.975)
#   
#   ad.fitted[i]<-quantile(TRALipm$sims.list$ann.surv[,2,i], 0.5)
#   ad.lower[i]<-quantile(TRALipm$sims.list$ann.surv[,2,i], 0.025)
#   ad.upper[i]<-quantile(TRALipm$sims.list$ann.surv[,2,i], 0.975)
#   
#   # ju.fitted[i]<-quantile(TRALipm$sims.list$ann.surv[,1,i], 0.5)
#   # ju.lower[i]<-quantile(TRALipm$sims.list$ann.surv[,1,i], 0.025)
#   # ju.upper[i]<-quantile(TRALipm$sims.list$ann.surv[,1,i], 0.975)
#   
#   }
# 
# 
# par(mfrow=c(2,2))
# 
# plot(l.fitted~ad.fitted, xlim=c(0.5,1), ylim=c(0.5,1.5), xlab="Adult survival probability",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
# segments(ad.lower,l.fitted,ad.upper,l.fitted ,col="gray", lty=1, lwd=0.5)
# segments(ad.fitted,l.lower,ad.fitted,l.upper ,col="gray", lty=1, lwd=0.5)
# test<-cor.test(l.fitted,ad.fitted,alternative = c("two.sided"),method = "spearman")
# text(0.5,0.55, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)
# 
# # plot(l.fitted~ju.fitted, xlim=c(0.5,1), ylim=c(0.5,1.5), xlab="Juvenile survival probability",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
# # segments(ju.lower,l.fitted,ju.upper,l.fitted ,col="gray", lty=1, lwd=0.5)
# # segments(ju.fitted,l.lower,ju.fitted,l.upper ,col="gray", lty=1, lwd=0.5)
# # test<-cor.test(l.fitted,ju.fitted,alternative = c("two.sided"),method = "spearman")
# # text(0.5,0.55, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)
# 
# plot(l.fitted~pr.fitted, xlim=c(0,1), ylim=c(0.5,1.5), xlab="Annual productivity",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
# segments(pr.lower,l.fitted,pr.upper,l.fitted ,col="gray", lty=1, lwd=0.5)
# segments(pr.fitted,l.lower,pr.fitted,l.upper ,col="gray", lty=1, lwd=0.5)
# test<-cor.test(l.fitted,pr.fitted,alternative = c("two.sided"),method = "spearman")
# text(0,0.55, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)
# 
# plot(l.fitted[1:18]~longlineICCAT, xlim=c(-3,3), ylim=c(0.5,1.5), xlab="Fishing effort",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
# segments(longlineICCAT,l.lower[1:18],longlineICCAT,l.upper[1:18],col="gray", lty=1, lwd=0.5)
# test<-cor.test(l.fitted[1:18],longlineICCAT,alternative = c("two.sided"),method = "spearman")
# text(0,0.55, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)
# 
# dev.off()
# 
# 
# 
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # CORRELATION OF DEMOGRAPHIC PARAMETERS WITH FISHING EFFORT
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# 
# par(mfrow=c(2,2))
# 
# plot(longlineICCAT~ad.fitted[1:18], xlim=c(0.5,1), ylim=c(-3,3), xlab="Adult survival probability",ylab="Fishing effort",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
# segments(ad.lower[1:18],longlineICCAT,ad.upper[1:18],longlineICCAT ,col="gray", lty=1, lwd=0.5)
# test<-cor.test(longlineICCAT,ad.fitted[1:18],alternative = c("two.sided"),method = "spearman")
# text(0.5,-2.5, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)
# 
# plot(longlineICCAT~ju.fitted[1:18], xlim=c(0.5,1), ylim=c(-3,3), xlab="Juvenile survival probability",ylab="Fishing effort",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
# segments(ju.lower[1:18],longlineICCAT,ju.upper[1:18],longlineICCAT ,col="gray", lty=1, lwd=0.5)
# test<-cor.test(longlineICCAT,ju.fitted[1:18],alternative = c("two.sided"),method = "spearman")
# text(0.5,-2.5, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)
# 
# plot(longlineICCAT~pr.fitted[1:18], xlim=c(0,1), ylim=c(-3,3), xlab="Annual productivity",ylab="Fishing effort",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
# segments(pr.lower[1:18],longlineICCAT,pr.upper[1:18],longlineICCAT ,col="gray", lty=1, lwd=0.5)
# test<-cor.test(longlineICCAT,pr.fitted[1:18],alternative = c("two.sided"),method = "spearman")
# text(0,-2.5, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)
# 
# plot(longlineICCAT~l.fitted[1:18], ylim=c(-3,3), xlim=c(0.5,1.5), xlab="Population growth rate",ylab="Fishing effort",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
# segments(l.lower[1:18],longlineICCAT,l.upper[1:18],longlineICCAT,col="gray", lty=1, lwd=0.5)
# test<-cor.test(l.fitted[1:18],longlineICCAT,alternative = c("two.sided"),method = "spearman")
# text(0.5,-2.5, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
