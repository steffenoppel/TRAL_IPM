# DATA PREPARATION FOR LTRE ANALYSIS ###
## prepare output data from IPM saved as individual csv files
## convert into a single table with years in rows and cohorts in columns


library(tidyverse)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select

### for most parameters the years are in separate rows
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
parameters <- c("Ntot","Ntot.breed","ann.fec","phi.ad","phi.juv","p.ad","breed.prop","agebeta","mean.p.juv") ## added IM and JUV to facilitate LTRE analysis
selrows<-list(seq(1,18,1),seq(1,18,1),seq(1,18,1),seq(27,44,1),seq(27,44,1),seq(27,44,1),seq(1,18,1),rep(1,18),rep(2,18))
LTRE_input<-data.frame(Year=seq(2004,2021,1))
for(p in 1:length(parameters)){
  input<-fread(sprintf("IPM_output_%s.csv",parameters[p]))
  LTRE_input[,p+1]<-input$Median[selrows[[p]]]
  names(LTRE_input)[p+1]<-parameters[p]
  fwrite(LTRE_input, "LTRE_input_extended.csv")
}

### for immature birds we need to split by age group
LTRE_input<-fread("IPM_output_IM.csv") %>% select(parameter,median) %>%
  mutate(Age= as.numeric(str_match(parameter, "\\,\\s*(.*?)\\s*\\,")[,2])) %>%
  mutate(Year= as.numeric(str_match(parameter, "\\[\\s*(.*?)\\s*\\,")[,2])) %>%
  arrange(Age,Year) %>%
  mutate(Cohort=paste("IM",Age,sep="")) %>%
  select(Cohort,Year, median) %>%
  spread(key=Cohort,value=median) %>%
  mutate(Ntot.IM = rowSums(across(where(is.numeric)))-Year) %>%
  filter(Year<19) %>%
  mutate(Year=Year+2003) %>%
  left_join(LTRE_input, by="Year") %>%
  mutate(Ntot.nonbreed=Ntot-Ntot.breed-Ntot.IM)


fwrite(LTRE_input, "LTRE_input_extended.csv")




#########################################################################
# DO THE ABOVE FOR ALL MCMC SAMPLES
#########################################################################
load("TRAL_IPM_output_REV2022_FINAL.RData")
str(TRALipm$mcmc)
retain<-parameters[c(8,15,4,11,12,14,9,13,18)]

### need the following parameters from model
#which(dimnames(TRALipm$mcmc[[1]])[[2]]=="lambda[17]") # lambda: 8-24
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="phi.ad[43]") # phi.ad: 177-194
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="phi.juv[43]") # phi.juv: 220-237
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="ann.fec[1]") # ann.fec: 256 - 273
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="breed.prop[18]") # breed.prop: 4-21
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot.breed[18]") # Ntot.breed: 238-255
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="Ntot[18]") # Ntot: 44-61
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="IM[1,1,1]") # IM: 321-860
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="IM[18,30,1]") # IM: 321-860
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="agebeta")
which(dimnames(TRALipm$mcmc[[1]])[[2]]=="mean.p.juv[2]")

retain
retaincols<-c(43, #agebeta
              275, #mean.p.juv[2]
              4:21, # breed.prop
              177:194, #phi.ad
              220:238, # phi.juv
              256:273, # ann.fec
              44:61, #Ntot
              238:255, #Ntot.breed
              321:860) # IM year 1:18 for ages 1:30

### EXTRACT ALL FROM THE MODEL OUTPUT AND SAVE IN DIFFERENT LIST

LTRE_input_mcmc<-as.matrix(TRALipm$mcmc[[1]])[,retaincols]
str(LTRE_input_mcmc)

for(ch in 2:nc){
  LTRE_input_mcmc<-rbind(LTRE_input_mcmc,as.matrix(TRALipm$mcmc[[ch]])[,retaincols])
}

rm(list= ls()[!(ls() %in% c('LTRE_input_mcmc','parameters'))])
save.image("TRAL_LTRE_input.RData")

str(LTRE_input_mcmc)
