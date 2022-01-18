# DATA PREPARATION FOR LTRE ANALYSIS ###
## prepare output data from IPM saved as individual csv files
## convert into a single table with years in rows and cohorts in columns


library(tidyverse)
library(data.table)
filter<-dplyr::filter
select<-dplyr::select

### for most parameters the years are in separate rows
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\TRAL_IPM")
parameters <- c("Ntot","Ntot.breed","ann.fec","phi.ad","phi.juv") ## added IM and JUV to facilitate LTRE analysis
LTRE_input<-data.frame(Year=seq(2004,2020,1))
for(p in 1:length(parameters)){
  input<-fread(sprintf("IPM_output_%s.csv",parameters[p]))
  LTRE_input[,p+1]<-input$Median[1:17]
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
  filter(Year<18) %>%
  mutate(Year=Year+2003) %>%
  left_join(LTRE_input, by="Year") %>%
  mutate(Ntot.nonbreed=Ntot-Ntot.breed-Ntot.IM)


fwrite(LTRE_input, "LTRE_input_extended.csv")

