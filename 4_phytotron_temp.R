#Script to import phytotron temp data

{
# Initialization ---------------------------------------------------------------
rm(list = ls()) # Clear environment

library(lubridate)
library(tidyverse)

# Set paths  -----------------------------------------------------------------
mainDir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
data_dir <- paste0(mainDir,"/data/")
raw_data_dir <- paste0(data_dir,"raw/")
fig_dir<-paste0(mainDir,"/figures/")
}

# find file list in "raw" data folder--------
files <- list.files(path=raw_data_dir,pattern="^4_.*\\.txt$")#find files for run 4
files


#loop through each file (if multiple) and tidy data and combine into data_comb---------
for(i in 1:length(files)) {
  
  data <- read.csv(paste0(raw_data_dir,files[i]))
  
  treat <- as.character(names(data[1]))
  treat<-substr(treat,start = 2, stop = 4)
  
  data <- data %>%
    rename(Temp = Celsius.C.,
           Humid=Humidity..rh.)
  
# create treatment label and datestamp and VPD
  data_2 <- data %>%
    mutate(Treat = treat) %>%
    select(-paste0("X",treat),-Serial.Number,-Dew.Point.C.)%>%
    mutate(Datetime = ymd_hms(Time)) %>% #check date, month, year, etc. function matches input file
    mutate(es = 0.611*exp((17.502*Temp)/(240.97+Temp))) %>% 
    #^above constants from book: an intro to environmental biosphysics, 2nd ed, pg 41
    mutate(VPD = es *(1- (Humid/100)))%>%
    select(-Time)

  if(i==1) {data_comb <- data_2 }
  if(i>1) {data_comb <- rbind(data_comb,data_2) }
  
  rm(data, data_2)
  
}

#summary table + sheet
data_summ <- data_comb %>%
  group_by(Treat) %>%
  summarize(MeanTemp = mean(Temp), seTemp = sd(Temp)/sqrt(n()),
            MeanRH = mean(Humid),seRH = sd(Humid)/sqrt(n()),
            MeanVPD = mean(VPD), seVPD = sd(VPD)/sqrt(n()))

#save data_summ----
  saveRDS(data_summ,paste0(data_dir,"4_metdata_summary.rds"))



#View data----
#view temp data
ggplot(data_comb, aes(x=Datetime, y=Temp,col=Treat))+
  geom_point()+
  geom_line()
#ggsave("temps_run4.jpeg", device = "jpeg", path = paste0(fig_dir,"metdata/"), width = 20, height = 20, units = "cm", dpi = 300)


#view humidity data
ggplot(data_comb, aes(x=Datetime, y=Humid,col=Treat))+
  geom_point()+
  geom_line()
#ggsave("humidity_run4.jpeg", device = "jpeg", path = paste0(fig_dir,"metdata/"), width = 20, height = 20, units = "cm", dpi = 300)

#view VPD data
ggplot(data_comb, aes(x=Datetime, y=VPD,col=Treat))+
  geom_point()+
  geom_line()
#ggsave("VPD_run4.jpeg", device = "jpeg", path = paste0(fig_dir,"metdata/"), width = 20, height = 20, units = "cm", dpi = 300)
