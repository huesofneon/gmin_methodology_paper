# create mass_data for purposes of calculating gmin

{ 
# Initialization ---------------------------------------------------------------
rm(list=ls())

library(tidyverse)
library(readxl)
library(lubridate)
library(dplyr)
library(ggpubr)
library(minpack.lm)
library(gridExtra)
library(multcompView)

# Set paths  -----------------------------------------------------------------
mainDir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
data_dir <- paste0(mainDir,"/data/")
raw_data_dir <- paste0(data_dir,"raw/")
fig_dir<-paste0(mainDir,"/figures/")

#formula for standard error
se <- function(x) {
  return(sd(x) / sqrt(length(x)))
}
}

# Part 1: tidy mass loss data----
# find file list in "raw" data folder
excel_files <- list.files(path=raw_data_dir,pattern="^v.*\\.xlsx$")
excel_files


#identify which file you want below (1 or 2)
run <- 2

#Clean up Run 4 Data-----
# List the sheets in the Excel file
sheet_names <- excel_sheets(paste0(raw_data_dir,excel_files[run]))
sheet_names
# Limit to required Temps
sheet_names <-sheet_names[1]

#create df
{
for (i in 1:length(sheet_names)) {
data <- read_xlsx(paste0(raw_data_dir,excel_files[run]), sheet = sheet_names[i])%>%
  na.omit()%>%
  group_by(ID_Code)%>%
  mutate(seq=as.factor(row_number()),Sample_no=as.factor(Sample_no))%>%
  arrange(ID_Code, Date_Time) %>%  # Optional: Sort the data by Species and Date_Time
  mutate(Drydown_time = as.numeric(difftime(Date_Time, first(Date_Time), units = "mins")))

if(i==1) data_combined <- data
else data_combined <- rbind(data_combined, data)
}
#remove NAs
data_combined <- na.omit(data_combined) %>%
  mutate(
    Weight = as.numeric(Weight),
    Drydown_time = as.numeric(Drydown_time)
  )
}

#remove incomplete/failed samples----
data_combined <- data_combined %>%
  filter(ID_Code != "DKF250911SELDOL_20C_sample_5")

#PERCENT MASS LOSS#####
data_percent<-data_combined %>%
  group_by(ID_Code) %>%
  mutate(percent_weight = Weight/max(Weight)) 

#buffer

#View data (uncomment to run)-------
# #Combined plot (wt/drydown time)###
# ggplot(data_combined,aes(x=Drydown_time, y=Weight, col=Sample_no)) +
#   geom_line()+
#   #geom_point()+
#   facet_wrap(~GENSPP, nrow=2, scales = "free")+
#   labs(y="Weight (g)",x="Drydown time (mins)")+
#   scale_y_continuous(sec.axis = sec_axis(~ ., name = "Temperature (C)", breaks =NULL))
# 
# #ggsave("wt_by_drydowntime_run2.jpeg", device = "jpeg", 
# #path = fig_dir, width = 40, height = 20, units = "cm", dpi = 300)
# 
# #Individual species plots (wt/drydown time)###
# spp_list <- unique(data_combined$Species)
# 
# test<-data_combined[data_combined$Species==spp_list[3],]
# 
# ggplot(test,aes(x=Drydown_time, y=Weight, col=Sample_no)) +
#   geom_line()+
#   geom_point()+
#   facet_wrap(~Treat,nrow = 2)+
#   labs(title = test$Species)+
#   theme(plot.title = element_text(hjust=0.5))
# 
# #Combined plot (%/t)###
# ggplot(data_percent,aes(x=Drydown_time, y=percent_weight, col=Sample_no)) +
#   geom_line()+
#   geom_point()+
#   facet_wrap(~GENSPP, nrow=2, scales="free")+
#   labs(y="Percent weight",x="Drydown time (mins)")+
#   scale_y_continuous(sec.axis = sec_axis(~ ., name = "Temperature (C)", breaks =NULL))
# 
# #Individual plots (%/t)###
# GENSPP_list <- unique(data_percent$GENSPP)
# 
# test<-data_percent[data_percent$GENSPP==GENSPP_list[1],]
# #select treatment temp
# test<-test%>%
#   filter(Treat=="20C")
# 
# ggplot(test,aes(x=Drydown_time, y=percent_weight, col=Sample_no)) +
#   geom_line()+
#   facet_wrap(~Treat,nrow = 2)+
#   labs(title = test$GENSPP)+
#   theme(plot.title = element_text(hjust=0.5))

# Part 2: combine raw data----
{
#read in data-----
#define atmospheric pressure at height above sea level 
atmosp <- 99891.70/1000#kPa (@120m = approx height of HW Pearson)

#massloss data
mass_data <- data_percent
#read in dryweights for RWC
dry_wt <- read_xlsx(paste0(raw_data_dir,"dryweights_run4.xlsx"))%>%
  select(-total_wt,-Sample,-SPP)

#correcting for LA table
LA_transf <- read_xlsx(paste0(raw_data_dir,"4_LA_Corrections_DK.xlsx"))%>%
  select(-description, -motivation)%>%
  na.omit()

#Leaf Mass Area ratio data
ratioDKF <- read.csv(paste0(raw_data_dir,"4_LMA_data_DKF.csv"))%>%
  select(-Count,-Average.Size,-X.Area,-Total.weight,-FreshLA.DryMass)%>%
  mutate(`FreshLA/DryMass`=Total.Area/Sample.weight,
         Sample=substr(Slice, 1,9),
         Slice=substr(Slice, 1, 6))%>%
  rename(GENSPP=Slice)
ratioDKR <- read.csv(paste0(raw_data_dir,"4_LMA_data_DKR.csv"))%>%
  select(-Count,-Average.Size,-X.Area,-Total.weight,-FreshLA.DryMass)%>%
  mutate(`FreshLA/DryMass`=Total.Area/Sample.weight,
         Sample=substr(Slice, 1,9),
         Slice=substr(Slice, 1, 6))%>%
  rename(GENSPP=Slice)
ratio <- bind_rows(ratioDKR, ratioDKF)#combine sites
# Move Sample column before Total.Area
id_pos2 <- which(names(ratio) == "Total.Area")
ratio <- ratio[, c(1:(id_pos2 - 1), ncol(ratio), id_pos2:(ncol(ratio) - 1))]

#Read in sample LEAF dry_wts
leaf_dw <- read_xlsx(paste0(raw_data_dir,"4_LEAF_dryweights_DK.xlsx"))%>%
  select(-`Total weight`,-...6,-Notes,-SPP,-Sample)%>%
  rename(Leaf.dw=`Sample weight`)

}
 
{
# confirm datetime+add RWC----
mass_data2 <- left_join(mass_data,dry_wt, by="ID_Code") %>%
  mutate(Datetime = ymd_hms(Date_Time),Sample_no=as.factor(Sample_no))%>%
  select(-Date_Time)%>%
  group_by(ID_Code)%>%
  mutate(RWC = (Weight - dry_wt)/(max(Weight)-dry_wt))

mass_data4 <- mass_data2 %>%
  group_by(ID_Code)%>%
  mutate(seq = row_number())

mass_data5 <- mass_data4%>%
  mutate(Form = substr(unique(ID_Code), start = 1, stop = 3))

#clean up datasets
rm(list=c("mass_data","mass_data2","mass_data4"))
}






#join with LA corrections table----
mass_data5 <- left_join(mass_data5,LA_transf, by="GENSPP")

#change variable names to match code below
mass_data5 <- mass_data5%>%
  rename(Temp = Treat)%>%
  mutate(RWC = RWC*100)

#use mean LMA ratios----
ratio.summ <- ratio%>%
  group_by(GENSPP)%>%
  summarise(`FreshLA/DryMass`=mean(`FreshLA/DryMass`))

ratio.summ_comb <- ratio.summ

#adjoin LMA ratios & Leaf dw----
mass_data5 <- left_join(mass_data5,leaf_dw, by="ID_Code")
mass_data5 <- left_join(mass_data5,ratio.summ_comb, by="GENSPP")
#calc & correct LA
mass_data5 <- mass_data5%>%
  mutate(Calcd.Proj_Area=`FreshLA/DryMass`*Leaf.dw,
         Corr.Proj_Area=Calcd.Proj_Area*trans_onesided)


#save files for use in other scripts
saveRDS(mass_data5,paste0(data_dir,"mass_data_run4.rds"))
  