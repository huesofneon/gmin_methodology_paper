#script for calculating gmin over RWC80-50 interval for run4

### modified Adam West 19 March 


{
# Initialization ---------------------------------------------------------------
rm(list = ls()) # Clear environment

library(lubridate)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(readxl)
library(minpack.lm)
library(gridExtra)
library(multcompView)
library(cowplot)
library(grid)
}

# Set paths, values & functions  -----------------------------------------------------------------
{
mainDir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
data_dir <- paste0(mainDir,"/data/")
raw_data_dir <- paste0(data_dir,"raw/")
fig_dir<-paste0(mainDir,"/figures/")

#read in data-----
#Metdata: Using ave VPD over measuring period
metdata_summ <- readRDS(paste0(data_dir,"4_metdata_summary.rds"))%>%
  mutate(esMean = MeanVPD/(1-(MeanRH/100)))



#define atmospheric pressure at height above sea level 
atmosp <- 99891.70/1000 #kPa (@120m = approx height of HW Pearson)

#formula for standard error
se <- function(x) {
  return(sd(x) / sqrt(length(x)))
}

modelcoeffs_gmin <- readRDS(paste0(data_dir,"modelcoeffs_gmin_run4.rds"))
df_comb <- readRDS(paste0(data_dir,"df_comb_run4_agw.rds"))

treat_names_gmin <- sort(unique(df_comb$ID_Code))
spp_list <- unique(df_comb$GENSPP)
temp_list <- unique(df_comb$Temp)
samples <- unique(df_comb$Sample_no)

}


# gmin between RWC80 & RWC50 (gmin_RWC80-50)#####
{
  for (v in 1:length(treat_names_gmin)){
    if (v == 0){next}#skips problematic samples
    indiv_df <-df_comb%>%#subset df_comb to that one sample only
      filter(ID_Code == treat_names_gmin[v])%>%
      mutate(Calcd.Proj_Area = as.character(Calcd.Proj_Area))%>%
      mutate(Corr.Proj_Area = as.character(Corr.Proj_Area))
    
    #Set bin start & end
  RWC_bin_variables <- tibble(
    bin_start = 80,
    bin_end = 50)
  
  rownames(modelcoeffs_gmin) <- NULL #remove row names
  

  for (z in 1:nrow(RWC_bin_variables)){
    w <- modelcoeffs_gmin%>%
      filter(ID_Code == treat_names_gmin[v])#subset modelcoeffs_gmin to that one sample only
    
    w <- bind_rows(w,w) %>%
      mutate(fit = c(RWC_bin_variables$bin_start[z], RWC_bin_variables$bin_end[z]))%>%
      mutate(fit = as.numeric(fit))%>%
      mutate(Drydown_time = -1/k * log((fit - (100 - A))/A))%>%  # formula to calculate Drydown_time from fit
      mutate(Weight = fit/100 * ((max(indiv_df$Weight)) 
                                 - unique(indiv_df$dry_wt)) 
                                 + unique(indiv_df$dry_wt))%>%#calc weight from fit (modelled RWC)
      mutate(Calcd.Proj_Area=unique(indiv_df$Calcd.Proj_Area))%>%
      mutate(Corr.Proj_Area=unique(indiv_df$Corr.Proj_Area))
    
    #stack w rows
    w <- w %>%
      mutate(loop_no = z, 
             bin_start = RWC_bin_variables$bin_start[z],
             bin_end = RWC_bin_variables$bin_end[z],
             RWC_bin = paste0(RWC_bin_variables$bin_start[z],
                              "-",
                              RWC_bin_variables$bin_end[z]),
             waterloss = max(indiv_df$Weight)-w$Weight)
    
    #set NA exception
    for (h in 1:nrow(w)){
      if (is.nan(w$Drydown_time[h])){
        w$slope[h] <- NA
      }
      else {
        w$slope[h] <- coef(lm(w$waterloss~w$Drydown_time))[2]#else get lm coef
      }
    }
    
    #bind w2 rows
    if (z == 1) {
      w2 <- w
    }
    
    else {
      w2 <- rbind(w2, w)
      w2 <- w2 %>%
        distinct(RWC_bin, .keep_all = TRUE)
    }
  }
  
  #bind *_bin_values rows
  if (v == 1) {
    RWC80_50_bin_values <- w2
  }
  else {
    RWC80_50_bin_values <- rbind(RWC80_50_bin_values, w2)
  }
  }
  
  #remove NAs
  RWC80_50_bin_values <- na.omit(RWC80_50_bin_values)
  RWC80_50_bin_values <- RWC80_50_bin_values %>%
    rename("RWC (fit)" = fit)%>%
    rename(Treat = Temp)%>%
    left_join(metdata_summ, by = "Treat")#add VPD data for each temp
  
  RWC80_50_bin_values$Calcd.Proj_Area <- as.numeric(RWC80_50_bin_values$Calcd.Proj_Area)
  RWC80_50_bin_values$Corr.Proj_Area <- as.numeric(RWC80_50_bin_values$Corr.Proj_Area)
  
  RWC80_50_bin_values <- RWC80_50_bin_values %>% 
    mutate(`slope.s-1.metdata`=((slope/60)*(atmosp/MeanVPD)))%>%
    #gmin without correcting LA
    mutate(gmin_CalcdPLA=(((slope/60)/(18.015))*1000)/(Calcd.Proj_Area/10000)*(atmosp/MeanVPD))%>%
    #gmin with LA correction for surface structure
    mutate(gmin_CorrPLA=(((slope/60)/(18.015))*1000)/(Corr.Proj_Area/10000)*(atmosp/MeanVPD))

#Add Site
RWC80_50_bin_values <- RWC80_50_bin_values %>%
  mutate(Site=substr(RWC80_50_bin_values$ID_Code, 1, 3))

for (i in seq_along(RWC80_50_bin_values$Site)) {
  if (RWC80_50_bin_values$Site[i] == "DKR") {
    RWC80_50_bin_values$Site[i] <- "Renosterveld"
  }
  if (RWC80_50_bin_values$Site[i] == "DKF") {
    RWC80_50_bin_values$Site[i] <- "Fynbos"
  }
}

}  

#format & save----

#summ
  RWC80_50_summ <- RWC80_50_bin_values %>%
  mutate(Temp = Treat)%>%
  rename(Treatment = Treat)%>%
  mutate(Treatment = paste0(SPP," - ",Temp))%>%
  group_by(SPP, Temp, Treatment, Site) %>%
  summarise(
    CalcdPLA_RWC80_50_gmin = mean(gmin_CalcdPLA),
    CalcdPLA_RWC80_50_gmin.se = se(gmin_CalcdPLA),
    CorrPLA_RWC80_50_gmin = mean(gmin_CorrPLA),
    CorrPLA_RWC80_50_gmin.se = se(gmin_CorrPLA)
  )
  
saveRDS(RWC80_50_summ, paste0(data_dir,"gmin_modelled_summ_run4_agw.rds")) 

#full
saveRDS(RWC80_50_bin_values, paste0(data_dir,"gmin_modelled_full_run4_agw.rds"))
