#script for calculating gmin over various intervals for run3+run4ELYRHI

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
metdata_summ3 <- readRDS(paste0(data_dir,"3_metdata_summary.rds")) %>%
  mutate(esMean = MeanVPD / (1 - (MeanRH / 100))) %>%
  slice(rep(1, 3)) %>%                       # duplicate row 3 times
  mutate(SPP = c("PROREP", "CANCON", "ERIMON"))  # add new column
metdata_summ4 <- readRDS(paste0(data_dir,"4_metdata_summary.rds"))%>%
  mutate(esMean = MeanVPD/(1-(MeanRH/100))) %>%
  mutate(SPP = c("ELYRHI")) # add new column
#combine summs
metdata_summ <-bind_rows(metdata_summ3, metdata_summ4) %>%
  select(-Treat)
rm(list = "metdata_summ3","metdata_summ4") 

#define atmospheric pressure at height above sea level 
atmosp <- 99891.70/1000#kPa (@120m = approx height of HW Pearson)

#formula for standard error
se <- function(x) {
  return(sd(x) / sqrt(length(x)))
}

standardise_bins <- function(df) {
  df %>%
    mutate(
      bin_start = as.numeric(bin_start),
      bin_end   = as.numeric(bin_end)
    )
}

modelcoeffs_gmin <- readRDS(paste0(data_dir,"modelcoeffs_gmin_run3+run4ELYRHI.rds"))
df_comb <- readRDS(paste0(data_dir,"df_comb_run3+run4ELYRHI.rds"))

treat_names_gmin <- sort(unique(df_comb$ID_Code))
spp_list <- unique(df_comb$GENSPP)
temp_list <- unique(df_comb$Temp)
samples <- unique(df_comb$Sample_no)

}

# Part 1: gmin between RWC_tlp & RWC_p50 (gmin_SSM50)#####
{
  for (v in 1:length(treat_names_gmin)){
    if (v == 0){next}#skips problematic samples
    indiv_df <-df_comb%>%#subset df_comb to that one sample only
      filter(ID_Code == treat_names_gmin[v])%>%
      mutate(Calcd.Proj_Area = as.character(Calcd.Proj_Area))%>%
      mutate(Corr.Proj_Area = as.character(Corr.Proj_Area))
    
    #Set bin start & end
  RWC_bin_variables <- tibble(
    bin_start = max(indiv_df$RWC_tlp, na.rm = TRUE),
    bin_end = max(indiv_df$RWC_p50, na.rm = TRUE))
  
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
    p50_bin_values <- w2
  }
  else {
    p50_bin_values <- rbind(p50_bin_values, w2)
  }
  }
  
  #remove NAs
  p50_bin_values <- na.omit(p50_bin_values)
  p50_bin_values <- p50_bin_values %>%
    rename("RWC (fit)" = fit)%>%
    rename(Treat = Temp)%>%
    left_join(metdata_summ, by = "SPP")#add VPD data for each species
  
  p50_bin_values$Calcd.Proj_Area <- as.numeric(p50_bin_values$Calcd.Proj_Area)
  p50_bin_values$Corr.Proj_Area <- as.numeric(p50_bin_values$Corr.Proj_Area)
  
  p50_bin_values <- p50_bin_values %>% 
    mutate(`slope.s-1.metdata`=((slope/60)*(atmosp/MeanVPD)))%>%
    #gmin without correcting LA
    mutate(gmin_CalcdPLA=(((slope/60)/(18.015))*1000)/(Calcd.Proj_Area/10000)*(atmosp/MeanVPD))%>%
    #gmin with LA correction for 3D structure
    mutate(gmin_CorrPLA=(((slope/60)/(18.015))*1000)/(Corr.Proj_Area/10000)*(atmosp/MeanVPD))

#summarise p50_gmin
p50_summ <- p50_bin_values %>%
  mutate(Temp = Treat)%>%
  rename(Treatment = Treat)%>%
  mutate(Treatment = paste0(SPP," - ",Temp))%>%
  group_by(SPP, Temp, Treatment) %>%
  summarise(p50_gmin = mean(gmin_CorrPLA),
            p50_gmin.se = se(gmin_CorrPLA))
}

# Part 2: gmin between RWC_tlp & RWC_p88 (gmin_SSM88)#####
{
  for (v in 1:length(treat_names_gmin)){
    if (v == 0){next}#skips problematic samples
    indiv_df <-df_comb%>%#subset df_comb to that one sample only
      filter(ID_Code == treat_names_gmin[v])%>%
      mutate(Calcd.Proj_Area = as.character(Calcd.Proj_Area))%>%
      mutate(Corr.Proj_Area = as.character(Corr.Proj_Area))
    
    #Set bin start & end
  RWC_bin_variables <- tibble(
    bin_start = max(indiv_df$RWC_tlp, na.rm = TRUE),
    bin_end = max(indiv_df$RWC_p88, na.rm = TRUE))
  
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
    p88_bin_values <- w2
  }
  else {
    p88_bin_values <- rbind(p88_bin_values, w2)
  }
  }
  
  #remove NAs
  p88_bin_values <- na.omit(p88_bin_values)
  p88_bin_values <- p88_bin_values %>%
    rename("RWC (fit)" = fit)%>%
    rename(Treat = Temp)%>%
    left_join(metdata_summ, by = "SPP")#add VPD data for each species
  
  p88_bin_values$Calcd.Proj_Area <- as.numeric(p88_bin_values$Calcd.Proj_Area)
  p88_bin_values$Corr.Proj_Area <- as.numeric(p88_bin_values$Corr.Proj_Area)
  
  p88_bin_values <- p88_bin_values %>% 
    mutate(`slope.s-1.metdata`=((slope/60)*(atmosp/MeanVPD)))%>%
    #gmin without correcting LA
    mutate(gmin_CalcdPLA=(((slope/60)/(18.015))*1000)/(Calcd.Proj_Area/10000)*(atmosp/MeanVPD))%>%
    #gmin with LA correction for 3D structure
    mutate(gmin_CorrPLA=(((slope/60)/(18.015))*1000)/(Corr.Proj_Area/10000)*(atmosp/MeanVPD))

#summarise p88_gmin
  p88_summ <- p88_bin_values %>%
  mutate(Temp = Treat)%>%
  rename(Treatment = Treat)%>%
  mutate(Treatment = paste0(SPP," - ",Temp))%>%
  group_by(SPP, Temp, Treatment) %>%
  summarise(p88_gmin = mean(gmin_CorrPLA),
            p88_gmin.se = se(gmin_CorrPLA))
}  

# Part 3: gmin after RWC_tlp (gmin_postTLP)#####
{
  for (v in 1:length(treat_names_gmin)){
    if (v == 0){next}#skips problematic samples
    indiv_df <-df_comb%>%#subset df_comb to that one sample only
      filter(ID_Code == treat_names_gmin[v])%>%
      mutate(Calcd.Proj_Area = as.character(Calcd.Proj_Area))%>%
      mutate(Corr.Proj_Area = as.character(Corr.Proj_Area))
    
    #Set bin start & end
  RWC_bin_variables <- tibble(
    bin_start = max(indiv_df$RWC_tlp, na.rm = TRUE),
    bin_end = min(indiv_df$RWC, na.rm = TRUE))
  
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
    tlp_bin_values <- w2
  }
  else {
    tlp_bin_values <- rbind(tlp_bin_values, w2)
  }
  }
  
  #remove NAs
  tlp_bin_values <- na.omit(tlp_bin_values)
  tlp_bin_values <- tlp_bin_values %>%
    rename("RWC (fit)" = fit)%>%
    rename(Treat = Temp)%>%
    left_join(metdata_summ, by = "SPP")#add VPD data for each species
  
  tlp_bin_values$Calcd.Proj_Area <- as.numeric(tlp_bin_values$Calcd.Proj_Area)
  tlp_bin_values$Corr.Proj_Area <- as.numeric(tlp_bin_values$Corr.Proj_Area)
  
  tlp_bin_values <- tlp_bin_values %>% 
    mutate(`slope.s-1.metdata`=((slope/60)*(atmosp/MeanVPD)))%>%
    #gmin without correcting LA
    mutate(gmin_CalcdPLA=(((slope/60)/(18.015))*1000)/(Calcd.Proj_Area/10000)*(atmosp/MeanVPD))%>%
    #gmin with LA correction for 3D structure
    mutate(gmin_CorrPLA=(((slope/60)/(18.015))*1000)/(Corr.Proj_Area/10000)*(atmosp/MeanVPD))

#summarise tlp_gmin
  tlp_summ <- tlp_bin_values %>%
  mutate(Temp = Treat)%>%
  rename(Treatment = Treat)%>%
  mutate(Treatment = paste0(SPP," - ",Temp))%>%
  group_by(SPP, Temp, Treatment) %>%
  summarise(tlp_gmin = mean(gmin_CorrPLA),
            tlp_gmin.se = se(gmin_CorrPLA))
}  

# Part 4: gmin after RWC_100 (gmin_no_threshold)#####
{
  for (v in 1:length(treat_names_gmin)){
    if (v == 0){next}#skips problematic samples
    indiv_df <-df_comb%>%#subset df_comb to that one sample only
      filter(ID_Code == treat_names_gmin[v])%>%
      mutate(Calcd.Proj_Area = as.character(Calcd.Proj_Area))%>%
      mutate(Corr.Proj_Area = as.character(Corr.Proj_Area))
    
    #Set bin start & end
  RWC_bin_variables <- tibble(
    bin_start = 100,
    bin_end = min(indiv_df$RWC, na.rm = TRUE))
  
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
    none_bin_values <- w2
  }
  else {
    none_bin_values <- rbind(none_bin_values, w2)
  }
  }
  
  #remove NAs
  none_bin_values <- na.omit(none_bin_values)
  none_bin_values <- none_bin_values %>%
    rename("RWC (fit)" = fit)%>%
    rename(Treat = Temp)%>%
    left_join(metdata_summ, by = "SPP")#add VPD data for each species
  
  none_bin_values$Calcd.Proj_Area <- as.numeric(none_bin_values$Calcd.Proj_Area)
  none_bin_values$Corr.Proj_Area <- as.numeric(none_bin_values$Corr.Proj_Area)
  
  none_bin_values <- none_bin_values %>% 
    mutate(`slope.s-1.metdata`=((slope/60)*(atmosp/MeanVPD)))%>%
    #gmin without correcting LA
    mutate(gmin_CalcdPLA=(((slope/60)/(18.015))*1000)/(Calcd.Proj_Area/10000)*(atmosp/MeanVPD))%>%
    #gmin with LA correction for 3D structure
    mutate(gmin_CorrPLA=(((slope/60)/(18.015))*1000)/(Corr.Proj_Area/10000)*(atmosp/MeanVPD))

#summarise none_gmin
  none_summ <- none_bin_values %>%
  mutate(Temp = Treat)%>%
  rename(Treatment = Treat)%>%
  mutate(Treatment = paste0(SPP," - ",Temp))%>%
  group_by(SPP, Temp, Treatment) %>%
  summarise(none_gmin = mean(gmin_CorrPLA),
            none_gmin.se = se(gmin_CorrPLA))
}  

# Part 5: gmin before p88 (gmin_preΨ88)#####
{
  for (v in 1:length(treat_names_gmin)){
    if (v == 0){next}#skips problematic samples
    indiv_df <-df_comb%>%#subset df_comb to that one sample only
      filter(ID_Code == treat_names_gmin[v])%>%
      mutate(Calcd.Proj_Area = as.character(Calcd.Proj_Area))%>%
      mutate(Corr.Proj_Area = as.character(Corr.Proj_Area))
    
    #Set bin start & end
  RWC_bin_variables <- tibble(
    bin_start = 100,
    bin_end = max(indiv_df$RWC_p88, na.rm = TRUE))
  
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
    p88only_bin_values <- w2
  }
  else {
    p88only_bin_values <- rbind(p88only_bin_values, w2)
  }
  }
  
  #remove NAs
  p88only_bin_values <- na.omit(p88only_bin_values)
  p88only_bin_values <- p88only_bin_values %>%
    rename("RWC (fit)" = fit)%>%
    rename(Treat = Temp)%>%
    left_join(metdata_summ, by = "SPP")#add VPD data for each species
  
  p88only_bin_values$Calcd.Proj_Area <- as.numeric(p88only_bin_values$Calcd.Proj_Area)
  p88only_bin_values$Corr.Proj_Area <- as.numeric(p88only_bin_values$Corr.Proj_Area)
  
  p88only_bin_values <- p88only_bin_values %>% 
    mutate(`slope.s-1.metdata`=((slope/60)*(atmosp/MeanVPD)))%>%
    #gmin without correcting LA
    mutate(gmin_CalcdPLA=(((slope/60)/(18.015))*1000)/(Calcd.Proj_Area/10000)*(atmosp/MeanVPD))%>%
    #gmin with LA correction for 3D structure
    mutate(gmin_CorrPLA=(((slope/60)/(18.015))*1000)/(Corr.Proj_Area/10000)*(atmosp/MeanVPD))

#summarise p88only_gmin
  p88only_summ <- p88only_bin_values %>%
  mutate(Temp = Treat)%>%
  rename(Treatment = Treat)%>%
  mutate(Treatment = paste0(SPP," - ",Temp))%>%
  group_by(SPP, Temp, Treatment) %>%
  summarise(p88only_gmin = mean(gmin_CorrPLA),
            p88only_gmin.se = se(gmin_CorrPLA))
}

# Part 6: gmin between RWC80 & RWC50 (gmin_RWC80-50)#####
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
    left_join(metdata_summ, by = "SPP")#add VPD data for each species
  
  RWC80_50_bin_values$Calcd.Proj_Area <- as.numeric(RWC80_50_bin_values$Calcd.Proj_Area)
  RWC80_50_bin_values$Corr.Proj_Area <- as.numeric(RWC80_50_bin_values$Corr.Proj_Area)
  
  RWC80_50_bin_values <- RWC80_50_bin_values %>% 
    mutate(`slope.s-1.metdata`=((slope/60)*(atmosp/MeanVPD)))%>%
    #gmin without correcting LA
    mutate(gmin_CalcdPLA=(((slope/60)/(18.015))*1000)/(Calcd.Proj_Area/10000)*(atmosp/MeanVPD))%>%
    #gmin with LA correction for 3D structure
    mutate(gmin_CorrPLA=(((slope/60)/(18.015))*1000)/(Corr.Proj_Area/10000)*(atmosp/MeanVPD))

#summarise RWC80_50_gmin
  RWC80_50_summ <- RWC80_50_bin_values %>%
  mutate(Temp = Treat)%>%
  rename(Treatment = Treat)%>%
  mutate(Treatment = paste0(SPP," - ",Temp))%>%
  group_by(SPP, Temp, Treatment) %>%
  summarise(RWC80_50_gmin = mean(gmin_CorrPLA),
            RWC80_50_gmin.se = se(gmin_CorrPLA))
}  

#Combine summ dataframes----
gmin_modelled_summ <- reduce(
  list(
    p50_summ,
    p88_summ,
    tlp_summ,
    none_summ,
    p88only_summ,
    RWC80_50_summ
  ),
  left_join,
  by = c("SPP", "Temp", "Treatment")
)

saveRDS(gmin_modelled_summ, paste0(data_dir,"gmin_modelled_summ_run3+run4ELYRHI.rds"))

#Combine full dataframes----
gmin_modelled_full <- bind_rows(
  P50        = standardise_bins(p50_bin_values),
  P88        = standardise_bins(p88_bin_values),
  TLP        = standardise_bins(tlp_bin_values),
  None       = standardise_bins(none_bin_values),
  P88_only   = standardise_bins(p88only_bin_values),
  RWC80_50   = standardise_bins(RWC80_50_bin_values),
  .id = "threshold_interval"
)

saveRDS(gmin_modelled_full, paste0(data_dir,"gmin_modelled_full_run3+run4ELYRHI.rds"))
