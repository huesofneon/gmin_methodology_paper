#script for finalising Table data & calculating+plotting Stomatal Margin retention Index (SMRI) for run3+run4ELRYHI

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
df_comb <- readRDS(paste0(data_dir,"df_comb_run3+run4ELYRHI.rds"))
gmin_modelled_summ <- readRDS(paste0(data_dir,"gmin_modelled_summ_run3+run4ELYRHI.rds"))
gmin_modelled_full <- readRDS(paste0(data_dir,"gmin_modelled_full_run3+run4ELYRHI.rds"))

#read in WP, etc. data for run3+run4ELYRHI####
RWC_limits <- read_xlsx(paste0(raw_data_dir,"Arens_RWC_limits_JK+run4ELYRHI.xlsx"))%>%
  select(-`-1/p50`, -`-1/p50.se`, -`-1/p88`, -`-1/p88.se`, -`100-RWC_p50`, -`100-RWC_p88`)
RWC_limits <- head(RWC_limits,4)%>%
  rename(SPP = GENSPP)%>%
  mutate(SPP = recode(SPP,
                      "Cannomois congesta" = "CANCON",
                      "Erica monsoniana" = "ERIMON",
                      "Protea repens" = "PROREP",
                      "Dicerothamnus rhinocerotis" = "ELYRHI"))%>%
  rename(
    SSM_p50 = HSM_p50,
    SSM_p50.se = HSM_p50.se,
    SSM_p88 = HSM_p88,
    SSM_p88.se = HSM_p88.se
  ) 

#formula for standard error
se <- function(x) {
  return(sd(x) / sqrt(length(x)))
}

# Function to replace old letters with new letters and sort them alphabetically
replace_letters <- function(old_letters, mapping) {
  new_letters <- mapping[unlist(strsplit(old_letters, ""))]
  sorted_letters <- sort(new_letters)
  paste0(sorted_letters, collapse = "")
}

}


#Calculating SSMs & SMRIs----
gmin_SSM_SMRI_data_summ <- gmin_modelled_summ%>%
  left_join(RWC_limits, by = "SPP")%>%
  mutate(SMRIΨ50 = SSM_p50/p50_gmin,#mean SMRI as based on means
         SMRIΨ50.se = abs(SMRIΨ50)*sqrt((SSM_p50.se/SSM_p50)^2 + (p50_gmin.se/p50_gmin)^2),
         SMRIΨ88 = SSM_p88/p88_gmin,#mean SMRI as based on means
         SMRIΨ88.se = abs(SMRIΨ88)*sqrt((SSM_p88.se/SSM_p88)^2 + (p88_gmin.se/p88_gmin)^2)
  )

#save data for Tables
saveRDS(RWC_limits, paste0(data_dir,"Table1_physiological_thresholds_run3+run4ELYRHI.rds"))
saveRDS(gmin_SSM_SMRI_data_summ, paste0(data_dir,"Table2_gmin+SMRI88_data_summ_run3+run4ELYRHI.rds"))

#Stats----
##error bars----
{
  p88_gmin.se <- gmin_modelled_full%>%
    filter(threshold_interval == "P88")%>%#filter for SSM_p88 gmins
    mutate(Temp = Treat)%>%
    rename(Treatment = Treat)%>%
    mutate(Treatment = paste0(SPP," - ",Temp))%>%
    group_by(Treatment) %>%
    summarise(p88_gmin.se = se(gmin_CorrPLA))
  
  p88_modelled_full <- gmin_modelled_full%>%
    filter(threshold_interval == "P88")%>%#filter for SSM_p88 gmins
    distinct(ID_Code, .keep_all = TRUE)%>%
    select(ID_Code,SPP,Treat,gmin_CorrPLA)%>%
    rename(Temp = Treat)%>%
    mutate(Treatment = paste0(SPP," - ",Temp)
           ,p88_gmin = gmin_CorrPLA)%>%
    select(-gmin_CorrPLA)%>%#reorder
    left_join(p88_gmin.se, by = "Treatment")%>%
    left_join(RWC_limits, by = "SPP")%>%
    mutate(SMRIΨ88 = SSM_p88/p88_gmin)
}


#ANOVA,TukeyHSD,Letter generation----
###20C
#Supplementary Figure? p88_gmin----
{
  #1-way ANOVA (gmin)
  p88_gmin_20C_aov <- aov(p88_gmin~SPP,data =p88_modelled_full)
  summary(p88_gmin_20C_aov)
  
  #post-hoc to use
  p88_gmin_20C_thsd <- TukeyHSD(p88_gmin_20C_aov)
  p88_gmin_20C_thsd
  
  # Generate the letters for Tukey HSD results
  tukey_letters <- multcompLetters4(p88_gmin_20C_aov, p88_gmin_20C_thsd)
  
  # Extract the letters into a data frame
  letters_df <- as.data.frame.list(tukey_letters$SPP)
  
  # Add SPP column to match with the main data
  letters_df$SPP <- rownames(letters_df)
  rownames(letters_df) <- NULL
  
  # Reorder the data frame based on the SPP column
  letters_df <- letters_df[order(letters_df$SPP), ]
  
  # Create a mapping of old letters to new letters
  unique_letters <- unique(unlist(strsplit(letters_df$Letters, "")))
  new_letters <- letters[1:length(unique_letters)]
  names(new_letters) <- unique_letters
  
  # Reassign the letters based on the mapping
  letters_df$Letters <- sapply(letters_df$Letters, replace_letters, mapping = new_letters)
  
  # Merge with plot data
  plot_data <- gmin_SSM_SMRI_data_summ %>%
    left_join(letters_df, by = "SPP")
  
  #gmin_p88
  plot_gmin_p88_20C <- 
    ggplot(plot_data, aes(x = SPP, y= p88_gmin, fill = SPP)) +
    geom_col() +
    geom_errorbar(aes(ymin = p88_gmin - p88_gmin.se, ymax = p88_gmin + p88_gmin.se), 
                  width = 0.2) +
    labs(title = "gmin_SSMp88", x = "Species",
         y = bquote(g[min]~"(mmol.m"^{-2}*".s"^{-1}*")")) +
    scale_fill_manual(name = "Species",
                      values = c("yellowgreen","#8DA36D",
                                 "wheat1","deeppink2"),
                      labels = c("Cannomois congesta","Dicerothamnus rhinocerotis", 
                                 "Erica monsoniana","Protea repens")
    ) +
    coord_cartesian(ylim = c(0, 25)) +
    geom_text(aes(x = SPP, y = p88_gmin + p88_gmin.se + 0.5, label = Letters), vjust = 0) +
    theme_bw() +
    theme(
      # Remove gridlines
      panel.grid = element_blank(),
      # Customize facet heading text
      strip.text = element_text(face = "bold"),
      # Customize facet heading background
      strip.background = element_rect(fill = NA, color = NA)
    )
  
  plot_gmin_p88_20C
  # ggsave("10.1.3.2-gmin_p88 by Temp_run3.jpeg", plot = plot_gmin_p88_20C, device = "jpeg", path = paste0(fig_dir,"PLA_Arens_data/20C only/"), width = 15, height = 10, units = "cm", dpi = 300)
}

#SMRIΨ88----
{
  #1-way ANOVA (gmin)
  SMRIΨ88_20C_aov <- aov(SMRIΨ88~SPP,data = p88_modelled_full)
  summary(SMRIΨ88_20C_aov)
  
  #post-hoc to use
  SMRIΨ88_20C_thsd <- TukeyHSD(SMRIΨ88_20C_aov)
  SMRIΨ88_20C_thsd
  
  # Generate the letters for Tukey HSD results
  tukey_letters <- multcompLetters4(SMRIΨ88_20C_aov, SMRIΨ88_20C_thsd)
  
  # Extract the letters into a data frame
  letters_df <- as.data.frame.list(tukey_letters$SPP)
  
  # Add SPP column to match with the main data
  letters_df$SPP <- rownames(letters_df)
  rownames(letters_df) <- NULL
  
  # Reorder the data frame based on the SPP column
  letters_df <- letters_df[order(letters_df$SPP), ]
  
  # Create a mapping of old letters to new letters
  unique_letters <- unique(unlist(strsplit(letters_df$Letters, "")))
  new_letters <- letters[1:length(unique_letters)]
  names(new_letters) <- unique_letters
  
  # Reassign the letters based on the mapping
  letters_df$Letters <- sapply(letters_df$Letters, replace_letters, mapping = new_letters)
  
  # Merge with plot data
  plot_data <- gmin_SSM_SMRI_data_summ %>%
    left_join(letters_df, by = "SPP")
  
  #SMRIΨ88
  plot_SMRIΨ88_20C <- 
    ggplot(plot_data, aes(x = SPP, y= SMRIΨ88, fill = SPP)) +
    geom_col()+
    geom_errorbar(aes(ymin = SMRIΨ88 - SMRIΨ88.se, ymax = SMRIΨ88 + SMRIΨ88.se), 
                  width = 0.2) +
    labs(title = "SMRIΨ88 (20C)",x = "Species",
         y = bquote(SMRI[Ψ88]~"(MPa s m"^{2}*" mmol"^{-1}*")"))+
    scale_fill_manual(name = "Species",
                      values = c("yellowgreen","#8DA36D",
                                 "wheat1","deeppink2"),
                      labels = c("Cannomois congesta","Dicerothamnus rhinocerotis", 
                                 "Erica monsoniana","Protea repens")
    )+
    coord_cartesian(ylim = c(0, 0.8)) +
    geom_text(aes(x = SPP, y = SMRIΨ88 + SMRIΨ88.se + 0.025, label = Letters), vjust = 0)+
    theme_bw() +
    theme(
      # Remove gridlines
      panel.grid = element_blank(),
      # Customize facet heading text
      strip.text = element_text(face = "bold"),
      # Customize facet heading background
      strip.background = element_rect(fill = NA, color = NA)
    )
  plot_SMRIΨ88_20C
  #ggsave("SMRIΨ88 (20C)_run3.jpeg", device = "jpeg", path = CHANGEfig_dir, width = 20, height = 20, units = "cm", dpi = 300)
}

#Fig5: SMRIΨ88-Horizontal----
plot_data2 <- plot_data %>%
  mutate(SPP = factor(SPP, levels = c("PROREP","ERIMON","ELYRHI","CANCON")))

{
  plot_SMRIΨ88_20C_horiz <- 
    ggplot(plot_data2, aes(x = SPP, y= SMRIΨ88, fill = SPP)) +
    geom_col()+
    coord_flip()+
    geom_errorbar(aes(ymin = SMRIΨ88 - SMRIΨ88.se, ymax = SMRIΨ88 + SMRIΨ88.se), 
                  width = 0.2) +
    labs(x = "Species",
         y = bquote(SMRI[Ψ88]~"(MPa s m"^{2}*" mmol"^{-1}*")"))+
    scale_fill_manual(name = "Species",
                      values = c("deeppink2","wheat1",
                                 "#8DA36D","yellowgreen"),
                      labels = c("Protea repens","Erica monsoniana",
                                 "Dicerothamnus rhinocerotis","Cannomois congesta")
    )+
    geom_text(aes(x = SPP, y = SMRIΨ88 + SMRIΨ88.se + 0.025, label = Letters), vjust = 0)+
    theme_bw() +
    theme(
      # Remove gridlines
      panel.grid = element_blank(),
      # Customize facet heading text
      strip.text = element_text(face = "bold"),
      legend.text = element_text(face = "italic"),
      # Customize facet heading background
      strip.background = element_rect(fill = NA, color = NA),
      legend.position = "bottom",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank()
    )
  
  plot_SMRIΨ88_20C_horiz
  }

ggsave("Fig5_SMRI.png", 
       plot = plot_SMRIΨ88_20C_horiz, device = "png", 
       path = fig_dir, 
       width = 18, height = 10, units = "cm", dpi = 300)
