#script for plotting gmin_RWC80-50 

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

#formula for standard error
se <- function(x) {
  return(sd(x) / sqrt(length(x)))
}

# Function to remove legend AND format x axis AND remove titles
remove_legend <- function(plot) {
  plot + theme(legend.position = "none",
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.x = element_blank(),
               plot.title = element_blank()
               )
}

# Function to replace old letters with new letters and sort them alphabetically
replace_letters <- function(old_letters, mapping) {
  new_letters <- mapping[unlist(strsplit(old_letters, ""))]
  sorted_letters <- sort(new_letters)
  paste0(sorted_letters, collapse = "")
}
}

#read in data for run 4-----
# {#Metdata: Using ave VPD over measuring period
# metdata_summ <- readRDS(paste0(data_dir,"4_metdata_summary.rds"))%>%
#   mutate(esMean = MeanVPD/(1-(MeanRH/100)))
#   
# modelcoeffs_gmin <- readRDS(paste0(data_dir,"PLA_modelcoeffs_gmin_run4.rds"))
# df_comb <- readRDS(paste0(data_dir,"PLA_df_comb_run4.rds"))
# 
# 
# treat_names_gmin <- sort(unique(df_comb$ID_Code))
# spp_list <- unique(df_comb$GENSPP)
# temp_list <- unique(df_comb$Temp)
# samples <- unique(df_comb$Sample_no)
# }

gmin_modelled_summ <- readRDS(paste0(data_dir,"gmin_modelled_summ_run4.rds"))
gmin_modelled_full <- readRDS(paste0(data_dir,"gmin_modelled_full_run4.rds"))
# df_comb <- readRDS(paste0(data_dir,"df_comb_run4.rds"))

**REMOVE `` VARS

#aov & Tukey HSD for Corr.PLA----
{
aov_data <- gmin_modelled_full %>%
  distinct(ID_Code, .keep_all = TRUE)%>% #delete duped rows
  filter(SPP != "ELYRHI") #remove ELYRHI
  

#1-way ANOVA (gmin)
gmin_aov <- aov(gmin_CorrPLA~SPP, data = aov_data)
summary(gmin_aov)

#post-hoc to use
gmin_thsd <- TukeyHSD(gmin_aov)
gmin_thsd

# Generate the letters for Tukey HSD results
tukey_letters <- multcompLetters4(gmin_aov, gmin_thsd)

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
plot_letters_df <- gmin_modelled_summ %>%
  left_join(letters_df, by = "SPP")
}

#gmin_RWC80_50_Corr.PLA----
ggplot(filter(plot_letters_df, SPP != "ELYRHI"), 
       aes(x = SPP, 
           y = CorrPLA_RWC80_50_gmin, 
           fill = Site)) +
  geom_col() +
  scale_fill_manual(
    name = "Biome",
    values = c(
      "Renosterveld" = "#1b9e77",
      "Fynbos" = "#d95f02"
    )
  )+
  scale_x_discrete(
    labels = c(
      "AGACAP" = "Agathosma\ncapensis",
      "ASPSHA" = "Aspalathus\nshawii",
      "CLURUB" = "Clutia\nrubricaulis",
      "ELYRHI" = "Dicerothamnus\nrhinocerotis",
      "LOBDEC" = "Lobostemon\ndecorus",
      "MICPOL" = "Microdon\npolygaloides",
      "OEDSQU" = "Oedera\nsquarrosa",
      "PASOBT" = "Passerina\nobtusifolia",
      "PROLAU" = "Protea\nlaurifolia",
      "PROLOR" = "Protea\nlorifolia",
      "RUSMUL" = "Ruschia\nmultiflora",
      "SELDOL" = "Selago\ndolosa",
      "WAHNOD" = "Wahlenbergia\nneorigida"
    )
  )+
  geom_errorbar(aes(ymin = CorrPLA_RWC80_50_gmin - CorrPLA_RWC80_50_gmin.se, 
                    ymax = CorrPLA_RWC80_50_gmin + CorrPLA_RWC80_50_gmin.se), 
                width = 0.2) +
  labs(x = element_blank(),
       y = bquote(g[min]~"(mmol m"^{-2}*" s"^{-1}*")"),
       title = " ") +
  coord_cartesian(ylim = c(0, 75)) +
  geom_text(aes(x = SPP, 
                y = CorrPLA_RWC80_50_gmin + CorrPLA_RWC80_50_gmin.se + 1.5, 
                 label = Letters
                ), 
            vjust = 0) +
  theme_bw() +
  theme(
    # Remove gridlines
    panel.grid = element_blank(),
    # Customize facet heading text
    strip.text = element_text(face = "bold"),
    # Customize facet heading background
    strip.background = element_rect(fill = NA, color = NA),
    # Make x-axis labels italic + diagonal
    axis.text.x = element_text(
      face = "italic",
      angle = 30,
      hjust = 1,  # right-align
      vjust = 1
    )
  )

#save
ggsave("Fig6_gmin_comp_run4.png", device = "png", 
       path = fig_dir, 
       width = 30, height = 15, units = "cm", dpi = 300)
