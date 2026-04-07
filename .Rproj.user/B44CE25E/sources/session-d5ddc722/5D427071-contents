#script for plotting gmins over differing intervals for run3+run4ELRYHI

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
RWC_limits <- read_xlsx(paste0(raw_data_dir,"Arens_RWC_limits_JK+run4ELYRHI.xlsx"))%>%
  select(-`-1/p50`, -`-1/p50.se`, -`-1/p88`, -`-1/p88.se`, -`100-RWC_p50`, -`100-RWC_p88`)
RWC_limits <- head(RWC_limits,4)%>%
  rename(SPP = GENSPP)%>%
  mutate(Species = recode(SPP,
                          "CANCON" = "Cannomois congesta",
                          "ELYRHI" = "Dicerothamnus rhinocerotis",
                          "ERIMON" = "Erica monsoniana",
                          "PROREP" = "Protea repens"))

gmin_modelled_summ <- readRDS(paste0(data_dir,"gmin_modelled_summ_run3+run4ELYRHI.rds"))
}

# Reshape data for plotting
thresholds <- RWC_limits%>%
  select(Species, SPP, RWC_tlp, RWC_p50, RWC_p88)
#combine in one df
bigplot_data <- gmin_modelled_summ%>%
  left_join(thresholds, by = "SPP")

#######Big interval gmin comparison plot----------------------
ggplot() +
  ##setup visuals----
  facet_wrap(~ SPP, 
             labeller = labeller(SPP = c("ELYRHI" = "Dicerothamnus rhinocerotis",
                                         "PROREP" = "Protea repens",
                                         "CANCON" = "Cannomois congesta",
                                         "ERIMON" = "Erica monsoniana")),
             scales = "free_y") +
  labs(
       x = "RWC (%)",
       y = bquote(~g[min]~"(mmol m"^{-2}*" s"^{-1}*")")) +
  scale_colour_manual(
    name = "\nThreshold used in\ngmin calculation",
    values = c(
      "gmin_SSMΨ88"      = "#66a61e",
      "gmin_SSMΨ50"      = "#7570b3",
      "gmin_postTLP"   = "#1b9e77",
      "gmin_no_threshold" = "#e7298a",
      "gmin_preΨ88"   = "#e6ab02",
      "gmin_RWC80-50"  = "#117733"
    )
  )+
  #reverse reverse!
  scale_x_reverse(
    limits = c(100, 0),
    breaks = seq(100, 0, -20),
    labels = seq(100, 0, -20)
  )+
  theme_bw() +
  theme(
    # Remove gridlines
    panel.grid = element_blank(),
    # Customize facet heading text
    strip.text = element_text(face = "italic"),
    # Customize facet heading background
    strip.background = element_rect(fill = NA, color = NA,
    ),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(margin = margin(b = 20, unit = "pt"),# Adjust 'b' for bottom margin
                                #legend.title = element_text(margin = margin(b = -20, unit = "pt"),# Adjust 'b' for bottom margin
                                hjust = 0) # Left-align the legend title
  )+

  ##80-50 shaded interval----
  geom_rect(
    data = bigplot_data,
    aes(xmin = 50, 
        xmax = 80,
        ymin = -Inf,
        ymax = Inf,
        fill = NULL),  # Fill specified below
    inherit.aes = FALSE,
    alpha = 0.25,
    fill = "gray50"
  )+  
  ##tlp vert line----
  geom_vline(data = bigplot_data, aes(xintercept = as.numeric(RWC_tlp)), 
             linetype = "dashed", color = "black",
             linewidth = 1, alpha = 1
  )+
  ##p88 vert line----
  geom_vline(data = bigplot_data, aes(xintercept = as.numeric(RWC_p88)), 
             linetype = "dotted", color = "black",
             linewidth = 1, alpha = 1
  )+
  
  
  ##gmin_SSMΨ88 line----
geom_segment(
  data = bigplot_data,
  aes(
    x = as.numeric( RWC_tlp), 
    xend = as.numeric( RWC_p88), 
    y =  p88_gmin, 
    yend =  p88_gmin,
    color = "gmin_SSMΨ88"), 
  linewidth = 1.5, alpha = 1, linetype = "solid", inherit.aes = FALSE
)+
  
  ##gmin_SSMΨ50 line----
geom_segment(
  data = bigplot_data,
  aes(
    x = as.numeric( RWC_tlp), 
    xend = as.numeric( RWC_p50), 
    y =  p50_gmin, 
    yend =  p50_gmin,
    color = "gmin_SSMΨ50"), 
  linewidth = 1.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
)+
  
 ##gmin_tlp line----
geom_segment(
  data = bigplot_data,
  aes(
    x = 0, 
    xend = as.numeric( RWC_tlp), 
    y =  tlp_gmin, 
    yend =  tlp_gmin,
    color = "gmin_postTLP"), 
  linewidth = 1.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
)+
  
  ##gmin_none line----
geom_segment(
  data = bigplot_data,
  aes(
    x = 0, 
    xend = 100, 
    y =  none_gmin, 
    yend =  none_gmin,
    color = "gmin_no_threshold"), 
  linewidth = 1.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
)+
  
  ##gmin_p88only line----
geom_segment(
  data = bigplot_data,
  aes(
    x = as.numeric( RWC_p88), 
    xend = 100, 
    y =  p88only_gmin, 
    yend =  p88only_gmin,
    color = "gmin_preΨ88"), 
  linewidth = 1.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
)+
  
  ##gmin_RWC80-50 line----
geom_segment(
  data = bigplot_data,
  aes(
    x = 50, 
    xend = 80, 
    y =  RWC80_50_gmin, 
    yend =  RWC80_50_gmin,
    color = "gmin_RWC80-50"), 
  linewidth = 1.5, alpha = 1, linetype = "solid", inherit.aes = FALSE
)+
  ##error bars----
#p88 error bar
geom_segment(
  data = bigplot_data,
  aes(
    x = (as.numeric( RWC_tlp) + as.numeric( RWC_p88)) / 2, 
    xend = (as.numeric( RWC_tlp) + as.numeric( RWC_p88)) / 2, 
    y = p88_gmin - p88_gmin.se, 
    yend = p88_gmin + p88_gmin.se
  ),
  color = "#66a61e", linewidth = 0.5, alpha = 1, linetype = "solid", inherit.aes = FALSE
)+
  # bottom cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (as.numeric( RWC_tlp) + as.numeric( RWC_p88)) / 2 - 1, 
      xend = (as.numeric( RWC_tlp) + as.numeric( RWC_p88)) / 2 + 1,   # adjust 2 to control cap width
      y = p88_gmin - p88_gmin.se, 
      yend = p88_gmin - p88_gmin.se
    ),
    color = "#66a61e", linewidth = 0.5
  ) +
  # top cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (as.numeric( RWC_tlp) + as.numeric( RWC_p88)) / 2 - 1, 
      xend = (as.numeric( RWC_tlp) + as.numeric( RWC_p88)) / 2 + 1,
      y = p88_gmin + p88_gmin.se, 
      yend = p88_gmin + p88_gmin.se
    ),
    color = "#66a61e", linewidth = 0.5
  )+
  #p50 error bar
  geom_segment(
    data = bigplot_data,
    aes(
      x = (as.numeric( RWC_tlp) + as.numeric( RWC_p50)) / 2, 
      xend = (as.numeric( RWC_tlp) + as.numeric( RWC_p50)) / 2, 
      y = p50_gmin - p50_gmin.se, 
      yend = p50_gmin + p50_gmin.se
    ),
    color = "#7570b3", linewidth = 0.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
  )+
  # bottom cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (as.numeric( RWC_tlp) + as.numeric( RWC_p50)) / 2 - 1, 
      xend = (as.numeric( RWC_tlp) + as.numeric( RWC_p50)) / 2 + 1,   # adjust 2 to control cap width
      y = p50_gmin - p50_gmin.se, 
      yend = p50_gmin - p50_gmin.se
    ),
    color = "#7570b3", linewidth = 0.5, alpha = 0.5
  ) +
  # top cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (as.numeric( RWC_tlp) + as.numeric( RWC_p50)) / 2 - 1, 
      xend = (as.numeric( RWC_tlp) + as.numeric( RWC_p50)) / 2 + 1,
      y = p50_gmin + p50_gmin.se, 
      yend = p50_gmin + p50_gmin.se
    ),
    color = "#7570b3", linewidth = 0.5, alpha = 0.5
  )+
  #tlp error bar
  geom_segment(
    data = bigplot_data,
    aes(
      x = (as.numeric( RWC_tlp) + 0) / 2, 
      xend = (as.numeric( RWC_tlp) + 0) / 2, 
      y = tlp_gmin - tlp_gmin.se, 
      yend = tlp_gmin + tlp_gmin.se
    ),
    color = "#1b9e77", linewidth = 0.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
  )+
  # bottom cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (as.numeric( RWC_tlp) + 0) / 2 - 1, 
      xend = (as.numeric( RWC_tlp) + 0) / 2 + 1,   # adjust 2 to control cap width
      y = tlp_gmin - tlp_gmin.se, 
      yend = tlp_gmin - tlp_gmin.se
    ),
    color = "#1b9e77", linewidth = 0.5, alpha = 0.5
  ) +
  # top cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (as.numeric( RWC_tlp) + 0) / 2 - 1, 
      xend = (as.numeric( RWC_tlp) + 0) / 2 + 1,
      y = tlp_gmin + tlp_gmin.se, 
      yend = tlp_gmin + tlp_gmin.se
    ),
    color = "#1b9e77", linewidth = 0.5, alpha = 0.5
  )+
  #none error bar
  geom_segment(
    data = bigplot_data,
    aes(
      x = (100 + 0) / 2, 
      xend = (100 + 0) / 2, 
      y = none_gmin - none_gmin.se, 
      yend = none_gmin + none_gmin.se
    ),
    color = "#e7298a", linewidth = 0.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
  )+
  # bottom cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (100 + 0) / 2 - 1, 
      xend = (100 + 0) / 2 + 1,   # adjust 2 to control cap width
      y = none_gmin - none_gmin.se, 
      yend = none_gmin - none_gmin.se
    ),
    color = "#e7298a", linewidth = 0.5, alpha = 0.5
  ) +
  # top cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (100 + 0) / 2 - 1, 
      xend = (100 + 0) / 2 + 1,
      y = none_gmin + none_gmin.se, 
      yend = none_gmin + none_gmin.se
    ),
    color = "#e7298a", linewidth = 0.5, alpha = 0.5
  )+
  #p88 only error bar
  geom_segment(
    data = bigplot_data,
    aes(
      x = (100 + as.numeric( RWC_p88)) / 2, 
      xend = (100 + as.numeric( RWC_p88)) / 2, 
      y = p88only_gmin - p88only_gmin.se, 
      yend = p88only_gmin + p88only_gmin.se
    ),
    color = "#e6ab02", linewidth = 0.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
  )+
  # bottom cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (100 + as.numeric( RWC_p88)) / 2 - 1, 
      xend = (100 + as.numeric( RWC_p88)) / 2 + 1,   # adjust 2 to control cap width
      y = p88only_gmin - p88only_gmin.se, 
      yend = p88only_gmin - p88only_gmin.se
    ),
    color = "#e6ab02", linewidth = 0.5, alpha = 0.5
  ) +
  # top cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (100 + as.numeric( RWC_p88)) / 2 - 1, 
      xend = (100 + as.numeric( RWC_p88)) / 2 + 1,
      y = p88only_gmin + p88only_gmin.se, 
      yend = p88only_gmin + p88only_gmin.se
    ),
    color = "#e6ab02", linewidth = 0.5, alpha = 0.5
  )+
  #RWC80-50 error bar
  geom_segment(
    data = bigplot_data,
    aes(
      x = (80 + 50) / 2, 
      xend = (80 + 50) / 2, 
      y = RWC80_50_gmin - RWC80_50_gmin.se, 
      yend = RWC80_50_gmin + RWC80_50_gmin.se
    ),
    color = "#117733", linewidth = 0.5, alpha = 1, linetype = "solid", inherit.aes = FALSE
  )+
  # bottom cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (80 + 50) / 2 - 1, 
      xend = (80 + 50) / 2 + 1,   # adjust 2 to control cap width
      y = RWC80_50_gmin - RWC80_50_gmin.se,
      yend = RWC80_50_gmin - RWC80_50_gmin.se
    ),
    color = "#117733", linewidth = 0.5
  ) +
  # top cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (80 + 50) / 2 - 1, 
      xend = (80 + 50) / 2 + 1,
      y = RWC80_50_gmin + RWC80_50_gmin.se,
      yend = RWC80_50_gmin + RWC80_50_gmin.se
    ),
    color = "#117733", linewidth = 0.5
  )

#save plot----
ggsave("Fig4_gmin_interval_comp.png", device = "png", 
       path = fig_dir, 
       width = 22, height = 25, units = "cm", dpi = 300)
