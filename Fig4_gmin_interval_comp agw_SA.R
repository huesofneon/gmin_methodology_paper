#script for plotting gmins over differing intervals for run3+run4ELRYHI

### modified by Adam West 19 March 2026, to remove several metrics and replot using projected leaf area



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

gmin_modelled_summ <- readRDS(paste0(data_dir,"gmin_modelled_summ_run3+run4ELYRHI_agw.rds"))
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
    name = "\nThreshold used in\ngmin calculation (Projected)",
    values = c(
      "gmin_SSMΨ88"      = "#66a61e",
      "gmin_SSMΨ50"      = "#7570b3",
      #"gmin_postTLP"   = "#1b9e77",
      #"gmin_no_threshold" = "#e7298a",
      #"gmin_preΨ88"   = "#e6ab02",
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
  
#  ##gmin_tlp line----
# geom_segment(
#   data = bigplot_data,
#   aes(
#     x = 0, 
#     xend = as.numeric( RWC_tlp), 
#     y =  tlp_gmin, 
#     yend =  tlp_gmin,
#     color = "gmin_postTLP"), 
#   linewidth = 1.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
# )+
  
#   ##gmin_none line----
# geom_segment(
#   data = bigplot_data,
#   aes(
#     x = 0, 
#     xend = 100, 
#     y =  none_gmin, 
#     yend =  none_gmin,
#     color = "gmin_no_threshold"), 
#   linewidth = 1.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
# )+
#   
#   ##gmin_p88only line----
# geom_segment(
#   data = bigplot_data,
#   aes(
#     x = as.numeric( RWC_p88), 
#     xend = 100, 
#     y =  p88only_gmin, 
#     yend =  p88only_gmin,
#     color = "gmin_preΨ88"), 
#   linewidth = 1.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
# )+
#   
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
  # #tlp error bar
  # geom_segment(
  #   data = bigplot_data,
  #   aes(
  #     x = (as.numeric( RWC_tlp) + 0) / 2, 
  #     xend = (as.numeric( RWC_tlp) + 0) / 2, 
  #     y = tlp_gmin - tlp_gmin.se, 
  #     yend = tlp_gmin + tlp_gmin.se
  #   ),
  #   color = "#1b9e77", linewidth = 0.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
  # )+
  # # bottom cap
  # geom_segment(
  #   data = bigplot_data,
  #   aes(
  #     x = (as.numeric( RWC_tlp) + 0) / 2 - 1, 
  #     xend = (as.numeric( RWC_tlp) + 0) / 2 + 1,   # adjust 2 to control cap width
  #     y = tlp_gmin - tlp_gmin.se, 
  #     yend = tlp_gmin - tlp_gmin.se
  #   ),
  #   color = "#1b9e77", linewidth = 0.5, alpha = 0.5
  # ) +
  # # top cap
  # geom_segment(
  #   data = bigplot_data,
  #   aes(
  #     x = (as.numeric( RWC_tlp) + 0) / 2 - 1, 
  #     xend = (as.numeric( RWC_tlp) + 0) / 2 + 1,
  #     y = tlp_gmin + tlp_gmin.se, 
  #     yend = tlp_gmin + tlp_gmin.se
  #   ),
  #   color = "#1b9e77", linewidth = 0.5, alpha = 0.5
  # )+
  # #none error bar
  # geom_segment(
  #   data = bigplot_data,
  #   aes(
  #     x = (100 + 0) / 2, 
  #     xend = (100 + 0) / 2, 
  #     y = none_gmin - none_gmin.se, 
  #     yend = none_gmin + none_gmin.se
  #   ),
  #   color = "#e7298a", linewidth = 0.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
  # )+
  # # bottom cap
  # geom_segment(
  #   data = bigplot_data,
  #   aes(
  #     x = (100 + 0) / 2 - 1, 
  #     xend = (100 + 0) / 2 + 1,   # adjust 2 to control cap width
  #     y = none_gmin - none_gmin.se, 
  #     yend = none_gmin - none_gmin.se
  #   ),
  #   color = "#e7298a", linewidth = 0.5, alpha = 0.5
  # ) +
  # # top cap
  # geom_segment(
  #   data = bigplot_data,
  #   aes(
  #     x = (100 + 0) / 2 - 1, 
  #     xend = (100 + 0) / 2 + 1,
  #     y = none_gmin + none_gmin.se, 
  #     yend = none_gmin + none_gmin.se
  #   ),
  #   color = "#e7298a", linewidth = 0.5, alpha = 0.5
  # )+
  # #p88 only error bar
  # geom_segment(
  #   data = bigplot_data,
  #   aes(
  #     x = (100 + as.numeric( RWC_p88)) / 2, 
  #     xend = (100 + as.numeric( RWC_p88)) / 2, 
  #     y = p88only_gmin - p88only_gmin.se, 
  #     yend = p88only_gmin + p88only_gmin.se
  #   ),
  #   color = "#e6ab02", linewidth = 0.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
  # )+
  # # bottom cap
  # geom_segment(
  #   data = bigplot_data,
  #   aes(
  #     x = (100 + as.numeric( RWC_p88)) / 2 - 1, 
  #     xend = (100 + as.numeric( RWC_p88)) / 2 + 1,   # adjust 2 to control cap width
  #     y = p88only_gmin - p88only_gmin.se, 
  #     yend = p88only_gmin - p88only_gmin.se
  #   ),
  #   color = "#e6ab02", linewidth = 0.5, alpha = 0.5
  # ) +
  # # top cap
  # geom_segment(
  #   data = bigplot_data,
  #   aes(
  #     x = (100 + as.numeric( RWC_p88)) / 2 - 1, 
  #     xend = (100 + as.numeric( RWC_p88)) / 2 + 1,
  #     y = p88only_gmin + p88only_gmin.se, 
  #     yend = p88only_gmin + p88only_gmin.se
  #   ),
  #   color = "#e6ab02", linewidth = 0.5, alpha = 0.5
  # )+
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
ggsave("Fig4_gmin_interval_comp_agw.png", device = "png", 
       path = fig_dir, 
       width = 22, height = 25, units = "cm", dpi = 300)





###### Surface area plot

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
    name = "\nThreshold used in\ngmin calculation (Surface)",
    values = c(
      "gmin_SSMΨ88"      = "#66a61e",
      "gmin_SSMΨ50"      = "#7570b3",
      #"gmin_postTLP"   = "#1b9e77",
      #"gmin_no_threshold" = "#e7298a",
      #"gmin_preΨ88"   = "#e6ab02",
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
    y =  p88_gmin_SA, 
    yend =  p88_gmin_SA,
    color = "gmin_SSMΨ88"), 
  linewidth = 1.5, alpha = 1, linetype = "solid", inherit.aes = FALSE
)+
  
  ##gmin_SSMΨ50 line----
geom_segment(
  data = bigplot_data,
  aes(
    x = as.numeric( RWC_tlp), 
    xend = as.numeric( RWC_p50), 
    y =  p50_gmin_SA, 
    yend =  p50_gmin_SA,
    color = "gmin_SSMΨ50"), 
  linewidth = 1.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
)+
  
   
##gmin_RWC80-50 line----
geom_segment(
  data = bigplot_data,
  aes(
    x = 50, 
    xend = 80, 
    y =  RWC80_50_gmin_SA, 
    yend =  RWC80_50_gmin_SA,
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
    y = p88_gmin_SA - p88_gmin_SA.se, 
    yend = p88_gmin_SA + p88_gmin_SA.se
  ),
  color = "#66a61e", linewidth = 0.5, alpha = 1, linetype = "solid", inherit.aes = FALSE
)+
  # bottom cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (as.numeric( RWC_tlp) + as.numeric( RWC_p88)) / 2 - 1, 
      xend = (as.numeric( RWC_tlp) + as.numeric( RWC_p88)) / 2 + 1,   # adjust 2 to control cap width
      y = p88_gmin_SA - p88_gmin_SA.se, 
      yend = p88_gmin_SA - p88_gmin_SA.se
    ),
    color = "#66a61e", linewidth = 0.5
  ) +
  # top cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (as.numeric( RWC_tlp) + as.numeric( RWC_p88)) / 2 - 1, 
      xend = (as.numeric( RWC_tlp) + as.numeric( RWC_p88)) / 2 + 1,
      y = p88_gmin_SA + p88_gmin_SA.se, 
      yend = p88_gmin_SA + p88_gmin_SA.se
    ),
    color = "#66a61e", linewidth = 0.5
  )+
  #p50 error bar
  geom_segment(
    data = bigplot_data,
    aes(
      x = (as.numeric( RWC_tlp) + as.numeric( RWC_p50)) / 2, 
      xend = (as.numeric( RWC_tlp) + as.numeric( RWC_p50)) / 2, 
      y = p50_gmin_SA - p50_gmin_SA.se, 
      yend = p50_gmin_SA + p50_gmin_SA.se
    ),
    color = "#7570b3", linewidth = 0.5, alpha = 0.5, linetype = "solid", inherit.aes = FALSE
  )+
  # bottom cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (as.numeric( RWC_tlp) + as.numeric( RWC_p50)) / 2 - 1, 
      xend = (as.numeric( RWC_tlp) + as.numeric( RWC_p50)) / 2 + 1,   # adjust 2 to control cap width
      y = p50_gmin_SA - p50_gmin_SA.se, 
      yend = p50_gmin_SA - p50_gmin_SA.se
    ),
    color = "#7570b3", linewidth = 0.5, alpha = 0.5
  ) +
  # top cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (as.numeric( RWC_tlp) + as.numeric( RWC_p50)) / 2 - 1, 
      xend = (as.numeric( RWC_tlp) + as.numeric( RWC_p50)) / 2 + 1,
      y = p50_gmin_SA + p50_gmin_SA.se, 
      yend = p50_gmin_SA + p50_gmin_SA.se
    ),
    color = "#7570b3", linewidth = 0.5, alpha = 0.5
  )+
 
  #RWC80-50 error bar
  geom_segment(
    data = bigplot_data,
    aes(
      x = (80 + 50) / 2, 
      xend = (80 + 50) / 2, 
      y = RWC80_50_gmin_SA - RWC80_50_gmin_SA.se, 
      yend = RWC80_50_gmin_SA + RWC80_50_gmin_SA.se
    ),
    color = "#117733", linewidth = 0.5, alpha = 1, linetype = "solid", inherit.aes = FALSE
  )+
  # bottom cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (80 + 50) / 2 - 1, 
      xend = (80 + 50) / 2 + 1,   # adjust 2 to control cap width
      y = RWC80_50_gmin_SA - RWC80_50_gmin_SA.se,
      yend = RWC80_50_gmin_SA - RWC80_50_gmin_SA.se
    ),
    color = "#117733", linewidth = 0.5
  ) +
  # top cap
  geom_segment(
    data = bigplot_data,
    aes(
      x = (80 + 50) / 2 - 1, 
      xend = (80 + 50) / 2 + 1,
      y = RWC80_50_gmin_SA + RWC80_50_gmin_SA.se,
      yend = RWC80_50_gmin_SA + RWC80_50_gmin_SA.se
    ),
    color = "#117733", linewidth = 0.5
  )


#save plot----
ggsave("Fig4_gmin_interval_comp_agw_SA.png", device = "png", 
       path = fig_dir, 
       width = 22, height = 25, units = "cm", dpi = 300)







##############
#### chat gpt combined figure

# script for plotting projected gmin with species-specific surface-area secondary y-axis
# Adam West, 19 March 2026

rm(list = ls())

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
library(patchwork)

# Set paths -------------------------------------------------------------------
mainDir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
data_dir <- paste0(mainDir, "/data/")
raw_data_dir <- paste0(data_dir, "raw/")
fig_dir <- paste0(mainDir, "/figures/")

# Read data -------------------------------------------------------------------
RWC_limits <- read_xlsx(paste0(raw_data_dir, "Arens_RWC_limits_JK+run4ELYRHI.xlsx")) %>%
  select(-`-1/p50`, -`-1/p50.se`, -`-1/p88`, -`-1/p88.se`, -`100-RWC_p50`, -`100-RWC_p88`) %>%
  head(4) %>%
  rename(SPP = GENSPP) %>%
  mutate(Species = recode(
    SPP,
    "CANCON" = "Cannomois congesta",
    "ELYRHI" = "Dicerothamnus rhinocerotis",
    "ERIMON" = "Erica monsoniana",
    "PROREP" = "Protea repens"
  ))

gmin_modelled_summ <- readRDS(paste0(data_dir, "gmin_modelled_summ_run3+run4ELYRHI_agw.rds"))

thresholds <- RWC_limits %>%
  select(Species, SPP, RWC_tlp, RWC_p50, RWC_p88)

bigplot_data <- gmin_modelled_summ %>%
  left_join(thresholds, by = "SPP")

# Species labels --------------------------------------------------------------
species_labs <- c(
  "ELYRHI" = "Dicerothamnus rhinocerotis",
  "PROREP" = "Protea repens",
  "CANCON" = "Cannomois congesta",
  "ERIMON" = "Erica monsoniana"
)

# Function to get species-specific linear transform ---------------------------
# Maps surface-area values onto projected-value plotting space
get_sa_transform <- function(df_sp) {
  
  proj_vals <- c(
    df_sp$p88_gmin - df_sp$p88_gmin.se,
    df_sp$p88_gmin + df_sp$p88_gmin.se,
    df_sp$p50_gmin - df_sp$p50_gmin.se,
    df_sp$p50_gmin + df_sp$p50_gmin.se,
    df_sp$RWC80_50_gmin - df_sp$RWC80_50_gmin.se,
    df_sp$RWC80_50_gmin + df_sp$RWC80_50_gmin.se
  )
  
  sa_vals <- c(
    df_sp$p88_gmin_SA - df_sp$p88_gmin_SA.se,
    df_sp$p88_gmin_SA + df_sp$p88_gmin_SA.se,
    df_sp$p50_gmin_SA - df_sp$p50_gmin_SA.se,
    df_sp$p50_gmin_SA + df_sp$p50_gmin_SA.se,
    df_sp$RWC80_50_gmin_SA - df_sp$RWC80_50_gmin_SA.se,
    df_sp$RWC80_50_gmin_SA + df_sp$RWC80_50_gmin_SA.se
  )
  
  proj_rng <- range(proj_vals, na.rm = TRUE)
  sa_rng   <- range(sa_vals, na.rm = TRUE)
  
  proj_pad <- diff(proj_rng) * 0.08
  sa_pad   <- diff(sa_rng) * 0.08
  
  if (proj_pad == 0) proj_pad <- max(abs(proj_rng[1]) * 0.08, 0.1)
  if (sa_pad == 0)   sa_pad   <- max(abs(sa_rng[1]) * 0.08, 0.1)
  
  proj_rng <- proj_rng + c(-proj_pad, proj_pad)
  sa_rng   <- sa_rng + c(-sa_pad, sa_pad)
  
  b <- diff(proj_rng) / diff(sa_rng)
  a <- proj_rng[1] - b * sa_rng[1]
  
  list(
    a = a,
    b = b,
    proj_rng = proj_rng,
    sa_rng = sa_rng
  )
}

# Plotting function for one species ------------------------------------------
make_species_plot <- function(df_sp, show_left_title = TRUE, show_right_title = TRUE) {
  
  tr <- get_sa_transform(df_sp)
  a <- tr$a
  b <- tr$b
  
  sa_to_proj <- function(y) a + b * y
  proj_to_sa <- function(y) (y - a) / b
  
  df_plot <- df_sp %>%
    mutate(
      p88_gmin_SA_plot       = sa_to_proj(p88_gmin_SA),
      p50_gmin_SA_plot       = sa_to_proj(p50_gmin_SA),
      RWC80_50_gmin_SA_plot  = sa_to_proj(RWC80_50_gmin_SA),
      
      p88_gmin_SA_lo_plot      = sa_to_proj(p88_gmin_SA - p88_gmin_SA.se),
      p88_gmin_SA_hi_plot      = sa_to_proj(p88_gmin_SA + p88_gmin_SA.se),
      p50_gmin_SA_lo_plot      = sa_to_proj(p50_gmin_SA - p50_gmin_SA.se),
      p50_gmin_SA_hi_plot      = sa_to_proj(p50_gmin_SA + p50_gmin_SA.se),
      RWC80_50_gmin_SA_lo_plot = sa_to_proj(RWC80_50_gmin_SA - RWC80_50_gmin_SA.se),
      RWC80_50_gmin_SA_hi_plot = sa_to_proj(RWC80_50_gmin_SA + RWC80_50_gmin_SA.se)
    )
  
  ggplot(df_plot) +
    # background interval
    geom_rect(
      aes(xmin = 50, xmax = 80, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      alpha = 0.25,
      fill = "gray50"
    ) +
    
    # threshold lines
    geom_vline(
      aes(xintercept = as.numeric(RWC_tlp)),
      linetype = "dashed", color = "black", linewidth = 1
    ) +
    geom_vline(
      aes(xintercept = as.numeric(RWC_p88)),
      linetype = "dotted", color = "black", linewidth = 1
    ) +
    
    # projected lines
    geom_segment(
      aes(
        x = as.numeric(RWC_tlp), xend = as.numeric(RWC_p88),
        y = p88_gmin, yend = p88_gmin
      ),
      color = "#66a61e", linewidth = 1.5
    ) +
    geom_segment(
      aes(
        x = as.numeric(RWC_tlp), xend = as.numeric(RWC_p50),
        y = p50_gmin, yend = p50_gmin
      ),
      color = "#7570b3", linewidth = 1.5, alpha = 0.5
    ) +
    geom_segment(
      aes(
        x = 50, xend = 80,
        y = RWC80_50_gmin, yend = RWC80_50_gmin
      ),
      color = "#117733", linewidth = 1.5
    ) +
    
    # projected error bars: p88
    geom_segment(
      aes(
        x = (as.numeric(RWC_tlp) + as.numeric(RWC_p88)) / 2,
        xend = (as.numeric(RWC_tlp) + as.numeric(RWC_p88)) / 2,
        y = p88_gmin - p88_gmin.se,
        yend = p88_gmin + p88_gmin.se
      ),
      color = "#66a61e", linewidth = 0.5
    ) +
    geom_segment(
      aes(
        x = (as.numeric(RWC_tlp) + as.numeric(RWC_p88)) / 2 - 1,
        xend = (as.numeric(RWC_tlp) + as.numeric(RWC_p88)) / 2 + 1,
        y = p88_gmin - p88_gmin.se,
        yend = p88_gmin - p88_gmin.se
      ),
      color = "#66a61e", linewidth = 0.5
    ) +
    geom_segment(
      aes(
        x = (as.numeric(RWC_tlp) + as.numeric(RWC_p88)) / 2 - 1,
        xend = (as.numeric(RWC_tlp) + as.numeric(RWC_p88)) / 2 + 1,
        y = p88_gmin + p88_gmin.se,
        yend = p88_gmin + p88_gmin.se
      ),
      color = "#66a61e", linewidth = 0.5
    ) +
    
    # projected error bars: p50
    geom_segment(
      aes(
        x = (as.numeric(RWC_tlp) + as.numeric(RWC_p50)) / 2,
        xend = (as.numeric(RWC_tlp) + as.numeric(RWC_p50)) / 2,
        y = p50_gmin - p50_gmin.se,
        yend = p50_gmin + p50_gmin.se
      ),
      color = "#7570b3", linewidth = 0.5, alpha = 0.5
    ) +
    geom_segment(
      aes(
        x = (as.numeric(RWC_tlp) + as.numeric(RWC_p50)) / 2 - 1,
        xend = (as.numeric(RWC_tlp) + as.numeric(RWC_p50)) / 2 + 1,
        y = p50_gmin - p50_gmin.se,
        yend = p50_gmin - p50_gmin.se
      ),
      color = "#7570b3", linewidth = 0.5, alpha = 0.5
    ) +
    geom_segment(
      aes(
        x = (as.numeric(RWC_tlp) + as.numeric(RWC_p50)) / 2 - 1,
        xend = (as.numeric(RWC_tlp) + as.numeric(RWC_p50)) / 2 + 1,
        y = p50_gmin + p50_gmin.se,
        yend = p50_gmin + p50_gmin.se
      ),
      color = "#7570b3", linewidth = 0.5, alpha = 0.5
    ) +
    
    # projected error bars: RWC80-50
    geom_segment(
      aes(
        x = 65, xend = 65,
        y = RWC80_50_gmin - RWC80_50_gmin.se,
        yend = RWC80_50_gmin + RWC80_50_gmin.se
      ),
      color = "#117733", linewidth = 0.5
    ) +
    geom_segment(
      aes(
        x = 64, xend = 66,
        y = RWC80_50_gmin - RWC80_50_gmin.se,
        yend = RWC80_50_gmin - RWC80_50_gmin.se
      ),
      color = "#117733", linewidth = 0.5
    ) +
    geom_segment(
      aes(
        x = 64, xend = 66,
        y = RWC80_50_gmin + RWC80_50_gmin.se,
        yend = RWC80_50_gmin + RWC80_50_gmin.se
      ),
      color = "#117733", linewidth = 0.5
    ) +
    
    # surface lines, transformed onto projected axis
    geom_segment(
      aes(
        x = as.numeric(RWC_tlp), xend = as.numeric(RWC_p88),
        y = p88_gmin_SA_plot, yend = p88_gmin_SA_plot
      ),
      color = "#66a61e", linewidth = 1.1, linetype = "22"
    ) +
    geom_segment(
      aes(
        x = as.numeric(RWC_tlp), xend = as.numeric(RWC_p50),
        y = p50_gmin_SA_plot, yend = p50_gmin_SA_plot
      ),
      color = "#7570b3", linewidth = 1.1, linetype = "22", alpha = 0.5
    ) +
    geom_segment(
      aes(
        x = 50, xend = 80,
        y = RWC80_50_gmin_SA_plot, yend = RWC80_50_gmin_SA_plot
      ),
      color = "#117733", linewidth = 1.1, linetype = "22"
    ) +
    
    # surface error bars
    geom_segment(
      aes(
        x = (as.numeric(RWC_tlp) + as.numeric(RWC_p88)) / 2,
        xend = (as.numeric(RWC_tlp) + as.numeric(RWC_p88)) / 2,
        y = p88_gmin_SA_lo_plot,
        yend = p88_gmin_SA_hi_plot
      ),
      color = "#66a61e", linewidth = 0.4
    ) +
    geom_segment(
      aes(
        x = (as.numeric(RWC_tlp) + as.numeric(RWC_p50)) / 2,
        xend = (as.numeric(RWC_tlp) + as.numeric(RWC_p50)) / 2,
        y = p50_gmin_SA_lo_plot,
        yend = p50_gmin_SA_hi_plot
      ),
      color = "#7570b3", linewidth = 0.4, alpha = 0.5
    ) +
    geom_segment(
      aes(
        x = 65, xend = 65,
        y = RWC80_50_gmin_SA_lo_plot,
        yend = RWC80_50_gmin_SA_hi_plot
      ),
      color = "#117733", linewidth = 0.4
    ) +
    
    scale_x_reverse(
      limits = c(100, 0),
      breaks = seq(100, 0, -20),
      labels = seq(100, 0, -20)
    ) +
    scale_y_continuous(
      limits = tr$proj_rng,
      name = expression(g[min]~"(mmol m"^{-2}*" s"^{-1}*") projected"),
      sec.axis = sec_axis(
        trans = ~ proj_to_sa(.),
        name = expression(g[min]~"(mmol m"^{-2}*" s"^{-1}*") surface")
      )
    ) +
    labs(
      x = "RWC (%)",
      title = species_labs[df_sp$SPP[1]]
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "italic", hjust = 0.5, size = 11),
      legend.position = "none",
      axis.title.x = element_text(size = 10),
      axis.text.y.left = element_text(size = 9),
      axis.text.y.right = element_text(size = 9),
      axis.ticks.y.left = element_line(),
      axis.ticks.y.right = element_line(),
      axis.title.y.left = if (show_left_title) element_text(size = 10) else element_blank(),
      axis.title.y.right = if (show_right_title) element_text(size = 10) else element_blank()
    )
}

# Make four species panels ----------------------------------------------------
p1 <- make_species_plot(
  filter(bigplot_data, SPP == "ELYRHI"),
  show_left_title = TRUE, show_right_title = FALSE
)

p2 <- make_species_plot(
  filter(bigplot_data, SPP == "PROREP"),
  show_left_title = FALSE, show_right_title = TRUE
)

p3 <- make_species_plot(
  filter(bigplot_data, SPP == "CANCON"),
  show_left_title = TRUE, show_right_title = FALSE
)

p4 <- make_species_plot(
  filter(bigplot_data, SPP == "ERIMON"),
  show_left_title = FALSE, show_right_title = TRUE
)

combined_plot <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(10, 10, 10, 10)
    )
  )

# Display ---------------------------------------------------------------------
combined_plot

# Save ------------------------------------------------------------------------
ggsave(
  filename = "Fig4_gmin_interval_comp_projected_with_surface_secondary_axis.png",
  plot = combined_plot,
  path = fig_dir,
  width = 22,
  height = 25,
  units = "cm",
  dpi = 300
)