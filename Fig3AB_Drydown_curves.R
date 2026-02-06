#script for plotting drydown curves + SSM & RWC80-50 regions

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
df_combA <- readRDS(paste0(data_dir,"df_comb_run3+run4ELYRHI.rds"))
df_combB <- readRDS(paste0(data_dir,"df_comb_run4.rds"))
#add field for r-squared display
df_combA$r_squared_display <- sprintf("R^2 == %.3f", df_combA$r_squared)
df_combB$r_squared_display <- sprintf("R^2 == %.3f", df_combB$r_squared)

gmin_modelled_summA <- readRDS(paste0(data_dir,"gmin_modelled_summ_run3+run4ELYRHI.rds"))
gmin_modelled_summB <- readRDS(paste0(data_dir,"gmin_modelled_summ_run4.rds"))

#read in WP, etc. data for run3+run4ELYRHI####
RWC_limits <- read_xlsx(paste0(raw_data_dir,"Arens_RWC_limits_JK+run4ELYRHI.xlsx"))%>%
  select(-`-1/p50`, -`-1/p50.se`, -`-1/p88`, -`-1/p88.se`, -`100-RWC_p50`, -`100-RWC_p88`)
RWC_limits <- head(RWC_limits,4)%>%
  rename(SPP = GENSPP)%>%
  select(SPP, RWC_tlp, RWC_p88)#select only relevant data for plots
#adjoin to gmin_modelled_summA
gmin_modelled_summA <- left_join(gmin_modelled_summA,RWC_limits, by="SPP")

}


####Fig3A####
intercepts <- data.frame(
  Treat = unique(gmin_modelled_summA$Temp),
  GENSPP = unique(gmin_modelled_summA$SPP),
  bin_start = unique(gmin_modelled_summA$RWC_tlp),
  bin_end = unique(gmin_modelled_summA$RWC_p88),
  tlp_label = "RWC_tlp", 
  p88_label = "RWC_p88"
)

ggplot(NULL) +
  geom_rect(
    data = filter(intercepts),
    aes(xmin = -Inf, xmax = Inf, ymin =50, ymax = 80),
    alpha = 1, fill = "#dfdfdf"
  ) +
  geom_hline(
    data = filter(intercepts),
    aes(yintercept = as.numeric(bin_start), linetype = tlp_label),
    color = "black", alpha = 1, linewidth = 0.8
  ) +
  geom_hline(
    data = filter(intercepts),
    aes(yintercept = as.numeric(bin_end), linetype = p88_label),
    color = "black", alpha = 1, linewidth = 0.8
  ) +
  
  geom_point(
    data = filter(df_combA),
    aes(x = Drydown_time/1440, y = RWC, col = Sample_no, shape = Sample_no), # Map shape to Sample_no
    size = 2
  ) +
  
  geom_line(
    data = filter(df_combA),
    aes(x = Drydown_time/1440, y = fit, col = Sample_no)
  ) +
  facet_wrap(
    ~ GENSPP, 
    scales = "free",
    nrow=2,
    labeller = labeller(
      GENSPP = c(
        "CANCON"	= "Cannomois congesta",
        "ERIMON"	= "Erica monsoniana",
        "PROREP"	= "Protea repens",
        "ELYRHI"	= "Dicerothamnus rhinocerotis"
      )
    )
  ) +
  labs(title = "A", 
       x= "Desiccation time (d)", 
       y= "RWC (%)")+
  scale_y_continuous(
    breaks = c(100, 80, 50, 20, 0)
  )+
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() +
  geom_label(
    data = filter(df_combA) %>% 
      distinct(Sample_no, .keep_all = TRUE),
    aes(x = -Inf, y = 100, label = r_squared_display, col = Sample_no),
    hjust = -5.5, vjust = 5.5, size = 2.5,
    position = position_nudge(y = seq(35, 0.05, length.out = n_distinct(df_combA$Sample_no))),
    show.legend = FALSE,
    parse = TRUE,
    label.size = 0
  )+
  scale_linetype_manual(
    name = "RWC Limits",
    values = c("dotted", "dashed"),
    labels = c(bquote(RWC[tlp]), bquote(RWC[Ψ88]))
  ) +
  scale_color_manual(
    name = "Sample No.",
    values = c("#d95f02", "#e6ab02", "#66a61e", "#1b9e77", "#7570b3", "#e7298a"),
    labels = c("1", "2", "3", "4", "5", "6")
  ) +
  scale_shape_manual(
    name = "Sample No.",
    values = c(16, 1, 15, 0, 17, 2), # Custom shapes
    labels = c("1", "2", "3", "4", "5", "6")
  ) +
  guides(
    color = guide_legend(title = "Sample No.", nrow = 1),
    shape = guide_legend(title = "Sample No."),
    linetype = guide_legend(title = "RWC Limits", 
                            nrow = 1,
                            override.aes = list(linetype = c("dashed", "dotted"), color = c("black", "black")))
  ) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "italic"),
    strip.background = element_rect(fill = NA, color = NA),
    legend.position = "bottom",
    legend.box = "horizontal",                    # Stack the legend title and items vertically
    legend.margin = margin(t = -5, unit = "pt")
  )

#save
ggsave("Fig3A_Drydown_curves_run3+run4ELYRHI.png", 
       device = "png", 
       path = fig_dir, 
       width = 20, height = 15, units = "cm", dpi = 300)

####Fig3B####
intercepts2 <- data.frame(
  Temp = unique(gmin_modelled_summB$Temp),
  GENSPP = unique(gmin_modelled_summB$SPP),
  bin_start = 80,
  bin_end = 50
)

#NB: ELYRHI removed!
ggplot(NULL) +
  geom_rect(
    data = filter(intercepts2, GENSPP != "ELYRHI"),
    aes(xmin = -Inf, xmax = Inf, ymin = as.numeric(bin_end), ymax = as.numeric(bin_start)),
    alpha = 1, fill = "#dfdfdf"
  ) +
  geom_point(
    data = filter(df_combB, GENSPP != "ELYRHI"),
    aes(x = Drydown_time/1440, y = RWC, col = Sample_no, shape = Sample_no), # Map shape to Sample_no
    size = 2
  ) +
  
  geom_line(
    data = filter(df_combB, GENSPP != "ELYRHI"),
    aes(x = Drydown_time/1440, y = fit, col = Sample_no)
  ) +
  facet_wrap(
    ~ GENSPP, 
    scales = "free",
    nrow=3,
    labeller = labeller(
      GENSPP = c(
        "AGACAP"	= "Agathosma capensis",
        "ASPSHA"	= "Aspalathus shawii",
        "CLURUB"	= "Clutia rubricaulis",
        "ELYRHI"	= "Dicerothamnus rhinocerotis",
        "LOBDEC"	= "Lobostemon decorus",
        "MICPOL"	= "Microdon polygaloides",
        "OEDSQU"	= "Oedera squarrosa",
        "PASOBT"	= "Passerina obtusifolia",
        "PROLAU"	= "Protea laurifolia",
        "PROLOR"	= "Protea lorifolia",
        "RUSMUL"	= "Ruschia multiflora",
        "SELDOL"	= "Selago dolosa",
        "WAHNOD"	= "Wahlenbergia neorigida"
      )
    )
  ) +
  labs(title = "B", 
       x= "Desiccation time (d)", 
       y= "RWC (%)")+
  scale_y_continuous(
    breaks = c(100, 80, 50, 20, 0)
  )+
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() +
  geom_label(
    data = filter(df_combB, GENSPP != "ELYRHI") %>% 
      distinct(Sample_no, .keep_all = TRUE),
    aes(x = -Inf, y = 100, label = r_squared_display, col = Sample_no),
    hjust = -2.75, vjust = 5.5, size = 2.5,
    position = position_nudge(y = seq(30, 0.05, length.out = n_distinct(df_combB$Sample_no))),
    show.legend = FALSE,
    parse = TRUE,
    linewidth = 0
  )+
  scale_linetype_manual(
    name = "RWC Limits",
    values = c("dashed", "dashed"),
    labels = c(bquote(RWC[80]), bquote(RWC[50]))
  ) +
  scale_color_manual(
    name = "Sample No.",
    values = c("#d95f02", "#e6ab02", "#66a61e", "#1b9e77", "#7570b3", "#e7298a"),
    labels = c("1", "2", "3", "4", "5", "6")
  ) +
  scale_shape_manual(
    name = "Sample No.",
    values = c(16, 1, 15, 0, 17, 2), # Custom shapes
    labels = c("1", "2", "3", "4", "5", "6")
  ) +
  guides(
    color = guide_legend(title = "Sample No.", nrow = 1),
    shape = guide_legend(title = "Sample No.")
    ) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "italic"),
    strip.background = element_rect(fill = NA, color = NA),
    legend.position = "bottom",
    legend.box = "horizontal",                    # Stack the legend title and items vertically
    legend.margin = margin(t = -5, unit = "pt")
  )

#save
ggsave("Fig3B_Drydown_curves_run4.png", 
       device = "png", 
       path = fig_dir, 
       width = 25, height = 25, units = "cm", dpi = 300)
