#script for plotting assorted demonstrative figures

# Initialization ---------------------------------------------------------------
rm(list = ls()) # Clear environment

{
  mainDir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
  data_dir <- paste0(mainDir,"/data/")
  raw_data_dir <- paste0(data_dir,"raw/")
  fig_dir<-paste0(mainDir,"/figures/")

#-------------------------------------------------------------------------
library(ggplot2)
library(grid)
library(tibble)
library(dplyr)
library(patchwork)
}

#setup----  
{# Define your parameters
A <- 80   # Example value for A
k <- 10    # Example value for k (increased to ensure the curve fits within the new x range)

# Create a data frame for x values
x_values <- seq(0, 1, length.out = 100)  # X values from 0 to 1
y_values <- A * exp(-k * x_values) + (100 - A)  # Calculate y values based on the equation

# Combine x and y into a data frame
curve_data <- data.frame(x = x_values, y = y_values)

#values for vertical lines Critical WP
x_vals_vlines <- c(
  0.05,  # Ψ[tlp]
  0.15,    # Ψ[50]
  0.275    # Ψ[88]
)
y_vals_vlines <- A * exp(-k * x_vals_vlines) + (100 - A)
vals_vlines <- data.frame(x = x_vals_vlines, y = y_vals_vlines)

# Define annotations
WP_labels <- c(bquote("Ψ"[TLP]),bquote("Ψ"[50]),bquote("Ψ"[88]))

slope_labels <- c(
  bquote(atop("SSM"[Ψ50], "(Ψ"[TLP]~"- Ψ"[50]~")")),
  bquote(atop("SSM"[Ψ88], "(Ψ"[TLP]~"- Ψ"[88]~")")),
  bquote("post-Ψ"[TLP]),
  "no threshold",
  bquote("pre-Ψ"[88])
)

slope_params<- data.frame(
  slope = c(
    "SSMΨ50",
    "SSMΨ88",
    "Ψtlp_only",
    "no_limits",
    "Ψ88_only"
  ),
  colour = c(
    "#7570b3",
    "#66a61e",
    "#1b9e77",
    "#e7298a",
    "#e6ab02"
  ),
  #100%=0, psitlp=0.05, psi50=0.15, psi88=0.275, end=0.8 (not really but for display purposes)
  starting_x = c(
    0.05,
    0.05,
    0.05,
    0,
    0
  ),
  #100%=0, psitlp=0.05, psi50=0.15, psi88=0.275, end=0.8 (not really but for display purposes)
  ending_x = c(
    0.15,
    0.275,
    0.8,
    0.8,
    0.275
  ),
  Ψtlplabel_shade = c("black",  "black",  "black",  "gray50", "gray50"),
  Ψ50label_shade2 = c("black",  "gray50", "gray50", "gray50", "gray50"),
  Ψ88label_shade3 = c("gray50", "black",  "gray50", "gray50", "black" ),
  stringsAsFactors = FALSE
)

#position for labels on RHS
x_info <-0.6
}

####SELECT GRAPH HERE####
#set graph: 1 = SSMΨ50, 2 = SSMΨ88, 3 = Ψtlp_only, 4 = no_limits", 5 = Ψ88_only
# graph_no <- 2 #testing purposes

plot_list <- list()

for (graph_no in 1:5) {
# Calculate x-values for y = 12%, 50%, 88% (scaled to 100) based on the inverse sigmoid for the yellow curve
points <- data.frame(
  x = c(slope_params[graph_no,3], slope_params[graph_no,4]),  
  y = A * exp(-k * c(slope_params[graph_no,3], slope_params[graph_no,4])) + (100 - A)) 

#limit the selected slope
selected_curve_data <- curve_data %>%
  filter(x >= points[1,1] & x <= points[2,1])



##### THE PLOT#####
current_plot <- ggplot() +

  # Add vertical dashed lines with reduced opacity
  geom_segment(data = vals_vlines, aes(x = x, xend = x, y = 0, yend = y), 
               linetype = "dotted", color = "gray50", size = 1, alpha = 0.5) +
  #horizontal dashed lines
  geom_segment(data = vals_vlines, aes(x = 0, xend = x_info+0.02, y = y, yend = y), 
               linetype = "dashed", 
               color = c(slope_params[graph_no,5],slope_params[graph_no,6],slope_params[graph_no,7]), #set color by interval 
               size = 1, alpha = 0.5) +

  # Add annotations at the end of each line
  geom_text(data = vals_vlines, aes(x = x_info+0.02, y = y, label = WP_labels), 
            color = c(slope_params[graph_no,5],slope_params[graph_no,6],slope_params[graph_no,7]), #set color by interval
            hjust = 0, vjust = 0.5, parse = TRUE, size = 5) +

  #Whole RWC line
  geom_line(data = subset(curve_data, x <= 0.8), 
            aes(x = x, y = y), 
            color = "gray75", size = 1.2, linetype = "solid") +
  #Selected RWC line
  geom_line(data = selected_curve_data, aes(x = x, y = y), color = slope_params[graph_no,2], size = 2) +
  # Points on the selected RWC curve
  geom_point(data = points, aes(x = x, y = y), size = 3) +
  
    # Y-axis and X-axis labels
  labs(x = "Desiccation time", y = "RWC") +  
  # Customizing axis limits: y-axis is now 0-100, x-axis still 0-1
  scale_x_continuous(
    limits = c(0, 1), 
    labels = NULL, 
    expand = c(0, 0)
    )+
  scale_y_continuous(
    limits = c(0, 103), 
    labels = NULL, 
    expand = c(0, 0)
    ) +  # y-axis goes from 0 to 103, unmarked
  coord_cartesian(clip = "off") +
  
  
  # slope name
  annotation_custom(grob = textGrob(
    label = as.expression(slope_labels[graph_no]),  # Ensure parsing
    gp = gpar(col = slope_params[graph_no, 2], fontsize = 15), 
    x = x_info+0.03, y = 0.52,  # Position 
    just = "left"
  )) +
  
  #RECTANGULAR SHADED AREA
   annotate("rect", xmin = 0, xmax = x_info+0.01, 
           ymin = points[1,2], ymax = points[2,2], alpha = 0.15, fill = slope_params[graph_no,2]) +  # shaded area
  
  #SELECTED SLOPE LABEL
  annotate("text", x = 0.25, y = 54, label = "Slope interval", color = slope_params[graph_no,2], size = 4) +
 
  annotate("segment",
           x = x_info, xend = x_info,
           y = points[1,2], yend = points[2,2],
           linetype = "dotted", size = 1, color = "black"
  ) +
  
  # Top arrow
  annotate("segment",
           x = x_info, xend = x_info,
           y = points[1,2]-1, yend = points[1,2],
           size = 0, color = "black",
           arrow = arrow(type = "closed", length = unit(0.15, "inches"))
  ) +
  
  # Bottom arrow
  annotate("segment",
           x = x_info, xend = x_info,
           y = points[2,2]+1, yend = points[2,2],
           size = 0, color = "black",
           arrow = arrow(type = "closed", length = unit(0.15, "inches"))
  )  +
  # #arrow for RWC
  # geom_segment(aes(x = 0.79, y = A * exp(-k * 1) + (100 - A), 
  #                  xend = 0.8, yend = A * exp(-k * 1) + (100 - A)),
  #              size = 0,
  #              arrow = arrow(type = "closed", length = unit(0.2, "inches")), color = "gray75") +  # Arrow at the end of the curve
  
  # Adding x- and y-axis lines
  theme_minimal(base_size = 15) +
  theme(
    panel.background = element_rect(fill = "white", color = NA), # White background
    panel.grid = element_blank(),  # No gridlines
    axis.line = element_line(color = "black",  size = 1.2),  # Add x and y axis lines
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.ticks = element_blank()  # No axis ticks
  )
current_plot

plot_list[[graph_no]] <- current_plot

}

wrap_plots(plot_list, nrow = 2)

# row1 <- wrap_plots(plot_list[1:3], nrow = 1)
# row2 <- plot_spacer() | plot_list[[4]] | plot_list[[5]] | plot_spacer()
# 
# row1 / row2

ggsave("FigS1_intervals.png", device = "png", path = fig_dir, width = 35, height = 25, units = "cm", dpi = 300)
