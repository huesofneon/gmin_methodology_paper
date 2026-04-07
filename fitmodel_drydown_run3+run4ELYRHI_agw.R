#script for fitting models to drydown data

### modeifed AGW 19 Marhc 2026


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

# Set paths, values & functions -----------------------------------------------------------------
mainDir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
data_dir <- paste0(mainDir,"/data/")
raw_data_dir <- paste0(data_dir,"raw/")
fig_dir<-paste0(mainDir,"/figures/")

#define atmospheric pressure at height above sea level 
atmosp <- 99891.70/1000#kPa (@120m = approx height of HW Pearson)

#formula for standard error
se <- function(x) {
  return(sd(x) / sqrt(length(x)))
}
}

{
#read in data-----
mass_data5 <- readRDS(paste0(data_dir,"mass_data_run3+run4ELYRHI_agw.rds"))%>%
    na.omit()

treat_names_gmin <- sort(unique(mass_data5$ID_Code))
spp_list <- unique(mass_data5$GENSPP)
temp_list <- unique(mass_data5$Temp)
samples <- unique(mass_data5$Sample_no)
}

#Exponential decay model####
# Setup data frame for coefficients from the model
{
modelcoeffs_gmin <- data.frame(A = 1:1:length(treat_names_gmin), 
                               k = 1:1:length(treat_names_gmin),
                               ID_Code = treat_names_gmin,
                               row.names = treat_names_gmin)

# Modify treatment names and create _
modelcoeffs_gmin <- modelcoeffs_gmin %>%
  mutate(
    SPP = paste0(substr(ID_Code, 10, 15)),
    Temp = substr(ID_Code, 17, 19),
    Sample_no = substr(ID_Code, 21, 28),
    ID_Code = as.factor(ID_Code))

model_formula <- RWC ~ A * exp(-k * Drydown_time) + (100 - A)#Drydown_time in mins
}

#fitting nlsLM() to RWC drydown curves####
for (i in 1:length(treat_names_gmin)) {
  q <- mass_data5 %>%
    filter(ID_Code == treat_names_gmin[i])
  
  if (i == 0)
    {next}#^set to skip problematic samples 
  else{
    k <- 0.001
    A <- 80
    #change k
    
    #exponential decay + translation term to always start at 100
    fitmodel <- nlsLM(model_formula,
                      data = q,
                      start = list(A = A, k = k),
                      control = nls.lm.control(maxiter = 1000, ftol = 1e-5))
    
    
    params <- coef(fitmodel)
    modelcoeffs_gmin[i, 1:2] <- params
    
    q <- q %>%
      mutate(fit = params[1] * exp(-params[2] * (Drydown_time))+(100-params[1])) %>%
      mutate(r_squared = 1 - sum(residuals(fitmodel)^2)/sum((q$RWC - mean(q$RWC))^2))
    
    # Merge data
    if (i == 1) {
      df_comb <- q
    }
    if (i > 1) {
      df_comb <- rbind(df_comb, q)
    }
  }
}

df_comb$r_squared_display <- sprintf("R^2 == %.3f", df_comb$r_squared)

#save files for use in other scripts
saveRDS(modelcoeffs_gmin,paste0(data_dir,"modelcoeffs_gmin_run3+run4ELYRHI.rds"))
saveRDS(df_comb,paste0(data_dir,"df_comb_run3+run4ELYRHI_agw.rds"))

##buffer

####View data-Plot with the exponential decay fits (uncomment to run)####
# ggplot(NULL) +
#   geom_point(data = df_comb, aes(x = Drydown_time, y = RWC, col = Sample_no), size = 2) +
#   geom_line(data = df_comb, aes(x = Drydown_time, y = fit, col = Sample_no)) +
#   facet_wrap(~GENSPP, 
#              scales = "free",
#              nrow=2) +
#   coord_cartesian(ylim = c(NA, 100)) +
#   theme_bw()+
#   geom_text(data = df_comb %>% distinct(Sample_no, .keep_all = TRUE), 
#             aes(x = 80, y = max(RWC, na.rm = TRUE), 
#                 label = r_squared_display, col = Sample_no), 
#             vjust = 6.5, hjust = -4.5, size = 3, 
#             position = position_nudge(y = seq(40, 0.05, length.out = n_distinct(df_comb$Sample_no))),
#             parse = TRUE,
#             show.legend = FALSE)+
#   labs(x= "Dry-down time (min)",
#        y="RWC (%)")+
#   scale_color_discrete (
#     name = "Sample No.",                # legend title
#     labels = c("1","2","3","4","5","6"))  # legend labels
# 