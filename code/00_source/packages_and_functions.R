
# Script to install and load packages, define functions from overlapping packages, and define some custom functions
# For project Urban Bat Diet
# NIOO-KNAW
# By Sander W. J. Buddendorf, 2025-2026

rm(list = ls())

# Packages ----
## Package load/install function
get_packages <- function(list_of_packages) {
  new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  invisible(lapply(list_of_packages, require, character.only = TRUE))
}

## Load some first packages 
get_packages(
  c("tidyverse", "ggplot2", "data.table", "here", "readxl", "glmmTMB", "DHARMa", "car", "emmeans",
    "janitor", "vegan", "showtext", "ggsignif", "BiocManager"))

## Loading Biostring separately
if (!"Biostrings" %in% installed.packages()[,"Package"]) BiocManager::install("Biostrings")
library(Biostrings)

# Define functions because some overlap in different packages
select <- dplyr::select
margin <- ggplot2::margin
rename <- dplyr::rename

# Set sort in count to T by default
count <- function(x, ..., sort = TRUE) {
  dplyr::count(x, ..., sort = sort)
}

# Colours for urban-rural
custom_scale_colour <- scale_color_manual(values = c("Urban" = "#795443",
                                                     "Rural" = "#9cd24d",
                                                     "Natural" = "#9cd24d",
                                                     "Both" = "#5FA5A0"))

custom_scale_fill <- scale_fill_manual(values = c("Urban" = "#795443",
                                                  "Rural" = "#9cd24d",
                                                  "Natural" = "#9cd24d",
                                                  "Both" = "#5FA5A0"))
# Helper function ----
anova_summary <- function(
    models, # named list
    subject # character string, e.g. "Diptera"
) {
  
  anova_list <- eval(parse(text=paste("anova(",paste("models[[",1:length(models),"]]",sep="",collapse=","),")")))
  
  anova_df <- as.data.frame(anova_list)
  anova_df$subject <- subject
  anova_df$formula <- sapply(models, function(x) paste(deparse(formula(x)), collapse = ""))
  anova_df <- anova_df[, c("subject", "formula", names(anova_df)[1:(ncol(anova_df)-2)])]
  
  return(anova_df)
}
