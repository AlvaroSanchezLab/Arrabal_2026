

library(readxl)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(scales)
library(tidyverse)
library(grid)
library(reshape2)
library(colorspace)
library(stringr)
library(readr)



# Resolve this script's directory so the relative paths below work anywhere
# (supports both `Rscript` and RStudio "Source").
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0)
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    return(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)))
  return(normalizePath("."))
}
setwd(get_script_dir())

# Folder containing the input csv files (relative to scripts/clean_data/)
input_folder <- "../../data/2_csv_clean_batch_effect/"

# Folder where the new csv files will be saved
output_folder <- "../../data/3_clean_data_module/"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# List all csv files in the input folder
datasets <- list.files(
  path = input_folder,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Process each dataset
for(file in datasets){
  
  # Load dataset
  df <- read.csv(file)
  
  # Detect replicate columns automatically
  rep_cols <- grep("rep", colnames(df), value = TRUE)
  
  # Calculate mean across replicates
  df$mean <- rowMeans(df[, rep_cols], na.rm = TRUE)
  
  # Calculate standard deviation across replicates
  df$sd <- apply(
    df[, rep_cols],
    1,
    sd,
    na.rm = TRUE
  )
  
  # Create output filename
  output_file <- file.path(
    output_folder,
    paste0(
      tools::file_path_sans_ext(basename(file)),
      "_module.csv"
    )
  )
  
  # Save updated dataset
  write.csv(df, output_file, row.names = FALSE)
  
  cat("Processed:", basename(file), "\n")
}


