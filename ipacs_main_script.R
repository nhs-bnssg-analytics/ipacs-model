# IPACS model -----------------------------------------------------------------

# Import required libraries
library(doSNOW)
library(foreach)
library(tidyverse, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(parallel)
library(readxl)
library(here)

# Clear workspace and set working directory
rm(list = ls())
wd <- setwd("/home/set_working_directory_here")

# Name of file with model parameters (should be stored in model_inputs folder)
input_filename <- "IPACS_params.xlsx"

# Set manual parameters
nruns <- as.integer(50) # number of runs for each simulation
sd_isr <- 0.5 # initial service/visit rate
sd_esr <- 0.5 # end service rate/final visit rate
temp_seed <- 1
warmup <- 0 #leave as zero unless recoding to not use initial conditions

# Choose method for estimation of mu and sigma (either 1 or 2)
# Use 1 if mean and median are available, but not std dev
# Use 2 if mean and std dev are available
# See set_up.R for the difference between the methods
est_method <- 2

# Run model - P1 (visit) then P2 and P3 (bed) simulation - and print time taken
start_time <- Sys.time()
source(here("functions", "set_up.R"))
source(here("functions", "visit_functions.R"))
source(here("functions", "visit_model.R"))
source(here("functions", "bed_functions.R"))
source(here("functions", "bed_model.R"))
print(difftime(Sys.time(), start_time), quote = FALSE)

# Produce word document report using RMarkdown in outputs/report
rmd_filename <- paste0("ipacs_report_using_",
                       gsub(".xlsx", "", input_filename))
suppressWarnings({
  rmarkdown::render("functions/ipacs_produce_report.Rmd",
                    output_file = here("outputs", "report",rmd_filename))
})
