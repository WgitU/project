#===============================================================
#  Integrative single-cell clustering analysis
# [Wu]
# April, 2018
#===============================================================

#===============================================================
# Reproducibility: MASTER script
#===============================================================

# This file contains instructions for reproducing the data, all
# analyses, and plots contained in the report.

# Download the repository from  GitHub [url blinded].
# The following script assumes the working directory has
# been set to this folder. Except raw data, all of the preprocessed 
# data and code are now available.

# As denoted below, Steps 2 and 3 are computationally 
# intensive and take a *very* long time to run. As such, pre-
# processed data files are available for these steps. Otherwise,
# code within an individual step (e.g., plotting in step 4) 
# assumes code for the previous steps has been run.

#===============================================================
# Step 0: Download the raw data and install necessary packages;
#===============================================================
## Download the raw data and decompress files into "data" folder.

## Necessary packages
library(readr)
library(dplyr)
library(rgeos)
library(reshape2)
library(plotly)
library(ggplot2)
library(Rtsne)
library(factoextra)
library(gridExtra)
library(Seurat)

#===============================================================
# Step 1: Take exploratory data analysis on the raw data.
#===============================================================

setwd("../code")
source("step1.R")

#===============================================================
# Step 2: Conduct data processing/preparation for the analyses.
#===============================================================

setwd("../code")
source("step2.R")

## Alternatively: load pre-processed data
## scRNA-seq data
## M <-read.table("DNAm_pre.txt",header =T)
## DNA methylation data
## X <- read.csv("RNA_500.csv")

#===============================================================
# Step 3: Fit a set of Bayesian models 
#===============================================================

setwd("../code")
source("functions.R")
source("2clus.R")
source("3clus.R")
source("4clus.R")
source("5clus.R")

# Alternatively: load pre-processed data
# load("../data/RNA_DNAm_2.RData")
# load("../data/RNA_DNAm_3.RData")
# load("../data/RNA_DNAm_4.RData")
# load("../data/RNA_DNAm_5.RData")

#===============================================================
# Step 4: Generate all plots in the report
#===============================================================
setwd("../code")
source("step4.R")
