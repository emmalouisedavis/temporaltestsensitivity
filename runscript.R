#' Run script for temporal sensitivity analysis
#' @author Emma Louise Davis
#' @return Model outputs for generating figures in plots.R

setwd("~/OneDrive - Nexus365/BDI/COVID19/TestingModel")

library(purrr)
library(tidyverse)
library(survival)
library(survminer)
library(flexsurv)
library(caret)
library(survMS)
library(STAR)
library(reshape2)
library(gridExtra)

source("functions.R")
source("getparams.R")
source("bootstrap.R")
source("meanvals.R")