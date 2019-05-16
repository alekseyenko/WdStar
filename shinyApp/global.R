#if(!requireNamespace("BiocManager")){
#  install.packages("BiocManager")
#}
#BiocManager::install("phyloseq")
#BiocManager::install("Rhdf5lib")
#BiocManager::install("permute")

# https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
list.of.packages <- c(#"phyloseq", 
                      "ade4", "dplyr", "stringr", "shiny", 
                      "shinydashboard",
                      "shinyjs", "DT", "futile.logger")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(phyloseq)
library(ade4)
library(dplyr)
library(stringr)
library(shiny)
library(shinydashboard)
library(shinyjs)
library(DT)
library(futile.logger)

ts <- substr(as.character(Sys.time()), 1, 10)
log_filename <- paste0("shinywdstar_",ts,".log")
flog.appender(appender.file(log_filename))

source("newutil.R")
source("module-loadPanel.R")
source("module-factorPanel.R")

#source("module-testsPanel.R")
source("module-analysisPanel.R")
source("module-statusPanel.R")

distance_choices <- phyloseq::distanceMethodList
data(soilrep)
physeq <- soilrep

