library(phyloseq)
library(ade4)
library(dplyr)
library(stringr)
library(shiny)
library(shinydashboard)
library(shinyjs)
library(DT)
library(shinywdstar)
library(futile.logger)
library(DiagrammeR)

flog.appender(appender.file("shinywdstar.log"))

source("newutil.R")
source("module-loadPanel.R")
source("module-factorPanel.R")
source("module-testsPanel.R")
source("module-analysisPanel.R")
source("module-statusPanel.R")

distance_choices <- phyloseq::distanceMethodList
data(physeq)

#other data set
data(phy)
data(phenotypes)

#physeq <- readRDS("physeq.rds")
#physeq <- readRDS("physeq.rds")
#load("data/physeq.rda")

