## Title: 'helping functions' for R script analyzing the roosting structure of flying-fox roosts in SE QLD and NE NSW
## Manuscript: Conventional wisdom on roosting behavior of Australian flying-foxes—A critical review, and evaluation using new data <https://doi.org/10.1002/ece3.8079>
## Author: Tamika Lunn, Griffith University
## Version: VSubmission, created 13 November 2021

## V1-7 - manuscript preparation
## VSubmission - Code transferred from V10 "Bat structure_helper functions_V10-revision.R"

rm(list=ls())

##############################################################################################
##---------------------------------- Dependency Management: --------------------------------##
##############################################################################################
library(renv)
## run renv::init() to initialize a new project-local environment with a private R library. Run this in the console only
## after working in the project, run renv::snapshot() to save the state of the project library to the lockfile (called renv.lock). Run this in the console only
## if needed (e.g update or installation of a package breaks code), run renv::restore() to revert to the previous state as encoded in the lockfile. Run this in the console only
## run devtools::session_info() to print exact versions of all the packages, libraries, and software used, the operating system as well as hardware. To save in code_README

##############################################################################################
##---------------------------------- Overview of functions: --------------------------------##
##############################################################################################

## The first tree blocks of code (interquartile ranges, extraction of occupied sites, and functions to wrangle mean & SD values) are for converting raw data input into processed data (used in '01_download-and-clean-data_Vsubmission.R') 
## The fourth block of code (plotting functions) are for visualizing raw data
## The fifth block of code (model outputs) are for interpreting and presenting fitted models. The code gives best ranking models (based on AIC) and extracts output of best ranked models in a styled kable table for presentation in the manuscript

################################################################################################
##-----------------------------------------Start code-----------------------------------------##
################################################################################################

## Load base packages:
library(plyr) #for ddply and rename
library(tidyverse)
library(readr) #for saving kables
library(anchors) #for replace.value
library(HH) #for horizontal stacked bar charts
library(reshape2) #for melt
library(cowplot) #pretty Ggplot
library(dplyr)
library(mosaic) #for more efficient ifelse with mutate
library(knitr)
library(kableExtra)
library(spatstat) #spatial package
library(RColorBrewer) #create custom colour ramps
library(egg) #for ggarrange
library(sp) #for spDistsN1 (calculate euclidian distances)
library(DescTools) #for Overlap() and Interval (), in calculating overlap and interval of species roosting (https://rdrr.io/cran/DescTools/man/overlaps.html)
library(ggpubr) #for functions in ggarrange - scale='free'


##----------------------Function(s): Calculate lower & upper interquartile range ----------------------##
Liqr<-function(x) { 
  return(round(quantile(x,0.25,na.rm=TRUE),4)) 
}
Uiqr<-function(x) { 
  return(round(quantile(x,0.75,na.rm=TRUE),4)) 
}

mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
##------------------------------------------------------------------------------------------------------##


##------------------------------- Function: extract occupied sites only --------------------------------##
## The following are repeated copy-paste until the second filter for species

#########################################
### Occupied by at least one species: ###
#########################################

occ.all <- function(treebat, by) {
  
  occ <- treebat %>%
    filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>%
    ddply(by, summarise,
          BFF.N = sum(BFF.index),
          BFF.occ = max(BFF.index), 
          GHFF.N = sum(GHFF.index),
          GHFF.occ = max(GHFF.index), 
          LRFF.N = sum(LRFF.index),
          LRFF.occ = max(LRFF.index)) %>%
    mutate(BFF.occ = replace(BFF.occ, BFF.occ >0, 1)) %>%
    mutate(GHFF.occ = replace(GHFF.occ, GHFF.occ >0, 1)) %>%
    mutate(LRFF.occ = replace(LRFF.occ, LRFF.occ >0, 1))
  
  occ.all.df <- occ %>%
    filter(BFF.occ>0|GHFF.occ>0|LRFF.occ>0) #occupied by at least one species 
  return(occ.all.df)
}

########################
### Occupied by BFF: ###
########################

occ.BFF <- function(treebat, by) {
  
  occ <- treebat %>%
    filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>%
    ddply(by, summarise,
          BFF.N = sum(BFF.index),
          BFF.occ = max(BFF.index), 
          GHFF.N = sum(GHFF.index),
          GHFF.occ = max(GHFF.index), 
          LRFF.N = sum(LRFF.index),
          LRFF.occ = max(LRFF.index)) %>%
    mutate(BFF.occ = replace(BFF.occ, BFF.occ >0, 1)) %>%
    mutate(GHFF.occ = replace(GHFF.occ, GHFF.occ >0, 1)) %>%
    mutate(LRFF.occ = replace(LRFF.occ, LRFF.occ >0, 1))
  
  occ.BFF.df <- occ %>%
    filter(BFF.occ>0) #occupied by at least one BFF 
  return(occ.BFF.df)
}

#########################
### Occupied by GHFF: ###
#########################

occ.GHFF <- function(treebat, by) {
  
  occ <- treebat %>%
    filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>%
    ddply(by, summarise,
          BFF.N = sum(BFF.index),
          BFF.occ = max(BFF.index), 
          GHFF.N = sum(GHFF.index),
          GHFF.occ = max(GHFF.index), 
          LRFF.N = sum(LRFF.index),
          LRFF.occ = max(LRFF.index)) %>%
    mutate(BFF.occ = replace(BFF.occ, BFF.occ >0, 1)) %>%
    mutate(GHFF.occ = replace(GHFF.occ, GHFF.occ >0, 1)) %>%
    mutate(LRFF.occ = replace(LRFF.occ, LRFF.occ >0, 1))
  
  occ.GHFF.df <- occ %>%
    filter(GHFF.occ>0) #occupied by at least one GHFF 
  return(occ.GHFF.df)
}

#########################
### Occupied by LRFF: ###
#########################

occ.LRFF <- function(treebat, by) {
  
  occ <- treebat %>%
    filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>%
    ddply(by, summarise,
          BFF.N = sum(BFF.index),
          BFF.occ = max(BFF.index), 
          GHFF.N = sum(GHFF.index),
          GHFF.occ = max(GHFF.index), 
          LRFF.N = sum(LRFF.index),
          LRFF.occ = max(LRFF.index)) %>%
    mutate(BFF.occ = replace(BFF.occ, BFF.occ >0, 1)) %>%
    mutate(GHFF.occ = replace(GHFF.occ, GHFF.occ >0, 1)) %>%
    mutate(LRFF.occ = replace(LRFF.occ, LRFF.occ >0, 1))
  
  occ.LRFF.df <- occ %>%
    filter(LRFF.occ>0) #occupied by at least one LRFF 
  return(occ.LRFF.df)
}
##------------------------------------------------------------------------------------------------------##

##-------------------------------- Functions: Wrangle mean & SD values ---------------------------------##
## There are two main data types for bat data: 'index-based' and 'count-based'. Index-based data were taken for all trees (and include index and weighted index values), whereas count was taken for a subset of trees (and include count, min and max height).
## 'pixel-based' data are kernel density estimates, calculated from the spatial point pattern of tree distributions within plots weighted by the 'index-based' data per tree. 
## There are two possible summary levels: plot-level and roost-level. Note - in the published paper, 'plot' is referred to as 'subplot', and 'site' as 'roost'.
## There are two possible calculations for roost-level summaries: with tree as the lowest replicate, and plot as the lowest replicate
## Summary data frames are returned in either a wide format (species in one column, then separate columns for each summary value) OR long format (species in one column, then one column for 'value' and one column for 'measure'). May also be an 'extra-wide' format which reflects the way the data was entered (separate columns for site*measure combinations)
## Wide format for when named columns are to be selected for use. Long format for when values need to be filtered from a single column

### Naming convention:
###      Functions: 'data type' '.' 'LOWEST REPLICATE' '.' 'returned df format'
###      Output: 'data type' '.' 'LOWEST REPLICATE' '.' 'returned df format' '_' 'summary level'
###      data type = index | count
###      LOWEST REPLICATE = TREE | PLOT
###      returned df format = wide |Ewide | long
###      summary level = site | plot

##------------------------------------ Tree as lowest replicate -------------------------------------##
## Functions below give summaries of measures at the plot or site level (specified with 'by')
## For the following, mean and sd are calculated with tree as the lowest replicate. Either across all trees, or occupied trees only ('occtree')

##########################
### Index based values ###
##########################

### Summary in wide format:
index.TREE.Ewide <- function(treebat, by) { #returns data in extra-wide format
  ## calculate summaries per species:
  index.TREE.Ewide_output <- treebat %>%
    filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed
    mutate(all.index.weight = BFF.index.weight+GHFF.index.weight+LRFF.index.weight) %>% #sum WEIGHTS, not index values, to give indication of overall bat load
    ddply(by, summarise, #by = c("site.code", "session", "site.accession", "subplot") for PLOT LEVEL summary
          #BFF
          BFF.N = sum(BFF.index.weight), #total sum of bats across trees
          BFF.occ = sum(!is.na(BFF.index.weight)[BFF.index.weight>0]), #the number of occupied trees
          BFF.density=BFF.N/400,#400 = area of 20x20 meter grid. Number of bats per m squared. Note - can't calculate Uiqr or Liqr for this density measure, as only one replicate per plot
          BFF.max.tree = max(BFF.index.weight), #maximum number of bats in an individual tree
          BFF.mean.tree = mean(BFF.index.weight), #mean of bats per tree across all trees
          BFF.mode.tree = mode(BFF.index.weight), #mode of bat index per tree across all trees
          BFF.sd.tree = sd(BFF.index.weight), #standard deviation of bats per tree across all trees
          BFF.Liqr.tree = Liqr(BFF.index.weight), #lower inter quartile range for bats per tree across all trees
          BFF.Uiqr.tree = Uiqr(BFF.index.weight), #upper inter quartile range for bats per tree across all trees
          BFF.mean.occtree = BFF.N/BFF.occ, #mean of bats per tree across occupied trees
          BFF.mode.occtree = mode(BFF.index.weight[BFF.index.weight>0]), #mode of bat index per tree across occupied trees
          BFF.sd.occtree = sd(BFF.index.weight[BFF.index.weight>0]), #standard deviation of bats per tree across occupied trees
          BFF.Liqr.occtree = Liqr(BFF.index.weight[BFF.index.weight>0]), #lower inter quartile range for bats per tree across occupied trees
          BFF.Uiqr.occtree = Uiqr(BFF.index.weight[BFF.index.weight>0]), #upper inter quartile range for bats per tree across occupied trees
          #GHFF
          GHFF.N = sum(GHFF.index.weight),
          GHFF.occ = sum(!is.na(GHFF.index.weight)[GHFF.index.weight>0]), 
          GHFF.density=GHFF.N/400,
          GHFF.max.tree = max(GHFF.index.weight),
          GHFF.mean.tree = mean(GHFF.index.weight),
          GHFF.mode.tree = mode(GHFF.index.weight), #mode of bat index per tree across all trees
          GHFF.sd.tree = sd(GHFF.index.weight),
          GHFF.Liqr.tree = Liqr(GHFF.index.weight),
          GHFF.Uiqr.tree = Uiqr(GHFF.index.weight),
          GHFF.mean.occtree = GHFF.N/GHFF.occ,
          GHFF.mode.occtree = mode(GHFF.index.weight[GHFF.index.weight>0]), 
          GHFF.sd.occtree = sd(GHFF.index.weight[GHFF.index.weight>0]),
          GHFF.Liqr.occtree = Liqr(GHFF.index.weight[GHFF.index.weight>0]),
          GHFF.Uiqr.occtree = Uiqr(GHFF.index.weight[GHFF.index.weight>0]),
          #LRFF
          LRFF.N = sum(LRFF.index.weight),
          LRFF.occ = sum(!is.na(LRFF.index.weight)[LRFF.index.weight>0]), 
          LRFF.density=LRFF.N/400,
          LRFF.max.tree = max(LRFF.index.weight), 
          LRFF.mean.tree = mean(LRFF.index.weight), 
          LRFF.mode.tree = mode(LRFF.index.weight), #mode of bat index per tree across all trees
          LRFF.sd.tree = sd(LRFF.index.weight),
          LRFF.Liqr.tree = Liqr(LRFF.index.weight),
          LRFF.Uiqr.tree = Uiqr(LRFF.index.weight),
          LRFF.mean.occtree = LRFF.N/LRFF.occ, 
          LRFF.mode.occtree = mode(LRFF.index.weight[LRFF.index.weight>0]), 
          LRFF.sd.occtree = sd(LRFF.index.weight[LRFF.index.weight>0]),
          LRFF.Liqr.occtree = Liqr(LRFF.index.weight[LRFF.index.weight>0]),
          LRFF.Uiqr.occtree = Uiqr(LRFF.index.weight[LRFF.index.weight>0]),
          #total
          tree.count = sum(!is.na(BFF.index.weight)), #total number of trees in plots
          all.N = sum(all.index.weight),
          all.occ = sum(!is.na(session)[BFF.index.weight>0|GHFF.index>0|LRFF.index>0]), #number of trees occupied by at least one bat, regardless of species
          all.density = all.N/400,
          all.max.tree = max(all.index.weight),
          all.mean.tree = mean(all.index.weight), #mean number of all bats per tree across all trees
          all.mode.tree = mode(all.index.weight), #mode of bat index per tree across all trees
          all.sd.tree = sd(all.index.weight),
          all.Liqr.tree = Liqr(all.index.weight),
          all.Uiqr.tree = Uiqr(all.index.weight),
          all.mean.occtree = all.N/all.occ, #mean number of all bats per tree across occupied trees. Max not appropriate, because will be occupying different subsets of trees. Sum also not appropriate because occupy some of the same trees
          all.mode.occtree = mode(all.index.weight[all.index.weight>0]), 
          all.sd.occtree = sd(all.index.weight[all.index.weight>0]),
          all.Liqr.occtree = Liqr(all.index.weight[all.index.weight>0]),
          all.Uiqr.occtree = Uiqr(all.index.weight[all.index.weight>0])) %>% 
    mutate(BFF.mean.tree = replace_na(BFF.mean.tree, 0)) %>%
    mutate(BFF.mean.occtree = replace_na(BFF.mean.occtree, 0)) %>%
    mutate(GHFF.mean.tree = replace_na(GHFF.mean.tree, 0)) %>%
    mutate(GHFF.mean.occtree = replace_na(GHFF.mean.occtree, 0)) %>%
    mutate(LRFF.mean.tree = replace_na(LRFF.mean.tree, 0)) %>%
    mutate(LRFF.mean.occtree = replace_na(LRFF.mean.occtree, 0)) %>% 
    mutate(all.mean.tree = replace_na(all.mean.tree, 0)) %>%
    mutate(all.mean.occtree = replace_na(all.mean.occtree, 0)) %>%
    mutate(BFF.mode.tree = replace_na(BFF.mode.tree, 0)) %>%
    mutate(GHFF.mode.tree = replace_na(GHFF.mode.tree, 0)) %>%
    mutate(LRFF.mode.tree = replace_na(LRFF.mode.tree, 0)) %>%
    mutate(all.mode.tree = replace_na(all.mode.tree, 0)) %>%
    mutate(BFF.mode.occtree = replace_na(BFF.mode.occtree, 0)) %>%
    mutate(GHFF.mode.occtree = replace_na(GHFF.mode.occtree, 0)) %>%
    mutate(LRFF.mode.occtree = replace_na(LRFF.mode.occtree, 0)) %>%
    mutate(all.mode.occtree = replace_na(all.mode.occtree, 0)) %>%
    mutate(BFF.Liqr.occtree = replace_na(BFF.Liqr.occtree, 0)) %>%
    mutate(BFF.Uiqr.occtree = replace_na(BFF.Uiqr.occtree, 0)) %>%
    mutate(GHFF.Liqr.occtree = replace_na(GHFF.Liqr.occtree, 0)) %>%
    mutate(GHFF.Uiqr.occtree = replace_na(GHFF.Uiqr.occtree, 0)) %>%
    mutate(LRFF.Liqr.occtree = replace_na(LRFF.Liqr.occtree, 0)) %>%
    mutate(LRFF.Uiqr.occtree = replace_na(LRFF.Uiqr.occtree, 0)) %>%
    mutate(all.Liqr.occtree = replace_na(all.Liqr.occtree, 0)) %>%
    mutate(all.Uiqr.occtree = replace_na(all.Uiqr.occtree, 0)) %>%
    mutate(BFF.sd.tree = replace_na(BFF.sd.tree, 0)) %>%
    mutate(BFF.sd.occtree = replace_na(BFF.sd.occtree, 0)) %>%
    mutate(GHFF.sd.tree = replace_na(GHFF.sd.tree, 0)) %>%
    mutate(GHFF.sd.occtree = replace_na(GHFF.sd.occtree, 0)) %>%
    mutate(LRFF.sd.tree = replace_na(LRFF.sd.tree, 0)) %>%
    mutate(LRFF.sd.occtree = replace_na(LRFF.sd.occtree, 0)) %>%
    mutate(all.sd.tree = replace_na(all.sd.tree, 0)) %>%
    mutate(all.sd.occtree = replace_na(all.sd.occtree, 0)) %>%
    mutate(BFF.max.tree = replace_na(BFF.max.tree, 0)) %>%
    mutate(GHFF.max.tree = replace_na(GHFF.max.tree, 0)) %>%
    mutate(LRFF.max.tree = replace_na(LRFF.max.tree, 0)) %>%
    mutate(all.max.tree = replace_na(all.max.tree, 0)) 
  return(index.TREE.Ewide_output)
}

### Summary in wide format:
## Reformat so that there is one column to identify species:
index.TREE.wide <- function(index.TREE.Ewide_output, by, list) { #returns data in wide format
  index.TREE.wide_output <- index.TREE.Ewide_output %>%
    melt(id.vars = by, measure.vars = c("BFF.N","BFF.occ","BFF.density","BFF.max.tree","BFF.mean.tree","BFF.mode.tree","BFF.sd.tree","BFF.Liqr.tree","BFF.Uiqr.tree","BFF.mean.occtree","BFF.mode.occtree","BFF.sd.occtree","BFF.Liqr.occtree","BFF.Uiqr.occtree",
                                        "GHFF.N","GHFF.occ","GHFF.density","GHFF.max.tree","GHFF.mean.tree","GHFF.mode.tree","GHFF.sd.tree","GHFF.Liqr.tree","GHFF.Uiqr.tree","GHFF.mean.occtree","GHFF.mode.occtree","GHFF.sd.occtree","GHFF.Liqr.occtree","GHFF.Uiqr.occtree",
                                        "LRFF.N","LRFF.occ","LRFF.density","LRFF.max.tree","LRFF.mean.tree","LRFF.mode.tree","LRFF.sd.tree","LRFF.Liqr.tree","LRFF.Uiqr.tree","LRFF.mean.occtree","LRFF.mode.occtree","LRFF.sd.occtree","LRFF.Liqr.occtree","LRFF.Uiqr.occtree", 
                                        "all.N","all.occ","all.density","all.max.tree","all.mean.tree","all.mode.tree","all.sd.tree","all.Liqr.tree","all.Uiqr.tree","all.mean.occtree","all.mode.occtree","all.sd.occtree","all.Liqr.occtree","all.Uiqr.occtree"
    ),
    variable.name = c("species.value"), value.name="value") %>%
    dplyr::mutate(species = str_extract(species.value, "GHFF|BFF|LRFF|all")) %>%
    dplyr::mutate(valuecat = str_extract(species.value, "N|occ|density|max.tree|mean.tree|mode.tree|sd.tree|Liqr.tree|Uiqr.tree|mean.occtree|mode.occtree|sd.occtree|Liqr.occtree|Uiqr.occtree")) %>%
    dplyr::select(-c(species.value)) %>%
    pivot_wider(names_from = valuecat, values_from = value) %>%
    full_join(list, by = c(by))
  return(index.TREE.wide_output)
}

### Summary in long format (copied from above, without pivot_wider:
index.TREE.long <- function(index.TREE.Ewide_output, by) { #return data in long format
  index.TREE.long_output <- index.TREE.Ewide_output %>%         
    melt(id.vars = c(by, "tree.count"), measure.vars = c("BFF.N", "BFF.occ", "BFF.density", "BFF.mean.tree", "BFF.mode.tree", "BFF.mean.occtree","BFF.mode.occtree", "BFF.Liqr.occtree","BFF.Uiqr.occtree","GHFF.N", "GHFF.occ", "GHFF.density","GHFF.mean.tree","GHFF.mode.tree","GHFF.mean.occtree","GHFF.mode.occtree", "GHFF.Liqr.occtree","GHFF.Uiqr.occtree","LRFF.N", "LRFF.occ", "LRFF.density","LRFF.mean.tree","LRFF.mode.tree","LRFF.mean.occtree","LRFF.mode.occtree","LRFF.Liqr.occtree","LRFF.Uiqr.occtree","all.N", "all.occ", "all.density", "all.mean.tree","all.mode.tree", "all.mean.occtree","all.mode.occtree","all.Liqr.occtree","all.Uiqr.occtree", "BFF.sd.tree", "BFF.sd.occtree","GHFF.sd.tree","GHFF.sd.occtree","LRFF.sd.tree","LRFF.sd.occtree","all.sd.tree","all.sd.occtree","BFF.max.tree","GHFF.max.tree","LRFF.max.tree","all.max.tree"),
         variable.name = c("species.meas"), value.name="value") %>%
    dplyr::mutate(species = str_extract(species.meas, "GHFF|BFF|LRFF|all")) %>%
    dplyr::mutate(measure = str_extract(species.meas, "N|occ|density|mean.tree|mode.tree|mean.occtree|mode.occtree|Liqr.tree|Uiqr.tree|Liqr.occtree|Uiqr.occtree|sd.tree|sd.occtree|max.tree")) %>%
    dplyr::select(-c(species.meas)) 
  return(index.TREE.long_output)
}

##########################
### Count based values ###
##########################

### Extract use-able height values from counts subset
calculate.heightdiff <- function(treebat, by) { #returns data in long format
  heights.subset <- treebat %>% 
    mutate(all.count = BFF.count + GHFF.count + LRFF.count) %>%
    mutate(all.min = pmin(BFF.min, GHFF.min, LRFF.min, na.rm = TRUE)) %>%
    mutate(all.max = pmax(BFF.max, GHFF.max, LRFF.max, na.rm = TRUE)) %>%
    melt(id.vars = by, measure.vars = c("BFF.min","GHFF.min","LRFF.min", "all.min", "BFF.max","GHFF.max","LRFF.max","all.max", "BFF.count", "GHFF.count", "LRFF.count", "all.count"),
         variable.name = c("species.value"), value.name="value") %>%
    dplyr::mutate(species = str_extract(species.value, "GHFF|BFF|LRFF|all")) %>%
    dplyr::mutate(valuecat = str_extract(species.value, "min|max|count")) %>%
    dplyr::select(-c(species.value)) %>%
    spread(valuecat,value) %>%
    filter(!is.na(min)|!is.na(max)) %>%
    mutate(height.diff = max-min)
  #mutate(vert.density = count/height.diff)
  return(heights.subset)
}

### Summary in wide format:
count.TREE.wide <- function(heights.subset, by, list) { #return data in wide format
  ## calculate mean & sd:
  temp <- heights.subset %>%
    filter(!vert.density=="Inf")  %>%
    ddply(c(by, "species"), summarise,
          mean.count = mean(count),
          max.count = max(count),
          mean.max = mean(max), 
          mean.min = mean(min),
          mean.height.diff = mean(height.diff), 
          mean.vert.density = mean(vert.density), 
          sd.count = sd(count),
          sd.max = sd(max), 
          sd.min = sd(min),
          sd.height.diff = sd(height.diff), 
          sd.vert.density = sd(vert.density), 
          sample.size = sum(!is.na(count))
    )
  #add missing sites back in (at the species level)
  
  temp.BFF <- temp %>%
    filter(species == "BFF")
  temp.BFF <- temp.BFF %>%
    full_join(list, by = c(by)) %>%
    mutate(sample.size = replace_na(sample.size, 0)) %>%
    mutate(species = replace_na(species, "BFF"))  %>%
    mutate(mean.count = replace_na(mean.count, 0))  %>%
    mutate(max.count = replace_na(max.count, 0))  %>%
    arrange(site.accession)
  
  temp.GHFF <- temp %>%
    filter(species == "GHFF")
  temp.GHFF <- temp.GHFF %>%
    full_join(list, by = c(by)) %>%
    mutate(sample.size = replace_na(sample.size, 0)) %>%
    mutate(species = replace_na(species, "GHFF"))  %>%
    mutate(mean.count = replace_na(mean.count, 0))  %>%
    mutate(max.count = replace_na(max.count, 0))  %>%
    arrange(site.accession)
  
  temp.LRFF <- temp %>%
    filter(species == "LRFF")
  temp.LRFF <- temp.LRFF %>%
    full_join(list, by = c(by)) %>%
    mutate(sample.size = replace_na(sample.size, 0)) %>%
    mutate(species = replace_na(species, "LRFF"))  %>%
    mutate(mean.count = replace_na(mean.count, 0))  %>%
    mutate(max.count = replace_na(max.count, 0))  %>%
    arrange(site.accession)
  
  temp.all <- temp %>%
    filter(species == "all")
  temp.all <- temp.all %>%
    full_join(list, by = c(by)) %>%
    mutate(sample.size = replace_na(sample.size, 0)) %>%
    mutate(species = replace_na(species, "all"))  %>%
    mutate(mean.count = replace_na(mean.count, 0))  %>%
    mutate(max.count = replace_na(max.count, 0))  %>%
    arrange(site.accession)
  
  count.TREE.wide_output <- rbind(temp.all, temp.BFF, temp.GHFF, temp.LRFF)
  return(count.TREE.wide_output)
}

### Summary in long format:
count.TREE.long <- function(count.TREE.wide_output, by) { #return data in long format
  count.TREE.long_output <- count.TREE.wide_output %>%
    melt(id.vars = c(by, "species", "tree.count"), measure.vars = c("mean.count", "max.count", "mean.max", "mean.min", "mean.height.diff", "mean.vert.density", "sd.count", "sd.max", "sd.min", "sd.height.diff", "sd.vert.density", "sample.size"),
         variable.name = c("valuecat"), value.name="value")
  return(count.TREE.long_output)
}


################################
###### Pixel based values ######
################################

### Transpose kernel density data to wide table format (data already averages by plot)
kernel.TREE.wide <- function(kerneldata,list) {
  kernel.TREE.wide_output <- kerneldata %>%
    dplyr::mutate(site.code = str_extract(site.accession, "D...")) %>%
    dplyr::mutate(session = str_extract(site.accession, "0..")) %>%
    mutate(session=as.factor(as.numeric(session))) %>% #can't pull out 10 with "9|10" and needs to be 1 not 001 to match count data
    melt(id.vars = c('site.accession', "site.code", "session", "subplot"), measure.vars = c("GHFF_Liqr_v","BFF_Liqr_v","LRFF_Liqr_v", "all_Liqr_v", "GHFF_Uiqr_v", "BFF_Uiqr_v", "LRFF_Uiqr_v","all_Uiqr_v", "GHFF_median_v", "BFF_median_v", "LRFF_median_v","all_median_v", "GHFF_mean_v","BFF_mean_v", "LRFF_mean_v","all_mean_v","GHFF_sd_v","BFF_sd_v", "LRFF_sd_v","all_sd_v"),
         variable.name = c("species.meas"), value.name="value") %>%
    dplyr::mutate(species = str_extract(species.meas, "..FF|.FF|.ll")) %>%
    dplyr::mutate(measure = str_extract(species.meas, "min_v|max_v|mean_v|median_v|sd_v|Liqr_v|Uiqr_v")) %>%
    dplyr::mutate(value = value*6.41) %>% #v.density is number of bats per 0.156 x 0.156 meter pixel. * by 6.41 (1/0.156) to get bats per m
    dplyr::select(-c(species.meas)) %>%
    pivot_wider(names_from = measure, values_from = value) %>%
    full_join(list, by = c("site.code", "session", "site.accession", "subplot"))
  return(kernel.TREE.wide_output)
}

### Transpose kernel density data to long table format (data already averages by plot)
kernel.TREE.long <- function(kerneldata, list) {
  kernel.TREE.long_output <- kerneldata %>%
    dplyr::mutate(site.code = str_extract(site.accession, "D...")) %>%
    dplyr::mutate(session = str_extract(site.accession, "0..")) %>%
    mutate(session=as.factor(as.numeric(session))) %>% #can't pull out 10 with "9|10" and needs to be 1 not 001 to match count data
    melt(id.vars = c('site.accession', "site.code", "session", "subplot"), measure.vars = c("GHFF_Liqr_v","BFF_Liqr_v","LRFF_Liqr_v", "all_Liqr_v", "GHFF_Uiqr_v", "BFF_Uiqr_v", "LRFF_Uiqr_v","all_Uiqr_v", "GHFF_median_v", "BFF_median_v", "LRFF_median_v","all_median_v", "GHFF_mean_v","BFF_mean_v", "LRFF_mean_v","all_mean_v","GHFF_sd_v","BFF_sd_v", "LRFF_sd_v","all_sd_v"),
         variable.name = c("species.meas"), value.name="value") %>%
    dplyr::mutate(species = str_extract(species.meas, "..FF|.FF|.ll")) %>%
    dplyr::mutate(measure = str_extract(species.meas, "min_v|max_v|mean_v|median_v|sd_v|Liqr_v|Uiqr_v")) %>%
    dplyr::mutate(value = value*6.41) %>% #v.density is number of bats per 0.156 x 0.156 meter pixel. * by 6.41 (1/0.156) to get bats per m
    dplyr::select(-c(species.meas)) %>%
    full_join(list, by = c("site.code", "session", "site.accession", "subplot"))
  return(kernel.TREE.long_output)
}
##------------------------------------------------------------------------------------------------------##

##------------------------------------ Plot as the lowest replicate ------------------------------------##
## Functions below give summaries of measures at the site level by default (still specified with 'by')
## For the following, mean and sd are calculated with plot as the lowest replicate
## Averages are calculated across occupied plots only. Mean in this case is the mean of the plot means. Standard deviation here is the standard deviation between plots. 
## To note: this is fine for some measures (total sum of bats, number of occupied trees, density of bats, max number of bats - i.e. anything where plot is the lowest possible replicate), and the other measures IF the interest is in the mean and variation between plot replicates (e.g. mean of bats per tree across all or occupied trees - does this vary by plot?).

################################
###### Index based values ######
################################

### Summary in long format 
index.PLOT.long.sp <- function(index.TREE.long_output, by, sp, occ.plots.sp, list) { 
  index.PLOT.long.sp_sp <- index.TREE.long_output %>%
    filter(species ==sp) %>%
    filter(rep %in% occ.plots.sp$rep) %>% #choose occupied subplots only. Want this to be species specific 
    ddply(c(by, "measure"), summarise,
          mean.value = mean(value),
          sd.value = sd(value)) %>% 
    as.data.frame() %>%
    right_join(list, by=c(by, "measure")) %>% 
    mutate(mean.value.x = replace_na(mean.value.x, 0)) %>% 
    mutate(mean.value = mean.value.x) %>% 
    mutate(species = sp) %>% 
    dplyr::select(-c(mean.value.x,mean.value.y))
  #nrow(index.PLOT.long.sp_sp) #should be 8 * 13 * 10 (N & occ & density & mean.tree & max.tree & mean.occtree & Liqr.occtree & Uiqr.occtree & sd.occtree & sd.tree)
  return(index.PLOT.long.sp_sp)
}

index.PLOT.long.merge <- function(index.PLOT.long_BFF, index.PLOT.long_GHFF, index.PLOT.long_LRFF, index.PLOT.long_all) {
  index.PLOT.long_merged <- index.PLOT.long_BFF %>%
    full_join(index.PLOT.long_GHFF, by = c("site.accession", "site.code", "session", "measure", "species", "mean.value", "sd.value", "tree.count")) %>%
    full_join(index.PLOT.long_LRFF, by = c("site.accession", "site.code", "session", "measure", "species", "mean.value","sd.value", "tree.count")) %>%
    full_join(index.PLOT.long_all, by = c("site.accession", "site.code", "session", "measure", "species", "mean.value","sd.value", "tree.count")) %>%
    mutate(species=as.factor(species))
  #nrow(index.PLOT.long_merged) #should be 4160 (1040 * 4)
  return(index.PLOT.long_merged)
}

### Summary in wide format
index.PLOT.wide <- function(index.PLOT.long_merged) {
  
  tempmean <- index.PLOT.long_merged  %>%
    dplyr::select(-c(sd.value)) %>%
    spread(measure, mean.value) %>%
    rename( #relabel to merge with sd values, to generate wide format dataframe
      c("mean.density"= "density",
        "mean.Liqr.occtree" = "Liqr.occtree",
        "mean.max.tree" = "max.tree",
        "mean.mean.occtree" = "mean.occtree",
        "mean.mean.tree" = "mean.tree",
        "mean.N" = "N",
        "mean.occ" = "occ",
        "mean.sd.occtree" = "sd.occtree",
        "mean.sd.tree" = "sd.tree",
        "mean.Uiqr.occtree" = "Uiqr.occtree"))
  
  tempsd <- index.PLOT.long_merged  %>%
    dplyr::select(-c(mean.value)) %>%
    spread(measure, sd.value) %>%
    rename( #relabel to merge with mean values, to generate wide format dataframe
      c("sd.density"="density",
        "sd.Liqr.occtree"="Liqr.occtree", 
        "sd.max.tree" = "max.tree",
        "sd.mean.occtree" = "mean.occtree",
        "sd.mean.tree" = "mean.tree",
        "sd.N" = "N",
        "sd.occ"= "occ",
        "sd.sd.occtree" = "sd.occtree",
        "sd.sd.tree" = "sd.tree",
        "sd.Uiqr.occtree" = "Uiqr.occtree"))
  
  index.PLOT.wide_output <- merge(tempmean, tempsd, by=c("site.accession", "site.code","session","species", "tree.count")) 
  
  return(index.PLOT.wide_output)
}


################################
###### Count based values ######
################################
#note that this code is meant for calculating tree density, so there are no zero values by default. Need to be careful about merging back with full list of site x plot x session combos, and that values are NA and not 0

### Summary in wide format:
## Difference is that we're using the pre-calculated mean and sd values, instead of the raw data. Is still in the same format
count.PLOT.wide <- function(heights.subset, list) { #return data in wide format
  ## calculate mean & sd per plot:
  temp_plot <- heights.subset %>%
    filter(!vert.density=="Inf")  %>%
    ddply(c("site.code", "session", "site.accession", "subplot", "rep", "species"), summarise,
          mean.count = mean(count),
          max.count = max(count),
          mean.max = mean(max), 
          mean.min = mean(min),
          mean.height.diff = mean(height.diff), 
          mean.vert.density = mean(vert.density), 
          sd.count = sd(count),
          sd.max = sd(max), 
          sd.min = sd(min),
          sd.height.diff = sd(height.diff), 
          sd.vert.density = sd(vert.density), 
          sample.size = sum(!is.na(count))
    )
  ## calculate mean & sd per site:
  temp_site <- temp_plot %>% 
    ddply(c("site.code", "session", "site.accession", "species"), summarise,
          sd.count = sd(mean.count), #sd calculations need to be first, because I've re-used names (if second, will calculate from the new mean.count values)
          sd.max = sd(mean.max), 
          sd.min = sd(mean.min),
          sd.height.diff = sd(mean.height.diff), 
          sd.vert.density = sd(mean.vert.density),
          mean.count = mean(mean.count),
          max.count = max(max.count),
          mean.max = mean(mean.max), 
          mean.min = mean(mean.min),
          mean.height.diff = mean(mean.height.diff), 
          mean.vert.density = mean(mean.vert.density),
          sample.size = sum(sample.size)
    )
  
  #add missing sites back in (at the species level)
  temp.BFF <- temp_site %>%
    filter(species == "BFF")
  temp.BFF <- temp.BFF %>%
    full_join(list, by = c("site.accession", "site.code", "session")) %>%
    mutate(sample.size = replace_na(sample.size, 0)) %>%
    mutate(species = replace_na(species, "BFF"))  %>%
    mutate(mean.count = replace_na(mean.count, 0))  %>%
    mutate(mean.max = replace_na(mean.max, 0))  %>%
    mutate(max.count = replace_na(max.count, 0))  %>%
    arrange(site.accession)
  
  temp.GHFF <- temp_site %>%
    filter(species == "GHFF")
  temp.GHFF <- temp.GHFF %>%
    full_join(list, by = c("site.accession", "site.code", "session")) %>%
    mutate(sample.size = replace_na(sample.size, 0)) %>%
    mutate(species = replace_na(species, "GHFF"))  %>%
    mutate(mean.count = replace_na(mean.count, 0))  %>%
    mutate(mean.max = replace_na(mean.max, 0))  %>%
    mutate(max.count = replace_na(max.count, 0))  %>%
    arrange(site.accession)
  
  temp.LRFF <- temp_site %>%
    filter(species == "LRFF")
  temp.LRFF <- temp.LRFF %>%
    full_join(list, by = c("site.accession", "site.code", "session")) %>%
    mutate(sample.size = replace_na(sample.size, 0)) %>%
    mutate(species = replace_na(species, "LRFF"))  %>%
    mutate(mean.count = replace_na(mean.count, 0))  %>%
    mutate(mean.max = replace_na(mean.max, 0))  %>%
    mutate(max.count = replace_na(max.count, 0))  %>%
    arrange(site.accession)
  
  count.PLOT.wide_output <- rbind(temp.BFF, temp.GHFF, temp.LRFF)
  return(count.PLOT.wide_output)
}

### Summary in long format:
count.PLOT.long <- function(count.PLOT.wide_output, by) { #return data in long format
  count.PLOT.long_output <- count.PLOT.wide_output %>%
    melt(id.vars = c(by, "species"), measure.vars = c("mean.count", "max.count", "mean.max", "mean.min", "mean.height.diff", "mean.vert.density", "sd.count", "sd.max", "sd.min", "sd.height.diff", "sd.vert.density", "sample.size"),
         variable.name = c("valuecat"), value.name="value")
  return(count.PLOT.long_output)
} #no difference to count.TREE.long 


################################
###### Pixel based values ######
################################

### Summary in long format:
kernel.PLOT.long <- function(kernel.TREE.long_plot, list, occ.BFF.output, occ.GHFF.output, occ.LRFF.output, occ.all.output) {
  ##split by species to identify occupied plots per species:
  #BFF
  kernel.PLOT.long_output.BFF <- kernel.TREE.long_plot %>%
    unite("rep", c(site.accession, subplot), remove = FALSE, sep = "-") %>%
    filter(species =="BFF") %>%
    filter(rep %in% occ.BFF.output$rep) %>% #choose occupied subolots only. Want this to be species specific 
    ddply(c('site.accession', 'site.code', "session", "measure", "species"), summarise,
          mean.value = mean(value), 
          sd.value = sd(value)) %>%
    as.data.frame() %>%
    right_join(list, by=c("site.accession", "site.code", "session", "measure")) %>%
    mutate(mean.value.x = replace_na(mean.value.x, 0)) %>%
    mutate(mean.value = mean.value.x) %>%
    mutate(species = "BFF") %>%
    dplyr::select(-c(mean.value.x,mean.value.y))
  #nrow(kernel.PLOT.long_output) #should be 520 = 8 x 13 x 5 (mean|median|sd|Liqr|Uiqr)
  #GHFF
  kernel.PLOT.long_output.GHFF <- kernel.TREE.long_plot %>%
    unite("rep", c(site.accession, subplot), remove = FALSE, sep = "-") %>%
    filter(species =="GHFF") %>%
    filter(rep %in% occ.GHFF.output$rep) %>% #choose occupied subolots only. Want this to be species specific 
    ddply(c('site.accession', 'site.code', "session", "measure", "species"), summarise,
          mean.value = mean(value), 
          sd.value = sd(value)) %>%
    as.data.frame() %>%
    right_join(list, by=c("site.accession", "site.code", "session", "measure")) %>%
    mutate(mean.value.x = replace_na(mean.value.x, 0)) %>%
    mutate(mean.value = mean.value.x) %>%
    mutate(species = "GHFF") %>%
    dplyr::select(-c(mean.value.x,mean.value.y))
  # LRFF
  kernel.PLOT.long_output.LRFF <- kernel.TREE.long_plot %>%
    unite("rep", c(site.accession, subplot), remove = FALSE, sep = "-") %>%
    filter(species =="LRFF") %>%
    filter(rep %in% occ.LRFF.output$rep) %>% #choose occupied subolots only. Want this to be species specific 
    ddply(c('site.accession', 'site.code', "session", "measure", "species"), summarise,
          mean.value = mean(value), 
          sd.value = sd(value)) %>%
    as.data.frame() %>%
    right_join(list, by=c("site.accession", "site.code", "session", "measure")) %>%
    mutate(mean.value.x = replace_na(mean.value.x, 0)) %>%
    mutate(mean.value = mean.value.x) %>%
    mutate(species = "LRFF") %>%
    dplyr::select(-c(mean.value.x,mean.value.y))
  #all
  kernel.PLOT.long_output.all <- kernel.TREE.long_plot %>%
    unite("rep", c(site.accession, subplot), remove = FALSE, sep = "-") %>%
    filter(species =="all") %>%
    filter(rep %in% occ.all.output$rep) %>% #choose occupied subolots only. Want this to be species specific 
    ddply(c('site.accession', 'site.code', "session", "measure", "species"), summarise,
          mean.value = mean(value), 
          sd.value = sd(value)) %>%
    as.data.frame() %>%
    right_join(list, by=c("site.accession", "site.code", "session", "measure")) %>%
    mutate(mean.value.x = replace_na(mean.value.x, 0)) %>%
    mutate(mean.value = mean.value.x) %>%
    mutate(species = "all") %>%
    dplyr::select(-c(mean.value.x,mean.value.y))
  #nrow(kernel.PLOT.long_output) #should be 520 = 8 x 13 x 5 (mean|median|sd|Liqr|Uiqr)
  
  kernel.PLOT.long_output <- kernel.PLOT.long_output.BFF %>%
    full_join(kernel.PLOT.long_output.GHFF, by = c("site.accession", "site.code", "session", "measure", "mean.value", "sd.value", "species", "tree.count")) %>%
    full_join(kernel.PLOT.long_output.LRFF, by = c("site.accession", "site.code", "session", "measure", "mean.value", "sd.value", "species", "tree.count")) %>%
    full_join(kernel.PLOT.long_output.all, by = c("site.accession", "site.code", "session", "measure", "mean.value", "sd.value", "species", "tree.count")) %>%
    mutate(species=as.factor(species))
  return(kernel.PLOT.long_output)
}


### Summary in wide format
kernel.PLOT.wide <- function(kernel.PLOT.long_site) {
  tempmean <- kernel.PLOT.long_site  %>%
    dplyr::select(-c(sd.value)) %>%
    spread(measure, mean.value) %>%
    rename( 
      c("mean_Liqr_v"="Liqr_v",
        "mean_mean_v"="mean_v",
        "mean_median_v"="median_v",
        "mean_sd_v"="sd_v",
        "mean_Uiqr_v"="Uiqr_v"))
  
  tempsd <- kernel.PLOT.long_site  %>%
    dplyr::select(-c(mean.value)) %>%
    spread(measure, sd.value) %>%
    rename( 
      c("sd_Liqr_v"="Liqr_v",
        "sd_mean_v"="mean_v",
        "sd_median_v"="median_v",
        "sd_sd_v"="sd_v",
        "sd_Uiqr_v"="Uiqr_v"))
  
  kernel.PLOT.wide_output <- merge(tempmean, tempsd, by=c("site.accession", "site.code","session","species", "tree.count"))
  return(kernel.PLOT.wide_output)
}
##------------------------------------------------------------------------------------------------------##

##----------------------------------------- Functions: plotting ----------------------------------------##

###############################################################################################################
### Create new column 'Site' to order site facets based on roost type (more contemporary to morehistorical) ###
###############################################################################################################

create.Site <- function(data) { #will work within piping functions, don't need to specify data to the function
  newdata <- data %>%
    mutate(Site = ifelse(site.code=="DAVO","01-DAVO", 
                         ifelse(site.code=="DSUN","02-DSUN",
                                ifelse(site.code=="DBUR","03-DBUR",
                                       ifelse(site.code=="DRED","04-DRED",
                                              ifelse(site.code=="DTOW","05-DTOW",
                                                     ifelse(site.code=="DCAN","06-DCAN",
                                                            ifelse(site.code=="DCLU","07-DCLU",
                                                                   "08-DLIS"))))))))
  return(newdata)
}

########################
###### Basic plots #####
########################

### Blank frames ###
plot_blank <- function(data, x, y, ylab, xlab, title, clab, flab) {
  ggplot(data, aes(x=as.numeric(as.character(.data[[x]])))) + 
    theme_bw() +
    background_grid(major="x", colour.major = "grey95")+
    labs(y=ylab, x=xlab, colour=clab, fill=flab)+
    ggtitle(title) +
    theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
          axis.title.x = element_text(size=22),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=22),
          legend.text = element_text(size=20),
          legend.title = element_text(size=20),
          plot.title = element_text(size=22), 
          strip.text = element_text(size = 22)) 
}

### With line:
plot_line <- function(data, x, y, ylab, xlab, title, clab, flab) {
  ggplot(data, aes(x=as.numeric(as.character(.data[[x]])))) + 
    geom_line(aes(y=.data[[y]])) +
    theme_bw() +
    background_grid(major="x", colour.major = "grey95")+
    labs(y=ylab, x=xlab, colour=clab, fill=flab)+
    ggtitle(title) +
    theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
          axis.title.x = element_text(size=22),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=22),
          legend.text = element_text(size=20),
          legend.title = element_text(size=20),
          plot.title = element_text(size=22), 
          strip.text = element_text(size = 22)) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) 
}

### With smoother:
### With coloured lines:
plot_col_smooth <- function(data, x, y, colour, ylab, xlab, title, clab, flab) {
  ggplot(data, aes(x=as.numeric(as.character(.data[[x]])))) + 
    #geom_point(aes(y=.data[[y]], color=.data[[colour]]), size=1) +
    geom_smooth(aes(y=.data[[y]], colour=.data[[colour]]), method = "loess", se = TRUE, linetype = "solid", size=1) +
    theme_bw() +
    background_grid(major="x", colour.major = "grey95")+
    labs(y=ylab, x=xlab, colour=clab, fill=flab)+
    ggtitle(title) +
    theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
          axis.title.x = element_text(size=22),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=22),
          legend.text = element_text(size=20),
          legend.title = element_text(size=20),
          plot.title = element_text(size=22), 
          strip.text = element_text(size = 22)) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) 
}

### With dashed lines:
plot_dash_smooth <- function(data, x, y, dash, ylab, xlab, title, dlab, flab) {
  ggplot(data, aes(x=as.numeric(as.character(.data[[x]])))) + 
    #geom_point(aes(y=.data[[y]], linetype=.data[[dash]]), size=1) +
    stat_smooth(aes(y = .data[[y]], linetype=.data[[dash]]), method = "loess", se = TRUE, size=1, colour="black") +
    theme_bw() +
    background_grid(major="x", colour.major = "grey95")+
    labs(y=ylab, x=xlab, linetype=dlab, fill=flab)+
    ggtitle(title) +
    theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
          axis.title.x = element_text(size=22),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=22),
          legend.text = element_text(size=20),
          legend.title = element_text(size=20),
          plot.title = element_text(size=22), 
          strip.text = element_text(size = 22)) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) 
}

### With black lines:
plot_black_smooth <- function(data, x, y, colour, ylab, xlab, title, clab, flab) {
  ggplot(data, aes(x=as.numeric(as.character(.data[[x]])))) + ##Need the . infront of data[[]] for some reason. Any plot using all.index.weight as the y variable won't plot without this notation 
    #geom_point(aes(y=.data[[y]], color=.data[[colour]]), size=1) +
    geom_smooth(aes(y=.data[[y]]), method = "loess", se = TRUE, linetype = "solid", size=1, colour="black") +
    theme_bw() +
    background_grid(major="x", colour.major = "grey95")+
    labs(y=ylab, x=xlab, colour=clab, fill=flab)+
    ggtitle(title) +
    theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
          axis.title.x = element_text(size=22),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=22),
          legend.text = element_text(size=20),
          legend.title = element_text(size=20),
          plot.title = element_text(size=22), 
          strip.text = element_text(size = 22)) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) 
}

### Grouped smooth (e.g separate lines for core and peripheral groups):
plot_group_smooth <- function(data, x, y, colour, linegroup, ylab, xlab, title, clab, flab, slab) {
  ggplot(data, aes(x=as.numeric(as.character(.data[[x]])))) + 
    #geom_point(aes(y=.data[[y]], color=.data[[colour]], shape=.data[[linegroup]]), size=1.5) +
    geom_smooth(aes(y=.data[[y]], colour=.data[[linegroup]]), method = "loess", se = TRUE, linetype = "solid", size=1) +
    theme_bw() +
    background_grid(major="x", colour.major = "grey95")+
    labs(y=ylab, x=xlab, colour=clab, fill=flab, shape=slab)+
    ggtitle(title) +
    theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
          axis.title.x = element_text(size=22),
          axis.text.y = element_text(size=20),
          axis.title.y = element_text(size=22),
          legend.text = element_text(size=20),
          legend.title = element_text(size=20),
          plot.title = element_text(size=22), 
          strip.text = element_text(size = 22)) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) 
}


####################################
### Extract data from loess fits ###
####################################
## These funcions are not too flexible - they assume input of x=session, facet=Site, colour=species
extract_loess_fac <- function(plot) {
  gg1 <- ggplot_build(plot)
  df2 <- data.frame(session = gg1$data[[1]]$x,
                    min = gg1$data[[1]]$y,
                    max = gg1$data[[2]]$y, 
                    Site=gg1$data[[1]]$PANEL, 
                    species=gg1$data[[1]]$group) #[1] refers to min smoothing and [2] to max smoothing
  ## re-label facet groups:
  df2$Site <- ifelse(df2$Site == 1, "01-DAVO", # needs to be site as the input (not site.code) so the number order is correct for conversion. 1 = facet 1, etc. 
                     ifelse(df2$Site == 2, "02-DSUN",
                            ifelse(df2$Site == 3, "03-DBUR",
                                   ifelse(df2$Site == 4, "04-DRED",
                                          ifelse(df2$Site == 5, "05-DTOW",
                                                 ifelse(df2$Site == 6, "06-DCAN",
                                                        ifelse(df2$Site == 7, "07-DCLU",
                                                               "08-DLIS")))))))
  
  # re-label colour groups:
  df2$species <- ifelse(df2$species == 1, "BFF",
                        ifelse(df2$species == 2, "GHFF",
                               "LRFF"))
  return(df2)
}

extract_loess_all <- function(plot) {
  gg1 <- ggplot_build(plot)
  df2 <- data.frame(session = gg1$data[[1]]$x,
                    min = gg1$data[[1]]$y,
                    max = gg1$data[[2]]$y, 
                    Site=gg1$data[[1]]$PANEL, 
                    species=gg1$data[[1]]$group) #[1] refers to min smoothing and [2] to max smoothing
  ## re-label colour groups:
  df2$species <- ifelse(df2$species == 1, "BFF",
                        ifelse(df2$species == 2, "GHFF",
                               "LRFF"))
  return(df2)
}
##------------------------------------------------------------------------------------------------------##

##------------------------------- Functions: model outputs -------------------------------#

## Custom function to return values when available, and * when an error is given
converged <- function(x) { 
  value <- try(x, silent = TRUE) #If there is a p value present (i.e. the model converged) give the value. If not (i.e. an error is returned) do nothing
  output <- ifelse(is.numeric(value)==TRUE, paste(c(round(value, digits=3))), "*") #if p value is present (a numeric) will return the rounded p value, if no value is present (the error) will return an "*
  return(noquote(output)) #return output with no quotes
}

#### Model Comparison: ####
## Function evaluates the AIC value of fitted models and ranks from most to least parsimonious. Output is a styled table for direct input into manuscript or appendices
gamm_modelcomp <- function(modlist, modlist_str, title) { 
  nmodels <- length(modlist)
  modelcomp_output <- as.data.frame(array(dim=c(nmodels,4)))
  colnames(modelcomp_output) <- c("Model structure", "R2", "AIC", "Delta AIC")
  
  for(i in 1:length(modlist)){
    modelcomp_output[i,1] <- noquote(modlist_str[[i]])
    modelcomp_output[i,3] <- round(AIC(summary(modlist[[i]]$lme)), digits=1)
    modelcomp_output[i,2] <- round(summary(modlist[[i]]$gam)[["r.sq"]], digits=3) #coefficient of determination, R2, gives the percentage variation in y explained by x-variables. The range is 0 to 1 (i.e. 0% to 100% of the variation in y can be explained by the x-variables)
  }
  modelcomp_output <- modelcomp_output[order(modelcomp_output$AIC),]
  for(i in 1:length(modlist)){
    modelcomp_output[i,4] <- modelcomp_output[i,3] - modelcomp_output[1,3]
  }
  
  ## Create table
  table_output <- kable(modelcomp_output[,1:4], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F)
  #column_spec(c(3:11), width = "4cm")  #note - collumn_spec slow for long tables
  #scroll_box(height = "500px")
  #footnote(general = "* Indicates models that did not converge")
  #save_kable(paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}

#### Model outputs ####
## Function save model output of best ranked model, where best model includes two variables with interaction term. Output is a styled table for direct input into manuscript or appendices
## For plot-level and tree-level models only
save.output.gamm_interaction <- function(modelname, var1, var2, interact, title) {
  output <- as.data.frame(array(dim=c(7,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- var2
  output[4,"Variable"] <- interact
  output[5,"Variable"] <- "Session"
  output[6,"Variable"] <- "Site"
  output[7,"Variable"] <- "Subplot"
  output[c(1:4), 1] <- "Fixed"
  output[c(5:7), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[2]], digits=3), 
                                   ")", 
                                   sep="")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  
  output[3,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[3]], digits=3)), #Fixed effect value (third level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[3]], digits=3), 
                                   ")", 
                                   sep="")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[3]], digits=4)))
  output[3,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[3]], digits=4)))
  
  output[4,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[4]], digits=3)), #Fixed effect value (fourth level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[4]], digits=3), 
                                   ")", 
                                   sep="")
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[4]], digits=4)))
  output[4,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[4]], digits=4)))
  output[1:4,"F value"] <- paste("-")
  
  ## Random effects
  output[5:7,"Coef (? se)"] <- paste("-")
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[6,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[7,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[5,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[6,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  output[7,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[5:7,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 4) %>%
    pack_rows("Random effects", 5, 7)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}

## Function save model output of best ranked model, where best model includes two variables with interaction term. Output is a styled table for direct input into manuscript or appendices
## For roost-level models only
save.output.gamm_interaction_ROOST <- function(modelname, var1, var2, interact, title) {
  output <- as.data.frame(array(dim=c(6,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- var2
  output[4,"Variable"] <- interact
  output[5,"Variable"] <- "Session"
  output[6,"Variable"] <- "Site"
  #output[7,"Variable"] <- "Subplot"
  output[c(1:4), 1] <- "Fixed"
  output[c(5:6), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[2]], digits=3), 
                                   ")", 
                                   sep="")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  
  output[3,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[3]], digits=3)), #Fixed effect value (third level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[3]], digits=3), 
                                   ")", 
                                   sep="")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[3]], digits=4)))
  output[3,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[3]], digits=4)))
  
  output[4,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[4]], digits=3)), #Fixed effect value (fourth level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[4]], digits=3), 
                                   ")", 
                                   sep="")
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[4]], digits=4)))
  output[4,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[4]], digits=4)))
  output[1:4,"F value"] <- paste("-")
  
  ## Random effects
  output[5:6,"Coef (? se)"] <- paste("-")
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[6,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  #output[7,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[5,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[6,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  #output[7,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[5:6,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 4) %>%
    pack_rows("Random effects", 5, 6)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
  #return(output)
}


## Function save model output of best ranked model, where best model includes two variables with interaction term plus an additional variable (i.e. the full model). Output is a styled table for direct input into manuscript or appendices
## For plot-level and tree-level models only
save.output.gamm_interactionfull <- function(modelname, var1, var2, interact, var3, title) {
  output <- as.data.frame(array(dim=c(8,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- var2
  output[4,"Variable"] <- var3
  output[5,"Variable"] <- interact
  output[6,"Variable"] <- "Session"
  output[7,"Variable"] <- "Site"
  output[8,"Variable"] <- "Subplot"
  output[c(1:5), 1] <- "Fixed"
  output[c(6:8), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[2]], digits=3), 
                                   ")", 
                                   sep="")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  
  output[3,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[3]], digits=3)), #Fixed effect value (third level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[3]], digits=3), 
                                   ")", 
                                   sep="")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[3]], digits=4)))
  output[3,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[3]], digits=4)))
  
  output[4,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[4]], digits=3)), #Fixed effect value (fourth level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[4]], digits=3), 
                                   ")", 
                                   sep="")
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[4]], digits=4)))
  output[4,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[4]], digits=4)))
  
  output[5,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[5]], digits=3)), #Fixed effect value (fifth level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[5]], digits=3), 
                                   ")", 
                                   sep="")
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[5]], digits=4)))
  output[5,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[5]], digits=4)))
  output[1:5,"F value"] <- paste("-")
  
  output[1:5,"F value"] <- paste("-")
  
  ## Random effects
  output[6:8,"Coef (? se)"] <- paste("-")
  output[6,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[7,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[8,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  
  output[6,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[7,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  output[8,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[6:8,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 5) %>%
    pack_rows("Random effects", 6, 8)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}

## Function save model output of best ranked model, where best model includes two variables with interaction term plus an additional variable (i.e. the full model). Output is a styled table for direct input into manuscript or appendices
## For roost-level models only
save.output.gamm_interactionfull_ROOST <- function(modelname, var1, var2, interact, var3, title) {
  output <- as.data.frame(array(dim=c(7,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- var2
  output[4,"Variable"] <- var3
  output[5,"Variable"] <- interact
  output[6,"Variable"] <- "Session"
  output[7,"Variable"] <- "Site"
  #output[8,"Variable"] <- "Subplot"
  output[c(1:5), 1] <- "Fixed"
  output[c(6:7), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[2]], digits=3), 
                                   ")", 
                                   sep="")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  
  output[3,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[3]], digits=3)), #Fixed effect value (third level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[3]], digits=3), 
                                   ")", 
                                   sep="")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[3]], digits=4)))
  output[3,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[3]], digits=4)))
  
  output[4,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[4]], digits=3)), #Fixed effect value (fourth level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[4]], digits=3), 
                                   ")", 
                                   sep="")
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[4]], digits=4)))
  output[4,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[4]], digits=4)))
  
  output[5,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[5]], digits=3)), #Fixed effect value (fifth level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[5]], digits=3), 
                                   ")", 
                                   sep="")
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[5]], digits=4)))
  output[5,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[5]], digits=4)))
  output[1:5,"F value"] <- paste("-")
  
  output[1:5,"F value"] <- paste("-")
  
  ## Random effects
  output[6:7,"Coef (? se)"] <- paste("-")
  output[6,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[7,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  #output[8,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  
  output[6,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[7,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  #output[8,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[6:7,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 5) %>%
    pack_rows("Random effects", 6, 7)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}

## Function save model output of best ranked model, where best model includes one variable only. Output is a styled table for direct input into manuscript or appendices
## For plot-level and tree-level models only
save.output.gamm_single <- function(modelname, var1, title) {
  output <- as.data.frame(array(dim=c(5,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- "Session"
  output[4,"Variable"] <- "Site"
  output[5,"Variable"] <- "Subplot"
  output[c(1:2), 1] <- "Fixed"
  output[c(3:5), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[2]], digits=3), 
                                   ")", 
                                   sep="")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  output[1:2,"F value"] <- paste("-")
  
  
  ## Random effects
  output[3:5,"Coef (? se)"] <- paste("-")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[3,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[4,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  output[5,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[3:5,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 2) %>%
    pack_rows("Random effects", 3, 5)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}

## Function save model output of best ranked model, where best model includes one variable only. Output is a styled table for direct input into manuscript or appendices
## For roost-level models only
save.output.gamm_single_ROOST <- function(modelname, var1, title) {
  output <- as.data.frame(array(dim=c(4,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- "Session"
  output[4,"Variable"] <- "Site"
  #output[5,"Variable"] <- "Subplot"
  output[c(1:2), 1] <- "Fixed"
  output[c(3:4), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[2]], digits=3), 
                                   ")", 
                                   sep="")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  output[1:2,"F value"] <- paste("-")
  
  
  ## Random effects
  output[3:4,"Coef (? se)"] <- paste("-")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  #output[5,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[3,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[4,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  #output[5,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[3:4,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 2) %>%
    pack_rows("Random effects", 3, 4)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}

## Function save model output of best ranked model, where best model includes no variables (i.e. null model). Output is a styled table for direct input into manuscript or appendices
## For plot-level and tree-level models only
save.output.gamm_null <- function(modelname, title) {
  output <- as.data.frame(array(dim=c(4,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- "Session"
  output[3,"Variable"] <- "Site"
  output[4,"Variable"] <- "Subplot"
  output[c(1), 1] <- "Fixed"
  output[c(2:4), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  ## Random effects
  output[2:4,"Coef (? se)"] <- paste("-")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[2,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[3,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  output[4,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[2:4,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 1) %>%
    pack_rows("Random effects", 2, 4)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}

## Function save model output of best ranked model, where best model includes no variables (i.e. null model). Output is a styled table for direct input into manuscript or appendices
## For roost-level models only
save.output.gamm_null_ROOST <- function(modelname, title) {
  output <- as.data.frame(array(dim=c(3,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- "Session"
  output[3,"Variable"] <- "Site"
  #output[4,"Variable"] <- "Subplot"
  output[c(1), 1] <- "Fixed"
  output[c(2:3), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  ## Random effects
  output[2:3,"Coef (? se)"] <- paste("-")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  #output[4,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[2,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[3,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  #output[4,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[2:3,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 1) %>%
    pack_rows("Random effects", 2, 3)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}

############################################################################################################
################################ Save as df, not kable #####################################################
############################################################################################################

## As above but output saved as a dataframe instead of a styled kable table

save.output.gamm_interaction_ROOST_df <- function(modelname, var1, var2, interact, title) {
  output <- as.data.frame(array(dim=c(6,7)))
  colnames(output) <- c("Effect","Variable", "Coef", "se", "t", "F", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- var2
  output[4,"Variable"] <- interact
  output[5,"Variable"] <- "z Session"
  output[6,"Variable"] <- "z Site"
  #output[7,"Variable"] <- "z Subplot"
  output[c(1:4), 1] <- "Fixed"
  output[c(5:6), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3))) #Fixed effect intercept
  output[1,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[1]], digits=3)))
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)))
  output[2,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[2]], digits=3)))
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  
  output[3,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[3]], digits=3)))
  output[3,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[3]], digits=3)))
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[3]], digits=4)))
  output[3,"t"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[3]], digits=4)))
  
  output[4,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[4]], digits=3)))
  output[4,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[4]], digits=3)))
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[4]], digits=4)))
  output[4,"t"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[4]], digits=4)))
  output[1:4,"F"] <- paste("-")
  
  ## Random effects
  output[5:6,"Coef (? se)"] <- paste("-")
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[6,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  #output[7,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[5,"F"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[6,"F"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  #output[7,"F"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[5:6,"t"] <- paste("-")
  
  output$Coef <- as.numeric(output$Coef)
  output$se <- as.numeric(output$se)
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 4) %>%
    pack_rows("Random effects", 5, 6)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  #return(table_output)
  return(output)
}


save.output.gamm_interactionfull_df <- function(modelname, var1, var2, interact, var3, title) {
  output <- as.data.frame(array(dim=c(8,7)))
  colnames(output) <- c("Effect","Variable", "Coef", "se", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- var2
  output[4,"Variable"] <- var3
  output[5,"Variable"] <- interact
  output[6,"Variable"] <- "z Session"
  output[7,"Variable"] <- "z Site"
  output[8,"Variable"] <- "z Subplot"
  output[c(1:5), 1] <- "Fixed"
  output[c(6:8), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)))
  output[1,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[1]], digits=3)))
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)))
  output[2,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[2]], digits=3)))                          
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  
  output[3,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[3]], digits=3)))
  output[3,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[3]], digits=3)))      
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[3]], digits=4)))
  output[3,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[3]], digits=4)))
  
  output[4,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[4]], digits=3)))
  output[4,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[4]], digits=3))) 
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[4]], digits=4)))
  output[4,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[4]], digits=4)))
  
  output[5,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[5]], digits=3)))
  output[5,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[5]], digits=3))) 
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[5]], digits=4)))
  output[5,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[5]], digits=4)))
  output[1:5,"F value"] <- paste("-")
  
  output[1:5,"F value"] <- paste("-")
  
  ## Random effects
  output[6:8,"Coef (? se)"] <- paste("-")
  output[6,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[7,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[8,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  
  output[6,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[7,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  output[8,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[6:8,"t value"] <- paste("-")
  
  output$Coef <- as.numeric(output$Coef)
  output$se <- as.numeric(output$se)
  
  return(output)
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 5) %>%
    pack_rows("Random effects", 6, 8)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  #return(table_output)
}


save.output.gamm_null_df <- function(modelname, title) {
  output <- as.data.frame(array(dim=c(4,7)))
  colnames(output) <- c("Effect","Variable", "Coef", "se", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- "z Session"
  output[3,"Variable"] <- "z Site"
  output[4,"Variable"] <- "z Subplot"
  output[c(1), 1] <- "Fixed"
  output[c(2:4), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)))
  output[1,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[1]], digits=3)))
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  ## Random effects
  output[2:4,"Coef (? se)"] <- paste("-")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[2,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[3,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  output[4,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[2:4,"t value"] <- paste("-")
  
  output$Coef <- as.numeric(output$Coef)
  output$se <- as.numeric(output$se)
  
  return(output)
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 1) %>%
    pack_rows("Random effects", 2, 4)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  #return(table_output)
}


save.output.gamm_single_df <- function(modelname, var1, title) {
  output <- as.data.frame(array(dim=c(5,7)))
  colnames(output) <- c("Effect","Variable", "Coef", "se", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- "z Session"
  output[4,"Variable"] <- "z Site"
  output[5,"Variable"] <- "z Subplot"
  output[c(1:2), 1] <- "Fixed"
  output[c(3:5), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)))
  output[1,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[1]], digits=3)))
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)))
  output[2,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[2]], digits=3)))
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  output[1:2,"F value"] <- paste("-")
  
  
  ## Random effects
  output[3:5,"Coef"] <- paste("-")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[3,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[4,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  output[5,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[3:5,"t value"] <- paste("-")
  
  output$Coef <- as.numeric(output$Coef)
  output$se <- as.numeric(output$se)
  
  return(output)
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 2) %>%
    pack_rows("Random effects", 3, 5)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  #return(table_output)
}

#################################################################################################################
#############################################  GAMs  ############################################################
#################################################################################################################

save.output.gamm_2lev <- function(modelname,factorname, title) {
  output <- as.data.frame(array(dim=c(5,4)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- factorname
  output[3,"Variable"] <- "Session"
  output[4,"Variable"] <- "Site"
  output[5,"Variable"] <- "Subplot"
  output[c(1,2), 1] <- "Fixed"
  output[c(3:5), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[2,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[2]], digits=3), 
                                   ")", 
                                   sep="")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  
  ## Random effects
  output[3:5,"Coef (? se)"] <- paste("-")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  
  
  ## Create table
  table_output <- kable(output[,2:4], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 2) %>%
    pack_rows("Random effects", 3, 5)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}


save.output.gamm_distance <- function(modelname, title) {
  output <- as.data.frame(array(dim=c(9,4)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- "Distance from center"
  output[3,"Variable"] <- "GFF"
  output[4,"Variable"] <- "LRFF"
  output[5,"Variable"] <- "GHFF * Distance from center"
  output[6,"Variable"] <- "LRFF * Distance from center"
  
  output[7,"Variable"] <- "Session"
  output[8,"Variable"] <- "Site"
  output[9,"Variable"] <- "Subplot"
  output[c(1:6), 1] <- "Fixed"
  output[c(7:9), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=1)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=1), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[2,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=1)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[2]], digits=1), 
                                   ")", 
                                   sep="")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[3,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[3]], digits=1)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[3]], digits=1), 
                                   ")", 
                                   sep="")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[3]], digits=4)))
  output[4,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[4]], digits=1)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[4]], digits=1), 
                                   ")", 
                                   sep="")
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[4]], digits=4)))
  output[5,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[5]], digits=1)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[5]], digits=1), 
                                   ")", 
                                   sep="")
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[5]], digits=4)))
  output[6,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[6]], digits=1)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[6]], digits=1), 
                                   ")", 
                                   sep="")
  output[6,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[6]], digits=4)))
  
  ## Random effects
  output[7:9,"Coef (? se)"] <- paste("-")
  output[7,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[8,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[9,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  
  
  ## Create table
  table_output <- kable(output[,2:4], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 7) %>%
    pack_rows("Random effects", 7, 9)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}



#### Model Comparison: ####
gamm_modelcomp <- function(modlist, modlist_str, title) { 
  nmodels <- length(modlist)
  modelcomp_output <- as.data.frame(array(dim=c(nmodels,4)))
  colnames(modelcomp_output) <- c("Model structure", "R2", "AIC", "Delta AIC")
  
  for(i in 1:length(modlist)){
    modelcomp_output[i,1] <- noquote(modlist_str[[i]])
    modelcomp_output[i,3] <- round(AIC(summary(modlist[[i]]$lme)), digits=1)
    modelcomp_output[i,2] <- round(summary(modlist[[i]]$gam)[["r.sq"]], digits=3) #coefficient of determination, R2, gives the percentage variation in y explained by x-variables. The range is 0 to 1 (i.e. 0% to 100% of the variation in y can be explained by the x-variables)
  }
  modelcomp_output <- modelcomp_output[order(modelcomp_output$AIC),]
  for(i in 1:length(modlist)){
    modelcomp_output[i,4] <- modelcomp_output[i,3] - modelcomp_output[1,3]
  }
  
  ## Create table
  table_output <- kable(modelcomp_output[,1:4], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F)
  #column_spec(c(3:11), width = "4cm")  #note - collumn_spec slow for long tables
  #scroll_box(height = "500px")
  #footnote(general = "* Indicates models that did not converge")
  #save_kable(paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}


#### Model outputs ####
save.output.gamm_interaction <- function(modelname, var1, var2, interact, title) {
  output <- as.data.frame(array(dim=c(7,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- var2
  output[4,"Variable"] <- interact
  output[5,"Variable"] <- "Session"
  output[6,"Variable"] <- "Site"
  output[7,"Variable"] <- "Subplot"
  output[c(1:4), 1] <- "Fixed"
  output[c(5:7), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[2]], digits=3), 
                                   ")", 
                                   sep="")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  
  output[3,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[3]], digits=3)), #Fixed effect value (third level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[3]], digits=3), 
                                   ")", 
                                   sep="")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[3]], digits=4)))
  output[3,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[3]], digits=4)))
  
  output[4,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[4]], digits=3)), #Fixed effect value (fourth level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[4]], digits=3), 
                                   ")", 
                                   sep="")
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[4]], digits=4)))
  output[4,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[4]], digits=4)))
  output[1:4,"F value"] <- paste("-")
  
  ## Random effects
  output[5:7,"Coef (? se)"] <- paste("-")
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[6,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[7,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[5,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[6,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  output[7,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[5:7,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 4) %>%
    pack_rows("Random effects", 5, 7)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}

save.output.gamm_interaction_ROOST <- function(modelname, var1, var2, interact, title) {
  output <- as.data.frame(array(dim=c(6,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- var2
  output[4,"Variable"] <- interact
  output[5,"Variable"] <- "Session"
  output[6,"Variable"] <- "Site"
  #output[7,"Variable"] <- "Subplot"
  output[c(1:4), 1] <- "Fixed"
  output[c(5:6), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[2]], digits=3), 
                                   ")", 
                                   sep="")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  
  output[3,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[3]], digits=3)), #Fixed effect value (third level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[3]], digits=3), 
                                   ")", 
                                   sep="")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[3]], digits=4)))
  output[3,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[3]], digits=4)))
  
  output[4,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[4]], digits=3)), #Fixed effect value (fourth level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[4]], digits=3), 
                                   ")", 
                                   sep="")
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[4]], digits=4)))
  output[4,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[4]], digits=4)))
  output[1:4,"F value"] <- paste("-")
  
  ## Random effects
  output[5:6,"Coef (? se)"] <- paste("-")
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[6,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  #output[7,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[5,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[6,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  #output[7,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[5:6,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 4) %>%
    pack_rows("Random effects", 5, 6)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
  #return(output)
}



save.output.gamm_interactionfull <- function(modelname, var1, var2, interact, var3, title) {
  output <- as.data.frame(array(dim=c(8,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- var2
  output[4,"Variable"] <- var3
  output[5,"Variable"] <- interact
  output[6,"Variable"] <- "Session"
  output[7,"Variable"] <- "Site"
  output[8,"Variable"] <- "Subplot"
  output[c(1:5), 1] <- "Fixed"
  output[c(6:8), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[2]], digits=3), 
                                   ")", 
                                   sep="")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  
  output[3,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[3]], digits=3)), #Fixed effect value (third level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[3]], digits=3), 
                                   ")", 
                                   sep="")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[3]], digits=4)))
  output[3,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[3]], digits=4)))
  
  output[4,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[4]], digits=3)), #Fixed effect value (fourth level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[4]], digits=3), 
                                   ")", 
                                   sep="")
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[4]], digits=4)))
  output[4,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[4]], digits=4)))
  
  output[5,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[5]], digits=3)), #Fixed effect value (fifth level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[5]], digits=3), 
                                   ")", 
                                   sep="")
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[5]], digits=4)))
  output[5,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[5]], digits=4)))
  output[1:5,"F value"] <- paste("-")
  
  output[1:5,"F value"] <- paste("-")
  
  ## Random effects
  output[6:8,"Coef (? se)"] <- paste("-")
  output[6,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[7,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[8,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  
  output[6,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[7,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  output[8,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[6:8,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 5) %>%
    pack_rows("Random effects", 6, 8)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}


save.output.gamm_interactionfull_ROOST <- function(modelname, var1, var2, interact, var3, title) {
  output <- as.data.frame(array(dim=c(7,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- var2
  output[4,"Variable"] <- var3
  output[5,"Variable"] <- interact
  output[6,"Variable"] <- "Session"
  output[7,"Variable"] <- "Site"
  #output[8,"Variable"] <- "Subplot"
  output[c(1:5), 1] <- "Fixed"
  output[c(6:7), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[2]], digits=3), 
                                   ")", 
                                   sep="")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  
  output[3,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[3]], digits=3)), #Fixed effect value (third level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[3]], digits=3), 
                                   ")", 
                                   sep="")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[3]], digits=4)))
  output[3,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[3]], digits=4)))
  
  output[4,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[4]], digits=3)), #Fixed effect value (fourth level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[4]], digits=3), 
                                   ")", 
                                   sep="")
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[4]], digits=4)))
  output[4,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[4]], digits=4)))
  
  output[5,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[5]], digits=3)), #Fixed effect value (fifth level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[5]], digits=3), 
                                   ")", 
                                   sep="")
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[5]], digits=4)))
  output[5,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[5]], digits=4)))
  output[1:5,"F value"] <- paste("-")
  
  output[1:5,"F value"] <- paste("-")
  
  ## Random effects
  output[6:7,"Coef (? se)"] <- paste("-")
  output[6,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[7,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  #output[8,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  
  output[6,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[7,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  #output[8,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[6:7,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 5) %>%
    pack_rows("Random effects", 6, 7)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}


save.output.gamm_single <- function(modelname, var1, title) {
  output <- as.data.frame(array(dim=c(5,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- "Session"
  output[4,"Variable"] <- "Site"
  output[5,"Variable"] <- "Subplot"
  output[c(1:2), 1] <- "Fixed"
  output[c(3:5), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[2]], digits=3), 
                                   ")", 
                                   sep="")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  output[1:2,"F value"] <- paste("-")
  
  
  ## Random effects
  output[3:5,"Coef (? se)"] <- paste("-")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[3,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[4,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  output[5,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[3:5,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 2) %>%
    pack_rows("Random effects", 3, 5)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}

save.output.gamm_single_ROOST <- function(modelname, var1, title) {
  output <- as.data.frame(array(dim=c(4,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- "Session"
  output[4,"Variable"] <- "Site"
  #output[5,"Variable"] <- "Subplot"
  output[c(1:2), 1] <- "Fixed"
  output[c(3:4), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)), #Fixed effect value (second level)
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[2]], digits=3), 
                                   ")", 
                                   sep="")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  output[1:2,"F value"] <- paste("-")
  
  
  ## Random effects
  output[3:4,"Coef (? se)"] <- paste("-")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  #output[5,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[3,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[4,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  #output[5,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[3:4,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 2) %>%
    pack_rows("Random effects", 3, 4)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}

save.output.gamm_null <- function(modelname, title) {
  output <- as.data.frame(array(dim=c(4,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- "Session"
  output[3,"Variable"] <- "Site"
  output[4,"Variable"] <- "Subplot"
  output[c(1), 1] <- "Fixed"
  output[c(2:4), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  ## Random effects
  output[2:4,"Coef (? se)"] <- paste("-")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[2,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[3,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  output[4,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[2:4,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 1) %>%
    pack_rows("Random effects", 2, 4)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}

save.output.gamm_null_ROOST <- function(modelname, title) {
  output <- as.data.frame(array(dim=c(3,6)))
  colnames(output) <- c("Effect type","Variable", "Coef (? se)", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- "Session"
  output[3,"Variable"] <- "Site"
  #output[4,"Variable"] <- "Subplot"
  output[c(1), 1] <- "Fixed"
  output[c(2:3), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef (? se)"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)), #Fixed effect intercept
                                   " (",
                                   round(summary(modelname$gam)[["se"]][[1]], digits=3), 
                                   ")", 
                                   sep="")
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  ## Random effects
  output[2:3,"Coef (? se)"] <- paste("-")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  #output[4,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[2,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[3,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  #output[4,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[2:3,"t value"] <- paste("-")
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 1) %>%
    pack_rows("Random effects", 2, 3)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  return(table_output)
}


############################################################################################################
################################ Save as df, not kable #####################################################
############################################################################################################

save.output.gamm_interaction_ROOST_df <- function(modelname, var1, var2, interact, title) {
  output <- as.data.frame(array(dim=c(6,7)))
  colnames(output) <- c("Effect","Variable", "Coef", "se", "t", "F", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- var2
  output[4,"Variable"] <- interact
  output[5,"Variable"] <- "z Session"
  output[6,"Variable"] <- "z Site"
  #output[7,"Variable"] <- "z Subplot"
  output[c(1:4), 1] <- "Fixed"
  output[c(5:6), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3))) #Fixed effect intercept
  output[1,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[1]], digits=3)))
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)))
  output[2,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[2]], digits=3)))
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  
  output[3,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[3]], digits=3)))
  output[3,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[3]], digits=3)))
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[3]], digits=4)))
  output[3,"t"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[3]], digits=4)))
  
  output[4,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[4]], digits=3)))
  output[4,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[4]], digits=3)))
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[4]], digits=4)))
  output[4,"t"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[4]], digits=4)))
  output[1:4,"F"] <- paste("-")
  
  ## Random effects
  output[5:6,"Coef (? se)"] <- paste("-")
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[6,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  #output[7,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[5,"F"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[6,"F"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  #output[7,"F"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[5:6,"t"] <- paste("-")
  
  output$Coef <- as.numeric(output$Coef)
  output$se <- as.numeric(output$se)
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 4) %>%
    pack_rows("Random effects", 5, 6)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  #return(table_output)
  return(output)
}


save.output.gamm_interactionfull_df <- function(modelname, var1, var2, interact, var3, title) {
  output <- as.data.frame(array(dim=c(8,7)))
  colnames(output) <- c("Effect","Variable", "Coef", "se", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- var2
  output[4,"Variable"] <- var3
  output[5,"Variable"] <- interact
  output[6,"Variable"] <- "z Session"
  output[7,"Variable"] <- "z Site"
  output[8,"Variable"] <- "z Subplot"
  output[c(1:5), 1] <- "Fixed"
  output[c(6:8), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)))
  output[1,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[1]], digits=3)))
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)))
  output[2,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[2]], digits=3)))                          
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  
  output[3,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[3]], digits=3)))
  output[3,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[3]], digits=3)))      
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[3]], digits=4)))
  output[3,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[3]], digits=4)))
  
  output[4,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[4]], digits=3)))
  output[4,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[4]], digits=3))) 
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[4]], digits=4)))
  output[4,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[4]], digits=4)))
  
  output[5,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[5]], digits=3)))
  output[5,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[5]], digits=3))) 
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[5]], digits=4)))
  output[5,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[5]], digits=4)))
  output[1:5,"F value"] <- paste("-")
  
  output[1:5,"F value"] <- paste("-")
  
  ## Random effects
  output[6:8,"Coef (? se)"] <- paste("-")
  output[6,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[7,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[8,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  
  output[6,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[7,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  output[8,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[6:8,"t value"] <- paste("-")
  
  output$Coef <- as.numeric(output$Coef)
  output$se <- as.numeric(output$se)
  
  return(output)
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 5) %>%
    pack_rows("Random effects", 6, 8)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  #return(table_output)
}


save.output.gamm_null_df <- function(modelname, title) {
  output <- as.data.frame(array(dim=c(4,7)))
  colnames(output) <- c("Effect","Variable", "Coef", "se", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- "z Session"
  output[3,"Variable"] <- "z Site"
  output[4,"Variable"] <- "z Subplot"
  output[c(1), 1] <- "Fixed"
  output[c(2:4), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)))
  output[1,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[1]], digits=3)))
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  ## Random effects
  output[2:4,"Coef (? se)"] <- paste("-")
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[2,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[3,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  output[4,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[2:4,"t value"] <- paste("-")
  
  output$Coef <- as.numeric(output$Coef)
  output$se <- as.numeric(output$se)
  
  return(output)
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 1) %>%
    pack_rows("Random effects", 2, 4)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  #return(table_output)
}


save.output.gamm_single_df <- function(modelname, var1, title) {
  output <- as.data.frame(array(dim=c(5,7)))
  colnames(output) <- c("Effect","Variable", "Coef", "se", "t value", "F value", "p")
  output[1,"Variable"] <- "Intercept"
  output[2,"Variable"] <- var1
  output[3,"Variable"] <- "z Session"
  output[4,"Variable"] <- "z Site"
  output[5,"Variable"] <- "z Subplot"
  output[c(1:2), 1] <- "Fixed"
  output[c(3:5), 1] <- "Random"
  
  ## Fixed effects
  output[1,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[1]], digits=3)))
  output[1,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[1]], digits=3)))
  output[1,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[1]], digits=4)))
  output[1,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[1]], digits=4)))
  
  output[2,"Coef"] <- paste(c(round(summary(modelname$gam)[["p.coeff"]][[2]], digits=3)))
  output[2,"se"] <- paste(c(round(summary(modelname$gam)[["se"]][[2]], digits=3)))
  output[2,"p"] <- paste(c(round(summary(modelname$gam)[["p.pv"]][[2]], digits=4)))
  output[2,"t value"] <- paste(c(round(summary(modelname$gam)[["p.t"]][[2]], digits=4)))
  output[1:2,"F value"] <- paste("-")
  
  
  ## Random effects
  output[3:5,"Coef"] <- paste("-")
  output[3,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[1]], digits=4)))
  output[4,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[2]], digits=4)))
  output[5,"p"] <- paste(c(round(summary(modelname$gam)[["s.pv"]][[3]], digits=4)))
  output[3,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][1,3], digits=4)))
  output[4,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][2,3], digits=4)))
  output[5,"F value"] <- paste(c(round(summary(modelname$gam)[["s.table"]][3,3], digits=4)))
  output[3:5,"t value"] <- paste("-")
  
  output$Coef <- as.numeric(output$Coef)
  output$se <- as.numeric(output$se)
  
  return(output)
  
  ## Create table
  table_output <- kable(output[,2:6], caption = title,  escape = F) %>%
    kable_styling("striped", full_width = F) %>%
    pack_rows("Fixed effects", 1, 2) %>%
    pack_rows("Random effects", 3, 5)
  #column_spec(c(3:11), width = "4cm") %>% #note - collumn_spec slow for long tables
  #scroll_box(height = "500px") %>%
  #save_kable(table_output,paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
  #return(table_output)
}


