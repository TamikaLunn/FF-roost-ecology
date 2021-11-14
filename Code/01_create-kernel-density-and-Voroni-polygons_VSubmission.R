## Title: R script for creating kernel density estimates and Voronoi Polygons for analyzing roosting structure of flying-fox roosts in SE QLD and NE NSW
## Manuscript: Conventional wisdom on roosting behavior of Australian flying-foxes—A critical review, and evaluation using new data <https://doi.org/10.1002/ece3.8079>
## Author: Tamika Lunn, Griffith University
## Version: VSubmission, created 13 November 2021

## V1-21 - manuscript preparation
## VSubmission - code copied from 'Bat spatial analysis_V21.R'

rm(list=ls())

##############################################################################################
##------------------------------------ Overview of data: -----------------------------------##
##############################################################################################

## See data_README

##############################################################################################
##----------------------------------------Set functions-------------------------------------##
##############################################################################################

## Load packages
library(spatstat) #spatial package
library(dplyr) #subsetting and managing data
library(RColorBrewer) #create custom colour ramps
library(ggplot2)
library(tidyverse)
library(dismo)

## Create custom functions
Liqr<-function(x) { 
  return(round(quantile(x,0.25,na.rm=TRUE),4)) 
}
Uiqr<-function(x) { 
  return(round(quantile(x,0.75,na.rm=TRUE),4)) 
}

Numextract <- function(string){
  unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
}

merge_lists <- function(list_tiles, split_tiles, output) { ##To merge lists based on position i.e. merge listA[[1]] & listB[[1]] etc.
  for(k in 1:length(list_tiles)) { ##note - only works because lists should be of the same length
    value_dirichlet <- list_tiles[[k]]
    tree.accession <- split_tiles[[k]][["marks"]]
    n <- as.data.frame(cbind(value_dirichlet, tree.accession))
    output[k,] <- n
  }
  return(output)
}

read.list <- function(series) {
  data <- list()
  for (i in 1:length(series)) { #read in all data and attach values
    data[[i]] <- read.csv(series[[i]], row.names=1)
    data[[i]]$site.code <- noquote(str_sub(series[[i]], 14, 17))
    data[[i]]$subplot <- noquote(regmatches(series[[i]], gregexpr("[[:digit:]]+", series[[i]]))[[1]]) 
  }
  return(data)
}

################################################################################################
##-----------------------------------------Start code-----------------------------------------##
################################################################################################

##--------------------------------------------------------------------------------------------##
##-------------------------Calculate and save kernel density estimates------------------------##
##--------------------------------------------------------------------------------------------##

data <- read.csv("Data/Raw/spatial-bat-structure-data.csv")

## Sum index weights to give indication of combined per tree abundance
data <- data %>%
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed
  mutate(all.index.weight = BFF.index.weight+GHFF.index.weight+LRFF.index.weight) #sum WEIGHTS, not index values, to give indication of overall bat load

plot.win <- owin(c(0,20), c(0,20)) # specify the plot window as 20x20m (the size of subplots)

## Create dataframe to store density values in
output <- data.frame(matrix(ncol = 14, nrow = length(unique(data$rep))))
x <- c("site.accession", "subplot", "GHFF_max_v", "GHFF_min_v", "GHFF_mean_v", "BFF_max_v", "BFF_min_v", "BFF_mean_v", "LRFF_max_v", "LRFF_min_v", "LRFF_mean_v","all_max_v", "all_min_v", "all_mean_v")
colnames(output) <- x

## Apply loop to store density values per site x time x plot combination
## Calculate summaries including zero value pixels:
counter = 1
for(i in unique(data$site.accession)){
  df = subset(data,data$site.accession == i)
  for (j in unique(df$subplot)){
    bat <- ppp(x = df$x[which(df$subplot == j)], y = df$y[which(df$subplot == j)], window = plot.win)
    GHFF_density<-density(bat, weights = df$GHFF.index.weight[which(df$subplot == j)], adjust=0.6)
    BFF_density<-density(bat, weights = df$BFF.index.weight[which(df$subplot == j)], adjust=0.6)
    LRFF_density<-density(bat, weights = df$LRFF.index.weight[which(df$subplot == j)], adjust=0.6)
    all_density<-density(bat, weights = df$all.index.weight[which(df$subplot == j)], adjust=0.6)
    output$GHFF_max_v[counter]<-max(GHFF_density$v) 
    output$BFF_max_v[counter]<-max(BFF_density$v) 
    output$LRFF_max_v[counter]<-max(LRFF_density$v) 
    output$all_max_v[counter]<-max(all_density$v) 
    output$GHFF_min_v[counter]<-min(GHFF_density$v) 
    output$BFF_min_v[counter]<-min(BFF_density$v) 
    output$LRFF_min_v[counter]<-min(LRFF_density$v) 
    output$all_min_v[counter]<-min(all_density$v) 
    output$GHFF_mean_v[counter]<-mean(GHFF_density$v) 
    output$BFF_mean_v[counter]<-mean(BFF_density$v) 
    output$LRFF_mean_v[counter]<-mean(LRFF_density$v) 
    output$all_mean_v[counter]<-mean(all_density$v)
    output$GHFF_median_v[counter]<-median(GHFF_density$v) 
    output$BFF_median_v[counter]<-median(BFF_density$v) 
    output$LRFF_median_v[counter]<-median(LRFF_density$v) 
    output$all_median_v[counter]<-median(all_density$v) 
    output$GHFF_Liqr_v[counter]<-Liqr(GHFF_density$v) 
    output$BFF_Liqr_v[counter]<-Liqr(BFF_density$v) 
    output$LRFF_Liqr_v[counter]<-Liqr(LRFF_density$v) 
    output$all_Liqr_v[counter]<-Liqr(all_density$v) 
    output$GHFF_Uiqr_v[counter]<-Uiqr(GHFF_density$v) 
    output$BFF_Uiqr_v[counter]<-Uiqr(BFF_density$v) 
    output$LRFF_Uiqr_v[counter]<-Uiqr(LRFF_density$v) 
    output$all_Uiqr_v[counter]<-Uiqr(all_density$v) 
    output$GHFF_sd_v[counter]<-sd(GHFF_density$v) 
    output$BFF_sd_v[counter]<-sd(BFF_density$v) 
    output$LRFF_sd_v[counter]<-sd(LRFF_density$v) 
    output$all_sd_v[counter]<-sd(all_density$v)
    output$site.accession[counter]<-i
    output$subplot[counter]<-j
    counter <- counter +1
  }
}
tail(output)
output$site.accession<-as.factor(output$site.accession)

## Save output
write.csv(output,file="pixel-density-data.csv") #For use in later analyses
#saveRDS(output, "pixel-density-data.rds") #Save data with all data structure changes 

## Calculate summaries not including zero value pixels:
counter = 1
for(i in unique(data$site.accession)){
  df = subset(data,data$site.accession == i)
  for (j in unique(df$subplot)){
    bat <- ppp(x = df$x[which(df$subplot == j)], y = df$y[which(df$subplot == j)], window = plot.win)
    GHFF_density<-density(bat, weights = df$GHFF.index.weight[which(df$subplot == j)], adjust=0.6)
    BFF_density<-density(bat, weights = df$BFF.index.weight[which(df$subplot == j)], adjust=0.6)
    LRFF_density<-density(bat, weights = df$LRFF.index.weight[which(df$subplot == j)], adjust=0.6)
    all_density<-density(bat, weights = df$all.index.weight[which(df$subplot == j)], adjust=0.6)
    output$GHFF_max_v[counter]<-max(GHFF_density$v) 
    output$BFF_max_v[counter]<-max(BFF_density$v) 
    output$LRFF_max_v[counter]<-max(LRFF_density$v) 
    output$all_max_v[counter]<-max(all_density$v) 
    output$GHFF_min_v[counter]<-min(GHFF_density$v[GHFF_density$v>0.1])
    output$BFF_min_v[counter]<-min(BFF_density$v[BFF_density$v>0.1]) 
    output$LRFF_min_v[counter]<-min(LRFF_density$v[LRFF_density$v>0.1]) 
    output$all_min_v[counter]<-min(all_density$v[all_density$v>0.1]) 
    output$GHFF_mean_v[counter]<-mean(GHFF_density$v[GHFF_density$v>0.1]) 
    output$BFF_mean_v[counter]<-mean(BFF_density$v[BFF_density$v>0.1]) 
    output$LRFF_mean_v[counter]<-mean(LRFF_density$v[LRFF_density$v>0.1]) 
    output$all_mean_v[counter]<-mean(all_density$v[all_density$v>0.1]) 
    output$GHFF_median_v[counter]<-median(GHFF_density$v[GHFF_density$v>0.1]) 
    output$BFF_median_v[counter]<-median(BFF_density$v[BFF_density$v>0.1]) 
    output$LRFF_median_v[counter]<-median(LRFF_density$v[LRFF_density$v>0.1]) 
    output$all_median_v[counter]<-median(all_density$v[all_density$v>0.1]) 
    output$GHFF_Liqr_v[counter]<-Liqr(GHFF_density$v[GHFF_density$v>0.1]) 
    output$BFF_Liqr_v[counter]<-Liqr(BFF_density$v[BFF_density$v>0.1]) 
    output$LRFF_Liqr_v[counter]<-Liqr(LRFF_density$v[LRFF_density$v>0.1]) 
    output$all_Liqr_v[counter]<-Liqr(all_density$v[all_density$v>0.1]) 
    output$GHFF_Uiqr_v[counter]<-Uiqr(GHFF_density$v[GHFF_density$v>0.1]) 
    output$BFF_Uiqr_v[counter]<-Uiqr(BFF_density$v[BFF_density$v>0.1]) 
    output$LRFF_Uiqr_v[counter]<-Uiqr(LRFF_density$v[LRFF_density$v>0.1]) 
    output$all_Uiqr_v[counter]<-Uiqr(all_density$v[all_density$v>0.1]) 
    output$GHFF_sd_v[counter]<-sd(GHFF_density$v[GHFF_density$v>0.1]) 
    output$BFF_sd_v[counter]<-sd(BFF_density$v[BFF_density$v>0.1]) 
    output$LRFF_sd_v[counter]<-sd(LRFF_density$v[LRFF_density$v>0.1]) 
    output$all_sd_v[counter]<-sd(all_density$v[all_density$v>0.1]) 
    output$site.accession[counter]<-i
    output$subplot[counter]<-j
    counter <- counter +1
  }
}
tail(output) #NaN and Inf values generated when max=0 (no bats)

output <- output %>%
  mutate(GHFF_min_v = replace(GHFF_min_v, GHFF_min_v=="Inf", 0)) %>%
  mutate(GHFF_mean_v = replace(GHFF_mean_v, GHFF_mean_v=="NaN", 0)) %>%
  mutate(BFF_min_v = replace(BFF_min_v, BFF_min_v=="Inf", 0)) %>%
  mutate(BFF_mean_v = replace(BFF_mean_v, BFF_mean_v=="NaN", 0)) %>%
  mutate(LRFF_min_v = replace(LRFF_min_v, LRFF_min_v=="Inf", 0)) %>%
  mutate(LRFF_mean_v = replace(LRFF_mean_v, LRFF_mean_v=="NaN", 0)) %>%
  mutate(all_min_v = replace(all_min_v, all_min_v=="Inf", 0)) %>%
  mutate(all_mean_v = replace(all_mean_v, all_mean_v=="NaN", 0)) %>%
  mutate(GHFF_median_v = replace_na(GHFF_median_v, 0)) %>%
  mutate(BFF_median_v = replace_na(BFF_median_v, 0)) %>%
  mutate(LRFF_median_v = replace_na(LRFF_median_v, 0)) %>%
  mutate(all_median_v = replace_na(all_median_v, 0)) %>%
  mutate(GHFF_Liqr_v = replace_na(GHFF_Liqr_v, 0)) %>%
  mutate(BFF_Liqr_v = replace_na(BFF_Liqr_v, 0)) %>%
  mutate(LRFF_Liqr_v  = replace_na(LRFF_Liqr_v, 0)) %>%
  mutate(all_Liqr_v  = replace_na(all_Liqr_v, 0)) %>%
  mutate(GHFF_Uiqr_v = replace_na(GHFF_Uiqr_v, 0)) %>%
  mutate(BFF_Uiqr_v = replace_na(BFF_Uiqr_v, 0)) %>%
  mutate(LRFF_Uiqr_v = replace_na(LRFF_Uiqr_v, 0)) %>%
  mutate(all_Uiqr_v = replace_na(all_Uiqr_v, 0)) %>%
  mutate(GHFF_sd_v = replace_na(GHFF_sd_v, 0)) %>%
  mutate(BFF_sd_v = replace_na(BFF_sd_v, 0)) %>%
  mutate(LRFF_sd_v = replace_na(LRFF_sd_v, 0)) %>%
  mutate(all_sd_v = replace_na(all_sd_v, 0))

output$site.accession<-as.factor(output$site.accession)

## Save output
write.csv(output,file="Data/Raw/pixel-density-data-nonzero.csv") #For use in later analyses
#saveRDS(output, "Data/Raw/pixel-density-data-nonzero.rds") #Save data with all data structure changes 


##--------------------------------------------------------------------------------------------##
##------------------------------Create maps of density (intensity)----------------------------##
##--------------------------------------------------------------------------------------------##

data <- read.csv("Data/Raw/spatial-bat-structure-data.csv")
data <- data %>%
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed
  mutate(all.index.weight = BFF.index.weight+GHFF.index.weight+LRFF.index.weight) #sum WEIGHTS, not index values, to give indication of overall bat load
output<-read.csv("pixel-density-data.csv")
tail(output) #make sure all accessions are there
tail(data) #make sure all accessions are there

#########################################################################################
### Between sites fixed scale plots (fixed scale across all site x time combinations) ###
#########################################################################################

## Choose directory to save plots ##
setwd("Output/Figures/Cluster plots/Scale fixed between sites")

## Save maximum values to calculate scale for plots
GHFF_vmax<-max(output$GHFF_max_v) #save maximum v value from all sites = max value for scale
BFF_vmax<-max(output$BFF_max_v) #save maximum v value from all sites = max value for scale
LRFF_vmax<-max(output$LRFF_max_v) 
all_vmax<-max(output$all_max_v) 

vmax<-max(c(GHFF_vmax, BFF_vmax, LRFF_vmax, all_vmax)) #save maxiumum v value from all sites = max value for scale
vmax #check value is sensible

## Set plotting window etc outside of loop/function
plot.win <- owin(c(0,20), c(0,20)) # determining the plot window as 20x20m
spat.bat.frame <- NULL # create a frame to store values in

cols1<-colorRampPalette(brewer.pal(8,"Blues"))(100) #define colour ramp for intensity plots
cols1[1]<- "#FFFFFF00" #sets the lowest value as not white (white = #ffffff)

cols2<-colorRampPalette(brewer.pal(8,"OrRd"))(100) #define colour ramp for intensity plots
cols2[1]<- "#FFFFFF00"

cols3<-colorRampPalette(brewer.pal(4,"Greys"))(100) #define colour ramp for intensity plots
cols3[1] <- "#FFFFFF00"

## Loop for spatial plots 
for(i in unique(data$site.accession)){
  df = subset(data,data$site.accession == i)
  jpeg(filename = paste(i, "PLOT.jpg", sep = ""), res = 300, units = "cm", width = 7, height = 35)
  par(mfrow=c(10,1),oma=c(0,0,0,0),mar=c(0,5,1,0),lwd=1)
  #outer margin: oma=c(b,l,t,r), inner margin: mar=c(b,l,t,r). Separated by lines of space 
  for (j in unique(df$subplot)){
    bat <- ppp(x = df$x[which(df$subplot == j)], y = df$y[which(df$subplot == j)], window = plot.win)
    tree <- ppp(x = df$x[which(df$subplot == j)], y = df$y[which(df$subplot == j)], window = plot.win)
    plot(density(bat, weights = df$GHFF.index.weight[which(df$subplot == j)], adjust=0.6), sigma = 1, col=cols1, main=paste(i, 'plot', j), zlim = range(zlim = range(0,vmax))) 
    plot(density(bat, weights = df$BFF.index.weight[which(df$subplot == j)], adjust=0.6), sigma = 1, col=cols2, main=paste(i, 'plot', j), zlim = range(zlim = range(0,vmax)), add=TRUE)
    plot(density(bat, weights = df$LRFF.index.weight[which(df$subplot == j)], adjust=0.6), sigma = 1, col=cols3, main=paste(i, 'plot', j), zlim = range(zlim = range(0,vmax)), add=TRUE)
    #plot(density(bat, weights = df$all.index.weight[which(df$subplot == j)], adjust=0.6), sigma = 1, col=cols3, main=paste(i, 'plot', j), zlim = range(zlim = range(0,vmax))) #can't have all bats, because will override colours of individual species
    plot(tree, add=TRUE) #add trees to density plot
  }
  dev.off()
}

##########################################################################################
#### Within sites fixed scale plots (fixed scale across all site x time combinations) ####
##########################################################################################

## Set plotting window etc outside of loop/function
plot.win <- owin(c(0,20), c(0,20)) # determining the plot window as 20x20m
spat.bat.frame <- NULL # create a frame to store values in

cols1<-colorRampPalette(brewer.pal(8,"Blues"))(100) #define colour ramp for intensity plots
cols1[1]<- "#FFFFFF00" #sets the lowest value as not white (white = #ffffff)

cols2<-colorRampPalette(brewer.pal(8,"OrRd"))(100) #define colour ramp for intensity plots
cols2[1]<- "#FFFFFF00"

cols3<-colorRampPalette(brewer.pal(4,"Greys"))(100) #define colour ramp for intensity plots
cols3[1] <- "#FFFFFF00"

## Save maximum values to calculate scale for plots:
subset<-subset(output, regexpr("DAVO", output$site.accession) > 0) #pulls out rows from pixel output with site.accession beginning in specified site (e.g. "DAVO" will pull out "DAVO001", "DAVO002" etc). repeat below code with different subsets of data
sitesubset <- data[which(data$site.code=="DAVO"),] #subset actual data with specied site to run loop. repeat below code with different subsets of data
tail(subset) #Quick check only one site is in dataframe
GHFF_vmax<-max(subset$GHFF_max_v) #save maxiumum v value from all sites 
BFF_vmax<-max(subset$BFF_max_v) #save maxiumum v value from all sites 
LRFF_vmax<-max(subset$LRFF_max_v) #save maxiumum v value from all sites
all_vmax<-max(subset$all_max_v) 
vmax<-max(c(GHFF_vmax, BFF_vmax, LRFF_vmax, all_vmax))  #Save maxiumum v value from all sites from all species. This will be the value used to set the scale
vmax #check its an expected number
## Need to repeat above chunk manually for each roost site (DAVO, DBUR, DCAN, DCLU, DLIS, DRED, DSUN, DTOW), ahead of running loops below:

### Loop for spatial plots - 10 plots per image (repeat manually for each roost site) ###

## Choose directory to save plots ##
setwd("Output/Figures/Cluster plots/Scale fixed within sites")

for(i in unique(sitesubset$site.accession)){
  df = subset(sitesubset,sitesubset$site.accession == i)
  jpeg(filename = paste(i, "PLOT.jpg", sep = ""), res = 300, units = "cm", width = 7, height = 35)
  par(mfrow=c(10,1),oma=c(0,0,0,0),mar=c(0,5,1,0),lwd=1)
  #outer margin: oma=c(b,l,t,r), inner margin: mar=c(b,l,t,r). Separated by lines of space 
  for (j in unique(df$subplot)){
    bat <- ppp(x = df$x[which(df$subplot == j)], y = df$y[which(df$subplot == j)], window = plot.win)
    tree <- ppp(x = df$x[which(df$subplot == j)], y = df$y[which(df$subplot == j)], window = plot.win)
    plot(density(bat, weights = df$GHFF.index.weight[which(df$subplot == j)], adjust=0.6), sigma = 1, col=cols1, main=paste(i, 'plot', j), zlim = range(zlim = range(0,vmax)))
    plot(density(bat, weights = df$BFF.index.weight[which(df$subplot == j)], adjust=0.6), sigma = 1, col=cols2, main=paste(i, 'plot', j), zlim = range(zlim = range(0,vmax)), add=TRUE)
    plot(density(bat, weights = df$LRFF.index.weight[which(df$subplot == j)], adjust=0.6), sigma = 1, col=cols3, main=paste(i, 'plot', j), zlim = range(zlim = range(0,vmax)), add=TRUE)
    #plot(density(bat, weights = df$all.index.weight[which(df$subplot == j)], adjust=0.6), sigma = 1, col=cols4, main=paste(i, 'plot', j), zlim = range(zlim = range(0,vmax))) #can't have all bats, because will override colours of individual species
    plot(tree, add=TRUE) #add trees to density plot
  }
  dev.off()
}


### Loop for spatial plots - single image per plot (repeat manually for each roost site) ###
## Choose directory to save plots ##
setwd("Output/Figures/Cluster plots/Scale fixed within sites/Avondale")
#setwd("Output/Figures/Cluster plots/Scale fixed within sites/Burleigh")
#setwd("Output/Figures/Cluster plots/Scale fixed within sites/Canungra")
#setwd("Output/Figures/Cluster plots/Scale fixed within sites/Clunes")
#setwd("Output/Figures/Cluster plots/Scale fixed within sites/Lismore")
#setwd("Output/Figures/Cluster plots/Scale fixed within sites/Redcliffe")
#setwd("Output/Figures/Cluster plots/Scale fixed within sites/Sunnybank")
#setwd("Output/Figures/Cluster plots/Scale fixed within sites/Toowoomba")

for(i in unique(sitesubset$site.accession)){
  df = subset(sitesubset,sitesubset$site.accession == i)
  #par(mfrow=c(10,1),oma=c(0,0,0,0),mar=c(0,5,1,0),lwd=1)
  #outer margin: oma=c(b,l,t,r), inner margin: mar=c(b,l,t,r). Separated by lines of space 
  for (j in unique(df$subplot)){
    jpeg(filename = paste(i, ' plot ', j, ".jpg", sep = ""), res = 300, units = "cm", width = 10, height = 10)
    par(mar=c(0,0,0,0),oma=c(0,0,0,0),lwd=1)
    bat <- ppp(x = df$x[which(df$subplot == j)], y = df$y[which(df$subplot == j)], window = plot.win)
    tree <- ppp(x = df$x[which(df$subplot == j)], y = df$y[which(df$subplot == j)], window = plot.win)
    plot(density(bat, weights = df$GHFF.index.weight[which(df$subplot == j)], adjust=0.6), ribbon=F, main=NULL, sigma = 1, col=cols1, zlim = range(zlim = range(0,vmax)))
    plot(density(bat, weights = df$BFF.index.weight[which(df$subplot == j)], adjust=0.6), ribbon=F, main=NULL, sigma = 1, col=cols2, zlim = range(zlim = range(0,vmax)), add=TRUE)
    plot(density(bat, weights = df$LRFF.index.weight[which(df$subplot == j)], adjust=0.6), ribbon=F, main=NULL, sigma = 1, col=cols3, zlim = range(zlim = range(0,vmax)), add=TRUE)
    #plot(density(bat, weights = df$all.index.weight[which(df$subplot == j)], adjust=0.6, legend=FALSE), ribbon=F, main=NULL, col=cols4, sigma = 1, zlim = range(zlim = range(0,vmax))) #can't have all bats, because will override colours of individual species
    plot(tree, pch=19, size=1.5, add=TRUE) #add trees to density plot
    dev.off()
  }
  #Save a single plot per species, for the scale bar:
  jpeg(filename = paste(unique(sitesubset$site.code), "_scale_GHFF.jpg", sep = ""), res = 300, units = "cm", width = 10, height = 10)
  par(mar=c(1,1,1,1),oma=c(1,1,1,1),lwd=1)
  plot(density(bat, weights = df$GHFF.index.weight[which(df$subplot == j)], adjust=0.6), main=NULL, sigma = 1, col=cols1, zlim = range(zlim = range(0,vmax)))
  dev.off()
  
  jpeg(filename = paste(unique(sitesubset$site.code), "_scale_BFF.jpg", sep = ""), res = 300, units = "cm", width = 10, height = 10)
  par(mar=c(1,1,1,1),oma=c(1,1,1,1),lwd=1)
  plot(density(bat, weights = df$BFF.index.weight[which(df$subplot == j)], adjust=0.6), main=NULL, sigma = 1, col=cols2, zlim = range(zlim = range(0,vmax)))
  dev.off()
  
  jpeg(filename = paste(unique(sitesubset$site.code), "_scale_LRFF.jpg", sep = ""), res = 300, units = "cm", width = 10, height = 10)
  par(mar=c(1,1,1,1),oma=c(1,1,1,1),lwd=1)
  plot(density(bat, weights = df$LRFF.index.weight[which(df$subplot == j)], adjust=0.6), main=NULL, sigma = 1, col=cols3, zlim = range(zlim = range(0,vmax)))
  dev.off()
  }


######################################################################
#################### Histogram of pixel values #######################
######################################################################
plot.win <- owin(c(0,20), c(0,20)) # determining the plot window as 20x20m
spat.bat.frame <- NULL # create a frame to store values in

### Between sites fixed scale plots (fixed scale across all site x time combinations) ###
## Run separately for each species:
setwd("Output/Figures/Histograms of pixels/BFF")
#setwd("Output/Figures/Histograms of pixels/GHFF")
#setwd("Output/Figures/Histograms of pixels/LRFF")

# Loop for spatial plots (repeat for species, and change working directory)
for(i in unique(data$site.accession)){
  df = subset(data,data$site.accession == i)
  jpeg(filename = paste(i, "PLOT.jpg", sep = ""), res = 300, units = "cm", width = 7, height = 35)
  par(mfrow=c(10,1),oma=c(0,0,0,0),mar=c(2,5,1,0),lwd=1)
  #outer margin: oma=c(b,l,t,r), inner margin: mar=c(b,l,t,r). Separated by lines of space 
  for (j in unique(df$subplot)){
    bat <- ppp(x = df$x[which(df$subplot == j)], y = df$y[which(df$subplot == j)], window = plot.win)
    #pixels <- density(bat, weights = df$BFF.index.weight[which(df$subplot == j)], adjust=0.6) 
    pixels <- density(bat, weights = df$GHFF.index.weight[which(df$subplot == j)], adjust=0.6) 
    #pixels <- density(bat, weights = df$LRFF.index.weight[which(df$subplot == j)], adjust=0.6) 
    
    store <- hist(pixels$v, breaks=70, plot=FALSE)
    ymax <- max(store$density)
    xLiqr<-Liqr(pixels$v)
    xUiqr <- Uiqr(pixels$v)
    hist(pixels$v, breaks=70, freq=FALSE, xlab="pixel density value", ylab="density of pixel values", main=paste(i, 'plot', j), border = "black", col='grey', mgp=c(1.75,0.5,0)) #gives freq=FALSE gives probability density
    polygon(c(xLiqr,xUiqr, xUiqr,xLiqr),c(0,0,ymax,ymax), col=rgb(0.1,0.1,0.1,0.2), border = NA)
    curve(dnorm(x, mean=mean(pixels$v), sd=sd(pixels$v)), add=TRUE, lwd=1)
    abline(v = mean(pixels$v), col = "royalblue", lwd = 2)
    abline(v = median(pixels$v), col = "red", lwd = 2)
      }
  dev.off()
}

##--------------------------------------------------------------------------------------------##
##-----------------------------------Create Voronoi Polygons----------------------------------##
##--------------------------------------------------------------------------------------------##
data <- read.csv("Data/Raw/tree survey data.csv") ## Using initial tree survey data. This assumes that canopy doesn't change in response to the few cases where trees were cleared
data <- data %>%
  filter(!is.na(tree.accession)) %>% #removes missing trees
  filter(crown == "C" | crown == "D" | crown == "CI") ## remove trees that aren't apart of the canopy
unique(data$crown) 

## Set plotting window etc outside of loop/function
plot.win <- owin(c(0,20), c(0,20)) # determining the plot window as 20x20m
spat.bat.frame <- NULL # create a frame to store values in

## Loop for plots & save output per plot 
for(i in unique(data$site.code)) {
  df = subset(data,data$site.code == i)
  for (j in unique(df$subplot)){
    bat <- ppp(x = df$x[which(df$subplot == j)], y = df$y[which(df$subplot == j)], window = plot.win)
    tree <- ppp(x = df$x[which(df$subplot == j)], y = df$y[which(df$subplot == j)], window = plot.win, marks = df$tree.accession[which(df$subplot == j)])
    tile <- dirichlet(tree)
    split_tiles <- split(tree, tile) #Separates the point pattern (first arg) into sub-patterns delineated by the tiles of the tessellation (second argument). In output: x and y are xy locations from spp. marks gives tree.accession
    list_tiles <- lapply(split_tiles, dirichletAreas) ##Use lapply to perform operations on each element of the list
    
    ## Merge outputs to retain tree.accession and tessellation area:
    setwd("Data/Raw/Voronoi Polygons/Area of polygons per plot")
    output <- data.frame(matrix(ncol = 2, nrow = length(list_tiles)))
    x <- c("value_dirichlet", "tree.accession")
    colnames(output) <- x
    output <- merge_lists(list_tiles, split_tiles, output) 
    #sum(as.numeric(output$value_dirichlet)) #this should add to 400, as the area of the plot is 400 meters squared
    write.csv(output, paste("voronoi area_",i,"_subplot ",j,".csv",sep=""))
    
    ## Visualise plots
    setwd("Data/Raw/Voronoi Polygons/Visuals per plot")
    jpeg(paste(i,"_subplot ",j,".png",sep=""), width=600, height=600)
    plot(tree, use.marks=FALSE, main=paste(i)) #add trees to density plot
    plot(tile, add=TRUE)
    dev.off()
  }
}

## Read in and merge files
setwd("Data/Raw/Voronoi Polygons/Area of polygons per plot")
series <- list.files() #need for every read-in, because file names are different
output <- do.call(rbind.data.frame,read.list(series)) %>% #read from list and convert to dataframe
  mutate(tree.accession = as.factor(tree.accession))  %>%
  mutate(site.code = as.factor(site.code))  %>%
  mutate(subplot = as.factor(subplot))
summary(output)

data <- read.csv("Data/Raw/tree survey data.csv") ## Using initial tree survey data. This assumes that canopy doesn't change in response to the few cases where trees were cleared
data <- data %>% ## re-read data and retain non-canopy trees for merge
  filter(!is.na(tree.accession)) %>% #removes missing trees
  mutate(tree.accession = as.factor(tree.accession))  %>%
  mutate(site.code = as.factor(site.code))  %>%
  mutate(subplot = as.factor(subplot))
summary(data) #make sure all accessions are there

tree.tessellation <- dplyr::full_join(output,
                                      data[,c("tree.accession", "site.code", "subplot","crown","x","y")], #full join to retain non canopy trees that were filtered from output
                 by = c("tree.accession", "site.code", "subplot"))

## Save combined tessellation
write.csv(tree.tessellation, "Data/Raw/tree.tessellation.csv")

##--------------------------------------------------------------------------------------------##