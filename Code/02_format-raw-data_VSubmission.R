## Title: R script for analysing roosting structure of flying-fox roosts in SE QLD and NE NSW
## Manuscript: Conventional wisdom on roosting behavior of Australian flying-foxes—A critical review, and evaluation using new data <https://doi.org/10.1002/ece3.8079>
## Author: Tamika Lunn, Griffith University
## Version: VSubmission, created 13 November 2021

## V1-7 - manuscript preparation
## VSubmission - code copied from 'Bat quantitative analysis_FF#1_V7-revision.R'

rm(list=ls())

##############################################################################################
##------------------------------------ Overview of data: -----------------------------------##
##############################################################################################

## See data_README

##############################################################################################
##---------------------------------Load data & set functions--------------------------------##
##############################################################################################

source ("Code/00_functions_VSubmission.R") ## Read from relative path
library(binom)

## Tree survey data:
treesurvey <- read.csv("Data/Raw/tree survey data.csv")

## Extra centroid data
centroids <- read.csv("Data/Raw/Centroids.csv") %>%
  mutate(roost.centroid.N = ifelse(site.accession=="DLIS010",NA,roost.centroid.N)) %>% ## Remove May 2019 Lismore area point (DLIS010)
  mutate(roost.centroid.E = ifelse(site.accession=="DLIS010",NA,roost.centroid.E)) 
centroids$session <- as.factor(centroids$session)

## Bat structure data:
treebat <- read.csv("Data/Raw/spatial-bat-structure-data.csv", row.names=1)
treebat$subplot <- as.factor (treebat$subplot)
treebat$session <- as.factor (treebat$session)
head(treebat)[,c("site.code", "tree.accession", "site.accession", "BFF.index", "GHFF.index", "LRFF.index", "subplot", "crowngroup")]

## Pixel density data:
##with zero values:
kernelbat <- read.csv("Data/Raw/pixel-density-data.csv", row.names=1)
kernelbat$subplot <- as.factor (kernelbat$subplot)
kernelbat <- arrange(kernelbat, site.accession, subplot) #arrange makes customising kables easier
##without zero values:
kernelbatnz <- read.csv("Data/Raw/pixel-density-data-nonzero.csv", row.names=1)
kernelbatnz$subplot <- as.factor (kernelbatnz$subplot)
kernelbatnz <- arrange(kernelbatnz, site.accession, subplot) #arrange makes customising kables easier

## Roost level data:
roostbat<-read.csv("Data/Raw/ALL-roost-use-data.csv", row.names=1) %>% filter(!is.na(site.code)) %>%
  mutate(roost.area = ifelse(site.accession=="DLIS010",NA,roost.area)) %>%
  mutate(roost.perimeter = ifelse(site.accession=="DLIS010",NA,roost.perimeter)) ## Remove May 2019 Lismore area point (DLIS010)
roostbat$session <- as.factor (roostbat$session)
head(roostbat)[,c("site.code", "session", "site.accession", "pop.estimate.BFF", "pop.estimate.GHFF", "pop.estimate.LRFF", "pop.estimate.total")]

temp.tree <- treebat %>%
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed. This is important to do so that proportion is calculated correctly
  plyr::ddply(c("site.code", "session", "site.accession"), summarise,
              total.trees = sum(!is.na(tree.accession)))
roostbat <- dplyr::left_join(roostbat, temp.tree, by=c("site.code", "session", "site.accession"))

## Additional subplot data:
plotdata <- read.csv("Data/Raw/plot-data.csv", row.names=1) %>%
  mutate(site.code = as.factor(site.code)) %>%
  mutate(subplot = as.factor(subplot)) %>%
  mutate(roost.type = as.factor(roost.type)) %>%
  mutate(plotID = as.factor(plotID))

## Tree tessellation data:
tree.tessellation <- read.csv("Data/Raw/Voronoi Polygons/tree.tessellation.csv", row.names = 1) %>%
  mutate(tree.accession = as.factor(tree.accession))  %>%
  mutate(site.code = as.factor(site.code))  %>%
  mutate(subplot = as.factor(subplot))  %>%
  mutate(crown = as.factor(crown))  %>%
  mutate(value_dirichlet = ifelse(crown=="I",Liqr(value_dirichlet), ## Set area of intermediate trees as the first quantile of the canopy trees (5.8 m2)
                                  ifelse(crown=="E",199, 
                                         ifelse(crown=="OG",199,
                                                value_dirichlet)))) %>%
  mutate(value_dirichlet = ifelse(value_dirichlet>199,199,value_dirichlet)) ## set maximum value, based on mean crown projection area for eucalypts in NSW, reported in Verma (2014) An allometric model for estimating DBH of isolated and clustered Eucalyptus trees from measurements of crown projection area

################################################################################################
##-----------------------------------------Start code-----------------------------------------##
################################################################################################

##--------------------------------------------------------------- Format data ---------------------------------------------------------------##
## There are two main data types for bat data: 'index-based' and 'count-based'. Index-based data were taken for all trees (and include index and weighted index values), whereas count was taken for a subset of trees (and include count, min and max height).
## 'pixel-based' data are kernel density estimates, calculated from the spatial point pattern of tree distributions within plots weighted by the 'index-based' data per tree. 
## There are two possible summary levels: plot-level and roost-level. Note - in the published paper, 'plot' is referred to as 'subplot', and 'site' as 'roost'.
## There are two possible calculations for roost-level summaries: with tree as the lowest replicate ('bytree'), and plot as the lowest replicate ('byplot')
## Summary data frames are returned in either a wide format (species in one column, then separate columns for each summary value) OR long format (species in one column, then one column for 'value' and one column for 'measure'). May also be an 'extra-wide' format which reflects the way the data was entered (separate columns for site*measure combinations)
## Wide format for when named columns are to be selected for use. Long format for when values need to be filtered from a single column

### Naming convention for functions:
###      Functions: 'data type' '.' 'LOWEST REPLICATE' '.' 'returned df format'
###      Output: 'data type' '.' 'LOWEST REPLICATE' '.' 'returned df format' '_' 'summary level'
###      data type = index | count
###      LOWEST REPLICATE = TREE | PLOT
###      returned df format = wide |Ewide | long
###      summary level = site | plot

#############################################
#### Plot-level (bat structure) overview #### 
#############################################

### Calculate data summaries:
byplot <- c("site.code", "session", "site.accession", "subplot", "rep") ## Group by subplot
bytree <- c("site.code", "session", "site.accession", "tree.accession", "id", "crowngroup", "subplot", "rep") ## Group by tree
list.count <- treebat  %>% #get list of accessions with tree count 
  filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>% #remove rare cases where trees missed in survey OR were removed
  ddply(c("site.accession", "site.code", "session", "subplot", "rep"), summarise,
        tree.count = sum(!is.na(tree.accession)))  
list.kernel <- treebat  %>% #get list of accessions with tree count
  filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>% #remove rare cases where trees missed in survey OR were removed
  ddply(c("site.accession", "site.code", "session", "subplot"), summarise,
        tree.count = sum(!is.na(tree.accession)))  

## Calculate summaries & transpose dataframes into different formats
index.TREE.Ewide_plot <- index.TREE.Ewide(treebat, byplot) ## Dataframe with separate column per measure*species
index.TREE.wide_plot <- index.TREE.wide(index.TREE.Ewide_plot, byplot, list.count) ## Dataframe with separate column per measure
index.TREE.long_plot <- index.TREE.long(index.TREE.Ewide_plot, byplot) 
heights.subset <- calculate.heightdiff(treebat, bytree)  %>%
  dplyr::left_join(
    tree.tessellation[,c("tree.accession", "site.code", "subplot","value_dirichlet")], 
    by = c("tree.accession", "site.code", "subplot")) %>% 
  mutate(vert.density = count/(height.diff*value_dirichlet)) %>% 
  mutate(vert.density = ifelse(vert.density==Inf,NA,vert.density))
count.TREE.wide_plot <- count.TREE.wide(heights.subset, byplot, list.count) 
count.TREE.long_plot <- count.TREE.long(count.TREE.wide_plot, byplot)
kernel.TREE.wide_plot<- kernel.TREE.wide(kernelbat, list.kernel)
kernel.TREE.long_plot<- kernel.TREE.long(kernelbat, list.kernel)
kernelnz.TREE.wide_plot<- kernel.TREE.wide(kernelbatnz, list.kernel)
kernelnz.TREE.long_plot<- kernel.TREE.long(kernelbatnz, list.kernel)

#############################################
#### Site-level (bat structure) overview #### 
#############################################

### Calculate data summaries. 
## Note - There are two calculations for site summaries: with tree as the lowest replicate, and subplot as the lowest replicate. When subplot is the lowest replicate, averages are calculated across occupied subplots only. Mean in this case is the mean of the subplot means. Standard deviation here is the standard deviation between subplots. 

##---------------------------------------- Subplot as lowest replicate ---------------------------------------- 
byplot <- c("site.code", "session", "site.accession", "subplot", "rep")
bysite <- c("site.code", "session", "site.accession")

## Create list of accessions and measures for merging:
list.count_site <- treebat  %>% #get list of accessions with tree count, doesn't matter which dataframe 
  filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>% #remove rare cases where trees missed in survey OR were removed
  ddply(c("site.accession", "site.code", "session"), summarise,
        tree.count = sum(!is.na(tree.accession)))  
list.index_site <- index.TREE.long_plot  %>%
  distinct(site.accession, site.code, session, measure) %>%
  mutate(mean.value = as.numeric(0)) %>%
  left_join(list.count_site, by=c("site.accession", "site.code", "session")) ##Left join to get tree count per site
list.count_site <-index.TREE.long_plot  %>%
  distinct(site.accession, site.code, session) %>%
  left_join(list.count_site, by=c("site.accession", "site.code", "session"))
list.pixel_site <- kernel.TREE.long_plot  %>%
  distinct(site.accession, site.code, session, measure) %>%
  mutate(mean.value = as.numeric(0)) %>%
  left_join(list.count_site, by=c("site.accession", "site.code", "session"))

##identify occupied subplots:
occ.all.output_plot <- occ.all(treebat, byplot)
occ.BFF.output_plot <- occ.BFF(treebat, byplot)
occ.GHFF.output_plot <- occ.GHFF(treebat, byplot)
occ.LRFF.output_plot <- occ.LRFF(treebat, byplot)

##calculate:
index.PLOT.long_BFF.site <- index.PLOT.long.sp(index.TREE.long_plot, bysite, "BFF", occ.BFF.output_plot, list.index_site) #by subplot, because you want to select occupied subplots 
index.PLOT.long_GHFF.site <- index.PLOT.long.sp(index.TREE.long_plot, bysite, "GHFF", occ.GHFF.output_plot, list.index_site)
index.PLOT.long_LRFF.site <- index.PLOT.long.sp(index.TREE.long_plot, bysite, "LRFF", occ.LRFF.output_plot, list.index_site)
index.PLOT.long_all.site <- index.PLOT.long.sp(index.TREE.long_plot, bysite, "all", occ.all.output_plot, list.index_site)
index.PLOT.long_site <- index.PLOT.long.merge(index.PLOT.long_BFF.site, index.PLOT.long_GHFF.site, index.PLOT.long_LRFF.site, index.PLOT.long_all.site)
index.PLOT.wide_site <- index.PLOT.wide(index.PLOT.long_site)
count.PLOT.wide_site <- count.PLOT.wide(heights.subset, list.count_site)
count.PLOT.long_site <- count.PLOT.long(count.PLOT.wide_site, bysite)
kernel.PLOT.long_site <- kernel.PLOT.long(kernel.TREE.long_plot, list.pixel_site, occ.BFF.output_plot, occ.GHFF.output_plot, occ.LRFF.output_plot, occ.all.output_plot)
kernel.PLOT.wide_site <- kernel.PLOT.wide(kernel.PLOT.long_site)
kernelnz.PLOT.long_site <- kernel.PLOT.long(kernelnz.TREE.long_plot, list.pixel_site, occ.BFF.output_plot, occ.GHFF.output_plot, occ.LRFF.output_plot, occ.all.output_plot)
kernelnz.PLOT.wide_site <- kernel.PLOT.wide(kernelnz.PLOT.long_site)

##---------------------------------------- Tree as lowest replicate ---------------------------------------- 
bysite <- c("site.code", "session", "site.accession")
index.TREE.Ewide_site <- index.TREE.Ewide(treebat, bysite)
index.TREE.wide_site <- index.TREE.wide(index.TREE.Ewide_site, bysite, list.count_site) 
heights.subset <- calculate.heightdiff(treebat, bytree) %>% 
  dplyr::left_join(
    tree.tessellation[,c("tree.accession", "site.code", "subplot","value_dirichlet")], 
    by = c("tree.accession", "site.code", "subplot")) %>% 
  mutate(vert.density = count/(height.diff*value_dirichlet)) %>% 
  mutate(vert.density = ifelse(vert.density==Inf,NA,vert.density))
count.TREE.wide_site <- count.TREE.wide(heights.subset, bysite, list.count_site) 
count.TREE.long_site <- count.TREE.long(count.TREE.wide_site, bysite) 

################################################################################################
##-----------------------------------------Save data------------------------------------------##
################################################################################################

write.csv(roostbat, "Data/Processed/roostbat.csv")
write.csv(kernel.TREE.wide_plot, "Data/Processed/kernel.TREE.wide_plot.csv")
write.csv(kernelnz.TREE.wide_plot, "Data/Processed/kernelnz.TREE.wide_plot.csv")
write.csv(heights.subset, "Data/Processed/heights.subset.csv")
write.csv(treebat, "Data/Processed/treebat.csv")
write.csv(index.TREE.wide_plot, "Data/Processed/index.TREE.wide_plot.csv")
write.csv(index.TREE.long_plot, "Data/Processed/index.TREE.long_plot.csv")

saveRDS(roostbat, "Data/Processed/roostbat.RDS")
saveRDS(kernel.TREE.wide_plot, "Data/Processed/kernel.TREE.wide_plot.RDS")
saveRDS(kernelnz.TREE.wide_plot, "Data/Processed/kernelnz.TREE.wide_plot.RDS")
saveRDS(heights.subset, "Data/Processed/heights.subset.RDS")
saveRDS(treebat, "Data/Processed/treebat.RDS")
saveRDS(index.TREE.wide_plot, "Data/Processed/index.TREE.wide_plot.RDS")
saveRDS(index.TREE.long_plot, "Data/Processed/index.TREE.long_plot")

## Note that models will not run from these saved .csv files - formatting (e.g. as factors) will be incorrect when read in again
## Use saved .csv files for figure generation only. For running models, run this entire script and use generated files, or load RDS files
