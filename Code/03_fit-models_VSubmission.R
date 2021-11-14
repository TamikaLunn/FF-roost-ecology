## Title: R script for analysing roosting structure of flying-fox roosts in SE QLD and NE NSW
## Manuscript: Conventional wisdom on roosting behavior of Australian flying-foxes—A critical review, and evaluation using new data <https://doi.org/10.1002/ece3.8079>
## Author: Tamika Lunn, Griffith University
## Version: VSubmission, created 13 November 2021

## V1-7 - manuscript preparation
## VSubmission - code copied from 'Bat quantitative analysis_FF#1_V7-revision.R'

rm(list=ls())

##############################################################################################
##---------------------------------Load data & set functions--------------------------------##
##############################################################################################

source ("Code/00_functions_VSubmission.R") ## Read from relative path
library(mgcv) #for GAMs
library(binom)

treebat <-readRDS("Data/Processed/treebat.RDS")
heights.subset <- readRDS("Data/Processed/heights.subset.RDS")
##index.TREE.wide_plot <- readRDS("Data/Processed/index.TREE.wide_plot.RDS")
index.TREE.long_plot
centroids <- read.csv("Data/Raw/Centroids.csv") %>%
  mutate(roost.centroid.N = ifelse(site.accession=="DLIS010",NA,roost.centroid.N)) %>% ## Remove May 2019 Lismore area point (DLIS010)
  mutate(roost.centroid.E = ifelse(site.accession=="DLIS010",NA,roost.centroid.E)) 
centroids$session <- as.factor(centroids$session)

##############################################################################################
##----------------------------------- Overview of models: ----------------------------------##
##############################################################################################

## Description on fixed vs random factors here: https://www.theanalysisfactor.com/specifying-fixed-and-random-factors-in-mixed-models/
## Random effects in GAMMs here: https://stat.ethz.ch/R-manual/R-patched/library/mgcv/html/random.effects.html
## Modeling seasonal time series with GAM (with formal equation) + nested structures (with formal equation) with GAM : https://petolau.github.io/Analyzing-double-seasonal-time-series-with-GAM-in-R/
## Modeling seasonal time series (with formal equation) + nested structures with GAM: https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
## ^^ Specify non-independent (nested) random affects by inclusion of autoregressive model (AR) for errors in the model
## General overview of smoothers and GAMs here: file:///C:/Users/s5083936/Downloads/The_Effect_of_Concurvity_in_Generalized_Additive.9.pdf

## Spline choices (bs=):
## 'tp' = Thin plate regression splines. The default smooth for s terms, do not have 'knots' 
## 'ds' = Duchon splines. These generalize thin plate splines
## 'cc ' = Cyclic cubic regression splines. Splines whose ends line up (no discontinuity between start and end. Specify number of occilations with knots (k= )
## 'cr' or 'cs' = Cubic regression splines
## 'ps' = P-splines
## 're' = Random effects
## Knots = the places where the polynomial pieces connect
## A smoother is a mathematics technique for approximating an observed variable Y by a smooth function of one (or several) independent variable(s) X - Smooth functions define the relation between the transformed mean response and the independent variables.

## On specifying family term:
## Poisson with log link (for zero-inflated count data)
## Gamma with log link (for non zero-inflated count data, or continuous data where the variance increase with the square of the mean)
## Binomial with logit-link (for binary data or proportion data) but note - binomial for proportion data ONLY if specified in the correct way: must be genuine binomial distributions (as in, what you've recorded is the number of "successes" x out of "trials" n). For other proportional data (e.g. proportion of area occupied) use a gaussian distribution with a transformation (log or logit) 

## On fitting and checking models:
## gam or gamm is used to fit general additive models. gamm can fit AR models. 
## test significance of random affects by a likelihood ratio test, comparing the fit of full and reduced models. Can do this with anova (see p 21/43 from BBolker_GLM FAQ)
## Check output of fitted model with gam.check(modelname) or gam.check(modelname$gam). Can compare full and reduced models with anova in same way so long as they're nested models: https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
## from summary(gam) will return: EDF (estimated degrees of freedom) - can be interpreted like how much given variable is smoothed (higher EDF value implies more complex splines. If too high the model is overfitted). GCV score (Generalized Cross Validation score) is also good criterion to choose optimal model among a set of fitted models (i.e. to choose an optimal smoothing parameter) - low values are better. P-values: statistical significance of given variable to response variable, tested by F-test (lower is better). R2 (adjusted R-square) indicated how much of the variabiltiy in the data is explained by the model - higher is better. 
## summary(gam)$r.sq #returns R2 value
## summary(gam)$s.table #returns EDF value and p value
## summary(gam)$sp.criterion #returns GCV score 
## AIC(gam) #returns the AIC value of the fitted model


##############################################################################################
##------------------------------------ Overview of data: -----------------------------------##
##############################################################################################

## See data_README

################################################################################################
##-----------------------------------------Start code-----------------------------------------##
################################################################################################

####------------------------------------------------------------------------------------------##
##----------------------------- Density of core/peripheral plots -----------------------------##
####------------------------------------------------------------------------------------------##

##identify sessions there were ANY bats:
bysite <- c("site.code", "session", "site.accession")
occ.all.output_site <- occ.all(treebat, bysite)
#occ.BFF.output_site <- occ.BFF(treebat, bysite)
#occ.GHFF.output_site <- occ.GHFF(treebat, bysite)
#occ.LRFF.output_site <- occ.LRFF(treebat, bysite)

threshold <- 0.8 #set threhold
core.plots <- index.TREE.wide_plot %>%
  filter(species == "all") %>% #extract combined species measure only
  filter(site.accession %in% occ.all.output_site$site.accession) %>% #choose surveys when at least 1 bat was present, only
  mutate(plot.occ = ifelse(occ==0, 0,
                           ifelse(occ>0, 1,
                                  NA))) %>%
  ddply(c("site.code", "subplot"), summarise,
        prop.occ = sum(plot.occ)/sum(!is.na(plot.occ))) %>% #calculate proportion of times each plot was occupied. Use this to assign core vs peripheral rating
  mutate(occupancy.cat = ifelse(prop.occ>=threshold,"Core","Peripheral")) #set threshold for "core" occupancy

## Identify plots there were ANY bats
byplot <- c("site.code", "session", "site.accession", "subplot", "rep")
occ.all.output_plot <- occ.all(treebat, byplot)
#occ.BFF.output_plot <- occ.BFF(treebat, byplot)
#occ.GHFF.output_plot <- occ.GHFF(treebat, byplot)
#occ.LRFF.output_plot <- occ.LRFF(treebat, byplot)

###########################################
#### VS total number of bats per plot #####
###########################################

data <- index.TREE.wide_plot %>%
  full_join(core.plots, by = c("site.code", "subplot")) %>% #join with occupancy category
  filter(species == "all") %>% #extract combined species measure only
  filter(rep %in% occ.all.output_plot$rep) %>% #choose *plots* when at least 1 bat was present, only
  mutate(N = as.numeric(as.character(N))) %>%
  mutate(session = as.numeric(as.character(session))) %>%
  mutate(subplot = as.factor(subplot)) %>%
  mutate(occupancy.cat = as.factor(occupancy.cat)) %>%
  mutate(site.code = as.factor(site.code)) %>%
  #mutate(N = log(N)) %>% 
  create.Site() 
data <- data[!(data$site.code=="DCLU"),] #remove Clunes because it only has core plots, can't run core/peri comparison

## View data:
jpeg("Output/Model outputs/Density_C-P plots/vs bats per plot/data.jpg", width=1200, height=800)
p1 <- ggplot(data, aes(session, N)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,1000)) + facet_grid(.~occupancy.cat)  #overall there is a cyclical affect of session on N
p2 <- ggplot(data, aes(session, N)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,1000)) + facet_grid(site.code~occupancy.cat) #but this isn't seen at the site level 
ggarrange(p1, p2)
dev.off()

## Fit models:
gam_nested_full <- gamm(N ~ occupancy.cat + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot, p=1), method = "REML", data=data, family=Gamma(link="log"))
gam_nested_full_seasonal <- gamm(N ~ occupancy.cat + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot, p=1), method = "REML", data=data, family=Gamma(link="log")) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Check and save model fits:
jpeg("Output/Model outputs/Density_C-P plots/vs bats per plot/Fit_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density_C-P plots/vs bats per plot/Fit nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Density_C-P plots/vs bats per plot/Fit_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density_C-P plots/vs bats per plot/Fit nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full,"Peripheral subplots", "Nested General Additive Model with a-seasonal sessional term")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Density_C-P plots/vs bats per plot/output_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal,"Peripheral subplots", "Nested General Additive Model with seasonal sessional term")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Density_C-P plots/vs bats per plot/output_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Density_C-P plots/vs bats per plot/Random effects_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()

###########################################
##### VS proportion of occupied trees #####
###########################################

data <- index.TREE.wide_plot %>%
  full_join(core.plots, by = c("site.code", "subplot")) %>% #join with occupancy category
  filter(species == "all") %>% #extract combined species measure only
  filter(occ >0) %>% #filter unoccupied plots
  mutate(occ = as.numeric(as.character(occ))) %>%
  mutate(prop.occ = occ/tree.count) %>%
  mutate(session = as.numeric(as.character(session))) %>%
  mutate(subplot = as.factor(subplot)) %>%
  mutate(occupancy.cat = as.factor(occupancy.cat)) %>%
  mutate(site.code = as.factor(site.code)) %>%
  #mutate(prop.occ = sqrt(prop.occ))  %>%
  create.Site()
data <- data[!(data$site.code=="DCLU"),] #remove Clunes because it only has core plots, can't run core/peri comparison

####---------- Model comparison ---------- ####

## View data:
jpeg("Output/Model outputs/Density_C-P plots/vs prop tree occ/data.jpg", width=1200, height=800)
p1 <- ggplot(data, aes(session, prop.occ)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,1)) + facet_grid(.~occupancy.cat) #overall there is a cyclical affect of session on N
p2 <- ggplot(data, aes(session, prop.occ)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,1)) + facet_grid(site.code~occupancy.cat) #but this isn't seen at the site level 
ggarrange(p1, p2)
dev.off()

## Fit models:
gam_nested_full <- gamm(prop.occ ~ occupancy.cat + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot, p=1), method = "REML", data=data, family=gaussian(link="log")) #family=binomial(link="logit")
gam_nested_full_seasonal <- gamm(prop.occ ~ occupancy.cat + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot, p=1), method = "REML", data=data, family=gaussian(link="log")) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Check and save model fits:
jpeg("Output/Model outputs/Density_C-P plots/vs prop tree occ/Fit_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density_C-P plots/vs prop tree occ/Fit nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Density_C-P plots/vs prop tree occ/Fit_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density_C-P plots/vs prop tree occ/Fit nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full,"Peripheral subplots", "Nested General Additive Model with a-seasonal sessional term")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Density_C-P plots/vs prop tree occ/output_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal,"Peripheral subplots", "Nested General Additive Model with seasonal sessional term")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Density_C-P plots/vs prop tree occ/output_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Density_C-P plots/vs prop tree occ/Random effects_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()

####------------------------------------------------------------------------------------------##
##----------------------------- Density of core/peripheral trees -----------------------------##
####------------------------------------------------------------------------------------------##

##identify sessions there were ANY bats:
bysite <- c("site.code", "session", "site.accession")
occ.all.output_site <- occ.all(treebat, bysite)
#occ.BFF.output_site <- occ.BFF(treebat, bysite)
#occ.GHFF.output_site <- occ.GHFF(treebat, bysite)
#occ.LRFF.output_site <- occ.LRFF(treebat, bysite)

### Identify trees as core or peripheral:
threshold <- 0.8
core.trees <- treebat %>%
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed. This is important to do so that proportion is calculated correctly
  mutate(all.index.weight = BFF.index.weight+GHFF.index.weight+LRFF.index.weight) %>%
  filter(site.accession %in% occ.all.output_site$site.accession) %>% #choose surveys when at least 1 bat was present, only (I think this is more appropriate for identifying core/peri trees, over occupied plots
  mutate(tree.occ = ifelse(all.index.weight==0, 0, #Assign binary value to show if tree is occupied
                           ifelse(all.index.weight>0, 1,
                                  NA))) %>% #nrow of tree.occ should be 2522 trees (this is how many were tagged in total)
  ddply(c("site.code", "tree.accession"), summarise, #for each tree calculate the number of times it was occupied, across surveys where there was one bat present at roost
        tree.occ = sum(tree.occ)/sum(!is.na(tree.occ))) %>%
  mutate(occupancy.cat = ifelse(tree.occ>=threshold,"Core","Peripheral")) #set threshold for "core" occupancy


###########################################
### VS number of bats per occupied tree ###
###########################################

data <- treebat %>% 
  dplyr::select(-c(N)) %>% #remove northing because N also used for total trees
  full_join(index.TREE.wide_plot[which(index.TREE.wide_plot$species=="all"),], by = c("site.code", "session", "site.accession", "subplot", "rep")) %>% #join with plot level data
  full_join(core.trees, by = c("site.code", "tree.accession")) %>% #join with tree occupancy category
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed
  create.Site() %>%
  mutate(all.index.weight = BFF.index.weight+GHFF.index.weight+LRFF.index.weight) %>%
  filter(all.index.weight>0) %>% #select occupied trees only
  mutate(session = as.numeric(as.character(session))) %>%
  mutate(subplot = as.factor(subplot)) %>%
  mutate(occupancy.cat = as.factor(occupancy.cat)) %>%
  mutate(site.code = as.factor(site.code)) %>%
  #mutate(all.index.weight = log(all.index.weight))  %>%
  create.Site()

####---------- Model comparison ---------- ####
## View data:
jpeg("Output/Model outputs/Density_C-P trees/vs bats per occ tree/data.jpg", width=1200, height=800)
p1 <- ggplot(data, aes(session, all.index.weight)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,100)) + facet_grid(.~occupancy.cat) #overall there is a cyclical affect of session on N
p2 <- ggplot(data, aes(session, all.index.weight)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,100)) + facet_grid(site.code~occupancy.cat) #but this isn't seen at the site level 
ggarrange(p1, p2)
dev.off()

## Fit models:
gam_nested_full <- gamm(all.index.weight ~ occupancy.cat + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data, family=Gamma(link="log"))
gam_nested_full_seasonal <- gamm(all.index.weight ~ occupancy.cat + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data, family=Gamma(link="log")) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Check and save model fits:
jpeg("Output/Model outputs/Density_C-P trees/vs bats per occ tree/Fit_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density_C-P trees/vs bats per occ tree/Fit nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Density_C-P trees/vs bats per occ tree/Fit_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density_C-P trees/vs bats per occ tree/Fit nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full,"Peripheral trees", "Nested General Additive Model with a-seasonal sessional term")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Density_C-P trees/vs bats per occ tree/output_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal,"Peripheral trees", "Nested General Additive Model with seasonal sessional term")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Density_C-P trees/vs bats per occ tree/output_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Density_C-P trees/vs bats per occ tree/Random effects_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()

####------------------------------------------------------------------------------------------##
##------------------------- Density vs distance from the roost center ------------------------##
####------------------------------------------------------------------------------------------##

### Combine data summaries with additional centroid data
treebat <- full_join(centroids, treebat, by = c("site.code", "session", "site.accession"))

### Calculate euclidean distance between trees and roost centroids (note - centroids vary per session) (euclidean because E/N are a projected coordinate system, not geographical, so the curvature of the earth is accounted for within the UTM zone
temp <- treebat %>% 
  ddply(c("site.code", "session", "tree.accession","site.accession", "subplot", "rep"), summarise,
        eud.distance_center = sqrt(((E-roost.centroid.E)^2)+((N-roost.centroid.N)^2))) # **approximation of ED. Very close to values calculated manually in QGIS (e.g. DRED01001_01 was 29.94 in QGIS, 29.93 in R. DRED01002_01 was 33.94 and 33.93) (see https://math.stackexchange.com/questions/738529/distance-between-two-points-in-utm-coordinates)
temp_fac <- temp %>% 
  ddply(c("site.code", "session","site.accession"), summarise,
        max_eud.distance = max(eud.distance_center)) #note - values are scaled by the maximum distance value *observed* per *session*. Not the maximum possible value. Chose this because the extent of the roost will vary over sessions, so the meaningful measure should be the extent at the time of the survey (i.e. where are bats relative to the perimeter at the time of the survey)
treebat_firstjoin <- full_join(temp, treebat, by = c("site.code", "session", "tree.accession","site.accession", "subplot", "rep")) 
treebat <- full_join(treebat_firstjoin, temp_fac, by = c("site.code", "session","site.accession")) %>%
  mutate(scaled_eud.distance=eud.distance_center/max_eud.distance)

## Identify occupied site.accessions:
#byplot <- c("site.code", "session", "site.accession", "subplot", "rep")
bysite <- c("site.code", "session", "site.accession")
occ.all.output_site <- occ.all(treebat, bysite)
occ.BFF.output_site <- occ.BFF(treebat, bysite)
occ.GHFF.output_site <- occ.GHFF(treebat, bysite)
occ.LRFF.output_site <- occ.LRFF(treebat, bysite)

###########################################
#### VS total number of bats per plot #####
###########################################

scaled_treebat <- treebat %>% 
  filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>% #remove rare cases where trees missed in survey OR were removed
  ddply(c("site.code", "session","site.accession", "subplot", "rep"), summarise,
        scaled_eud.distance_mean = mean(scaled_eud.distance)) #average per tree across plot
index.TREE.long_plot_data <- full_join(scaled_treebat, index.TREE.long_plot, by = c("site.code", "session", "site.accession", "subplot")) ##long data format because need to include species as a variable in the model 

## For data divided by species, filter sessions where species didn't occur:
data_BFF <- index.TREE.long_plot_data %>% 
  filter(species == "BFF")  %>% 
  filter(site.accession %in% occ.BFF.output_site$site.accession) 
data_GHFF <- index.TREE.long_plot_data %>% 
  filter(species == "GHFF")  %>% 
  filter(site.accession %in% occ.GHFF.output_site$site.accession) 
data_LRFF <- index.TREE.long_plot_data %>% 
  filter(species == "LRFF")  %>%
  filter(site.accession %in% occ.LRFF.output_site$site.accession) 
data_spp <- rbind(data_BFF, data_GHFF, data_LRFF) %>% 
  filter(measure == "N") %>% #Note N is estimated from weighted index value per tree, not true count per tree
  mutate(N = round(value, digits=0)) %>% #rename with more informative column name. Round to nearest integer for Poisson distribution 
  filter(scaled_eud.distance_mean!= "NA")  %>% #remove sessions where bats weren't present in roost (i.e. no centroid)
  mutate(session = as.numeric(as.character(session))) %>%
  mutate(subplot = as.factor(subplot)) %>%
  mutate(species = as.factor(species)) %>%
  mutate(measure = as.factor(measure)) %>%
  mutate(site.code = as.factor(site.code)) 

## For data combined by species, filter sessions where bats didn't occur:
data_all <- index.TREE.long_plot_data %>% 
  filter(species == "all")  %>% 
  filter(site.accession %in% occ.all.output_site$site.accession) 

## Prep data for models (format for measure of interest)
data_BFF <- data_spp %>%
  filter(species == "BFF")
data_GHFF <- data_spp %>%
  filter(species == "GHFF")
data_LRFF <- data_spp %>%
  filter(species == "LRFF")

data_all <- data_all %>% #data with combined species
  filter(measure == "N") %>%  #Note N is estimated from weighted index value per tree, not true count per tree
  mutate(N = round(value, digits=0)) %>% #rename with more informative column name. Round to nearest integer for Poisson distribution 
  filter(scaled_eud.distance_mean!= "NA") %>% #remove sessions where bats weren't present in roost (i.e. no centroid)
  mutate(session = as.numeric(as.character(session))) %>%
  mutate(subplot = as.factor(subplot)) %>%
  mutate(species = as.factor(species)) %>%
  mutate(measure = as.factor(measure)) %>%
  mutate(site.code = as.factor(site.code))

####---------- Model fitting ---------- ####
## *note that error distribution has changed to Poisson, because we haven't filtered out unoccupied subplots (only unoccupied roosts). Data is inflated with zeros
## *note that species is split into separate models, because the affect of time has a different seasonal cycle for the different species

## View data:
jpeg("Output/Model outputs/Density v centre/vs bats per plot/data.jpg", width=1200, height=800)
p3 <- ggplot(data_spp, aes(session, N)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,250)) 
p1 <- ggplot(data_spp, aes(session, N)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,250)) + facet_grid(.~species) #overall there is a cyclical affect of session on N
p2 <- ggplot(data_spp, aes(session, N)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,500)) + facet_grid(site.code~species) #but this isn't seen at the site level 
ggarrange(p3, p1)
dev.off()
jpeg("Output/Model outputs/Density v centre/vs bats per plot/data_by site.jpg", width=1200, height=800)
plot(p2)
dev.off()

####---------- Combined species ---------- ####
gam_nested_full <- gamm(N ~ scaled_eud.distance_mean + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot, p=1), method = "REML", data=data_all, family=poisson(link=log))
gam_nested_full_seasonal <- gamm(N ~ scaled_eud.distance_mean + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot, p=1), method = "REML", data=data_all, family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Check and save model fits:
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Distance from centroid", "Nested General Additive Model with a-seasonal sessional term")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Density v centre/vs bats per plot/output_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Distance from centroid","Nested General Additive Model with seasonal sessional term")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Density v centre/vs bats per plot/output_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Random effects_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()

####---------- BFF only ---------- ####
gam_nested_full <- gamm(N ~ scaled_eud.distance_mean + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot, p=1), method = "REML", data=data_BFF, family=poisson(link=log))
gam_nested_full_seasonal <- gamm(N ~ scaled_eud.distance_mean + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot, p=1), method = "REML", data=data_BFF, family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change
## note that if species are run as a fixed factor, AR term would need to be different as not all species are present in each session|site/subplot breakdown (so model will not run)

## Check and save model fits:
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit BFF_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit BFF nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit BFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit BFF nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Distance from centroid","Nested General Additive Model with a-seasonal sessional term - Black flying-fox")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Density v centre/vs bats per plot/output BFF_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Distance from centroid","Nested General Additive Model with seasonal sessional term - Black flying-fox")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Density v centre/vs bats per plot/output BFF_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Random effects BFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()


####---------- GHFF only ---------- ####
gam_nested_full <- gamm(N ~ scaled_eud.distance_mean + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot, p=1), method = "REML", data=data_GHFF, family=poisson(link=log))
gam_nested_full_seasonal <- gamm(N ~ scaled_eud.distance_mean + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot, p=1), method = "REML", data=data_GHFF, family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change
## note that if species are run as a fixed factor, AR term would need to be different as not all species are present in each session|site/subplot breakdown (so model will not run)

## Check and save model fits:
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit GHFF_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit GHFF nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit GHFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit GHFF nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Distance from centroid","Nested General Additive Model with a-seasonal sessional term - Black flying-fox")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Density v centre/vs bats per plot/output GHFF_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Distance from centroid","Nested General Additive Model with seasonal sessional term - Grey-headed flying-fox")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Density v centre/vs bats per plot/output GHFF_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Random effects GHFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()

####---------- LRFF only ---------- ####
gam_nested_full <- gamm(N ~ scaled_eud.distance_mean + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot, p=1), method = "REML", data=data_LRFF, family=poisson(link=log))
gam_nested_full_seasonal <- gamm(N ~ scaled_eud.distance_mean + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot, p=1), method = "REML", data=data_LRFF, family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change
## note that if species are run as a fixed factor, AR term would need to be different as not all species are present in each session|site/subplot breakdown (so model will not run)

## Check and save model fits:
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit LRFF_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit LRFF nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit LRFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Fit LRFF nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Distance from centroid","Nested General Additive Model with a-seasonal sessional term - Little red flying-fox")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Density v centre/vs bats per plot/output LRFF_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Distance from centroid","Nested General Additive Model with seasonal sessional term - Little red flying-fox")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Density v centre/vs bats per plot/output LRFF_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Density v centre/vs bats per plot/Random effects LRFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()

###########################################
##### VS proportion of occupied trees #####
###########################################
scaled_treebat <- treebat %>% 
  filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>% #remove rare cases where trees missed in survey OR were removed
  ddply(c("site.code", "session","site.accession", "subplot", "rep"), summarise,
        scaled_eud.distance_mean = mean(scaled_eud.distance)) #average per tree across plot
index.TREE.long_plot_data <- full_join(scaled_treebat, index.TREE.long_plot, by = c("site.code", "session", "site.accession", "subplot")) ##long data format because need to include species as a variable in the model 

## For data divided by species, filter sessions where species didn't occur:
data_BFF <- index.TREE.long_plot_data %>% 
  filter(species == "BFF")  %>% 
  filter(site.accession %in% occ.BFF.output_site$site.accession) 
data_GHFF <- index.TREE.long_plot_data %>% 
  filter(species == "GHFF")  %>% 
  filter(site.accession %in% occ.GHFF.output_site$site.accession) 
data_LRFF <- index.TREE.long_plot_data %>% 
  filter(species == "LRFF")  %>%
  filter(site.accession %in% occ.LRFF.output_site$site.accession) 
data_spp <- rbind(data_BFF, data_GHFF, data_LRFF) %>% 
  filter(measure == "occ") %>% 
  mutate(occ = as.numeric(as.character(value))) %>% #rename with more informative column name. 
  mutate(prop.occ = occ/tree.count) %>%
  #filter(scaled_eud.distance_mean!= "NA")  %>% #remove sessions where bats weren't present in roost (i.e. no centroid)
  mutate(session = as.numeric(as.character(session))) %>%
  mutate(subplot = as.factor(subplot)) %>%
  mutate(species = as.factor(species)) %>%
  mutate(measure = as.factor(measure)) %>%
  mutate(site.code = as.factor(site.code)) 

## For data combined by species, filter sessions where bats didn't occur:
data_all <- index.TREE.long_plot_data %>% 
  filter(species == "all")  %>% 
  filter(site.accession %in% occ.all.output_site$site.accession) 

## Prep data for models (format for measure of interest)
data_BFF <- data_spp %>% 
  filter(species == "BFF") 
data_GHFF <- data_spp %>% 
  filter(species == "GHFF") 
data_LRFF <- data_spp %>% 
  filter(species == "LRFF") 

data_all <- data_all %>% #data with combined species
  filter(measure == "occ") %>% 
  mutate(occ = as.numeric(as.character(value))) %>% #rename with more informative column name. 
  mutate(prop.occ = occ/tree.count) %>%
  #filter(scaled_eud.distance_mean!= "NA") %>% #remove sessions where bats weren't present in roost (i.e. no centroid)
  mutate(session = as.numeric(as.character(session))) %>%
  mutate(subplot = as.factor(subplot)) %>%
  mutate(species = as.factor(species)) %>%
  mutate(measure = as.factor(measure)) %>%
  mutate(site.code = as.factor(site.code))

####---------- Model fitting ---------- ####
## *note that species is split into separate models, because the affect of time has a different seasonal cycle for the different species

## View data:
jpeg("data.jpg", width=1200, height=800)
p3 <- ggplot(data_spp, aes(session, prop.occ)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,1)) #overall there is a cyclical affect of session on N
p1 <- ggplot(data_spp, aes(session, prop.occ)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,1)) + facet_grid(.~species) #overall there is a cyclical affect of session on N
p2 <- ggplot(data_spp, aes(session, prop.occ)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,1)) + facet_grid(site.code~species) #but this isn't seen at the site level 
ggarrange(p3, p1)
dev.off()
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/data_by site.jpg", width=1200, height=800)
plot(p2)
dev.off()

####---------- Combined species ---------- ####
gam_nested_full <- gamm(prop.occ ~ scaled_eud.distance_mean + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot, p=1), method = "REML", data=data_all, family=gaussian(link="identity"))
gam_nested_full_seasonal <- gamm(prop.occ ~ scaled_eud.distance_mean + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot, p=1), method = "REML", data=data_all, family=gaussian(link="identity")) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Check and save model fits:
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Distance from centroid", "Nested General Additive Model with a-seasonal sessional term")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Density v centre/vs prop tree occ/output_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Distance from centroid","Nested General Additive Model with seasonal sessional term")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Density v centre/vs prop tree occ/output_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Random effects_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()


####---------- BFF only ---------- ####
gam_nested_full <- gamm(prop.occ ~ scaled_eud.distance_mean + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot + session, p=1), method = "REML", data=data_BFF, family=gaussian(link="identity"))
gam_nested_full_seasonal <- gamm(prop.occ ~ scaled_eud.distance_mean + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot + session, p=1), method = "REML", data=data_BFF, family=gaussian(link="identity")) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change
## note that if species are run as a fixed factor, AR term would need to be different as not all species are present in each session|site/subplot breakdown (so model will not run)

## Check and save model fits:
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit BFF_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit BFF nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit BFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit BFF nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Distance from centroid","Nested General Additive Model with a-seasonal sessional term - Black flying-fox")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Density v centre/vs prop tree occ/output BFF_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Distance from centroid","Nested General Additive Model with seasonal sessional term - Black flying-fox")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Density v centre/vs prop tree occ/output BFF_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Random effects BFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()


####---------- GHFF only ---------- ####
gam_nested_full <- gamm(prop.occ ~ scaled_eud.distance_mean + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot + session, p=1), method = "REML", data=data_GHFF, family=gaussian(link="identity"))
gam_nested_full_seasonal <- gamm(prop.occ ~ scaled_eud.distance_mean + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot + session, p=1), method = "REML", data=data_GHFF,family=gaussian(link="identity")) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change
## note that if species are run as a fixed factor, AR term would need to be different as not all species are present in each session|site/subplot breakdown (so model will not run)

## Check and save model fits:
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit GHFF_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit GHFF nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit GHFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit GHFF nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Distance from centroid","Nested General Additive Model with a-seasonal sessional term - Grey-headed flying-fox")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Density v centre/vs prop tree occ/output GHFF_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Distance from centroid","Nested General Additive Model with seasonal sessional term - Grey-headed flying-fox")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Density v centre/vs prop tree occ/output GHFF_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Random effects GHFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()

####---------- LRFF only ---------- ####
gam_nested_full <- gamm(prop.occ ~ scaled_eud.distance_mean + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot + session, p=1), method = "REML", data=data_LRFF, family=gaussian(link="identity"))
gam_nested_full_seasonal <- gamm(prop.occ ~ scaled_eud.distance_mean + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ session|site.code/subplot + session, p=1), method = "REML", data=data_LRFF, family=gaussian(link="identity")) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change
## note that if species are run as a fixed factor, AR term would need to be different as not all species are present in each session|site/subplot breakdown (so model will not run)

## Check and save model fits:
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit LRFF_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit LRFF nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit LRFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Fit LRFF nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Distance from centroid","Nested General Additive Model with a-seasonal sessional term - Little red flying-fox")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Density v centre/vs prop tree occ/output LRFF_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Distance from centroid","Nested General Additive Model with seasonal sessional term - Little red flying-fox")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Density v centre/vs prop tree occ/output LRFF_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Density v centre/vs prop tree occ/Random effects LRFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()


###########################################
########### VS mixing of sexes ############
###########################################
treebat_long <- treebat %>% 
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed
  mutate(BFF.prop.M = BFF.M/(BFF.M+BFF.F))  %>%  
  mutate(GHFF.prop.M = GHFF.M/(GHFF.M+GHFF.F))  %>%  
  mutate(LRFF.prop.M = LRFF.M/(LRFF.M+LRFF.F))  %>% 
  mutate(all.prop.M = (BFF.M+GHFF.M+LRFF.M)/(BFF.M+GHFF.M+LRFF.M+BFF.F+GHFF.F+LRFF.F))  %>% 
  melt(id.vars = c("site.code", "session", "tree.accession", "site.accession", "subplot"), measure.vars = c("BFF.prop.M", "GHFF.prop.M", "LRFF.prop.M", "all.prop.M"),
       variable.name = c("species.meas"), value.name="value") %>%
  dplyr::mutate(species = str_extract(species.meas, "GHFF|BFF|LRFF|all")) %>%
  dplyr::mutate(measure = str_extract(species.meas, "prop.M")) %>%
  dplyr::select(-c(species.meas)) 

treebat_sub <- treebat %>% #select collumns to keep from merge:
  dplyr::select(c(site.code, session, tree.accession, site.accession, subplot, eud.distance_center, max_eud.distance, scaled_eud.distance, year, month, day)) 
treebat_long_data <- full_join(treebat_long, treebat_sub, by = c("site.code", "session", "site.accession", "tree.accession", "subplot")) ##long data format because need to include species as a variable in the model 

# For data divided by species, filter sessions where species didn't occur:
data_BFF <- treebat_long_data %>% 
  filter(species == "BFF")  %>% 
  filter(site.accession %in% occ.BFF.output_site$site.accession) 
data_GHFF <- treebat_long_data %>% 
  filter(species == "GHFF")  %>% 
  filter(site.accession %in% occ.GHFF.output_site$site.accession) 
data_LRFF <- treebat_long_data %>% 
  filter(species == "LRFF")  %>%
  filter(site.accession %in% occ.LRFF.output_site$site.accession) 
data_spp <- rbind(data_BFF, data_GHFF, data_LRFF) %>% 
  filter(measure == "prop.M") %>% 
  mutate(prop.M = as.numeric(as.character(value))) %>% #rename with more informative column name. 
  mutate(session = as.numeric(as.character(session))) %>%
  mutate(subplot = as.factor(subplot)) %>%
  mutate(species = as.factor(species)) %>%
  mutate(measure = as.factor(measure)) %>%
  mutate(site.code = as.factor(site.code)) %>%
  filter(!is.na(prop.M))

## For data combined by species, filter sessions where bats didn't occur:
data_all <- treebat_long_data %>% 
  filter(species == "all")  %>% 
  filter(site.accession %in% occ.all.output_site$site.accession) 

## Prep data for models (format for measure of interest)
data_BFF <- data_spp  %>% 
  filter(species=="BFF")
data_GHFF <- data_spp  %>% 
  filter(species=="GHFF")
data_LRFF <- data_spp  %>% 
  filter(species=="LRFF")

data_all <- data_all %>% #data with combined species
  filter(measure == "prop.M") %>% 
  mutate(prop.M = as.numeric(as.character(value))) %>% #rename with more informative column name. 
  mutate(session = as.numeric(as.character(session))) %>%
  mutate(subplot = as.factor(subplot)) %>%
  mutate(species = as.factor(species)) %>%
  mutate(measure = as.factor(measure)) %>%
  mutate(site.code = as.factor(site.code)) %>%
  filter(!is.na(prop.M))

####---------- Model fitting ---------- ####
## *note that error distribution has changed to Poisson, because we haven't filtered out unoccupied subplots (only unoccupied roosts). Data is inflated with zeros
## note that AR term is different as not all session|site/subplot combinations have bats recorded (so model will not run)

## View data:
jpeg("Output/Model outputs/Density v centre/vs sex/data.jpg", width=1200, height=800)
p3 <- ggplot(data_spp, aes(session, prop.M)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,1)) #overall there is a cyclical affect of session on N
p1 <- ggplot(data_spp, aes(session, prop.M)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,1)) + facet_grid(.~species) #overall there is a cyclical affect of session on N
p2 <- ggplot(data_spp, aes(session, prop.M)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,1)) + facet_grid(site.code~species) #but this isn't seen at the site level 
ggarrange(p3, p1)
dev.off()
jpeg("Output/Model outputs/Density v centre/vs sex/data_by site.jpg", width=1200, height=800)
plot(p2)
dev.off()

####---------- Combined species ---------- ####
gam_nested_full <- gamm(prop.M ~ scaled_eud.distance + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_all, family=gaussian(link="identity"))
gam_nested_full_seasonal <- gamm(prop.M ~ scaled_eud.distance + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_all, family=gaussian(link="identity")) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Check and save model fits:
jpeg("Output/Model outputs/Density v centre/vs sex/Fit_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs sex/Fit nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Density v centre/vs sex/Fit_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs sex/Fit nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Distance from centroid", "Nested General Additive Model with a-seasonal sessional term")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Density v centre/vs sex/output_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Distance from centroid","Nested General Additive Model with seasonal sessional term")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Density v centre/vs sex/output_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Density v centre/vs sex/Random effects_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()


####---------- BFF only ---------- ####
gam_nested_full <- gamm(prop.M ~ scaled_eud.distance + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_BFF, family=gaussian(link="identity"))
gam_nested_full_seasonal <- gamm(prop.M ~ scaled_eud.distance + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_BFF, family=gaussian(link="identity")) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change
## note that if species are run as a fixed factor, AR term would need to be different as not all species are present in each session|site/subplot breakdown (so model will not run)

## Check and save model fits:
jpeg("Output/Model outputs/Density v centre/vs sex/Fit BFF_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs sex/Fit BFF nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Density v centre/vs sex/Fit BFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs sex/Fit BFF nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Distance from centroid","Nested General Additive Model with a-seasonal sessional term - Black flying-fox")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Density v centre/vs sex/output BFF_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Distance from centroid","Nested General Additive Model with seasonal sessional term - Black flying-fox")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Density v centre/vs sex/output BFF_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Density v centre/vs sex/Random effects BFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()


####---------- GHFF only ---------- ####
gam_nested_full <- gamm(prop.M ~ scaled_eud.distance + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_GHFF, family=gaussian(link="identity"))
gam_nested_full_seasonal <- gamm(prop.M ~ scaled_eud.distance + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_GHFF,family=gaussian(link="identity")) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change
## note that if species are run as a fixed factor, AR term would need to be different as not all species are present in each session|site/subplot breakdown (so model will not run)

## Check and save model fits:
jpeg("Output/Model outputs/Density v centre/vs sex/Fit GHFF_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs sex/Fit GHFF nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Density v centre/vs sex/Fit GHFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs sex/Fit GHFF nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Distance from centroid","Nested General Additive Model with a-seasonal sessional term - Grey-headed flying-fox")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Density v centre/vs sex/output GHFF_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Distance from centroid","Nested General Additive Model with seasonal sessional term - Grey-headed flying-fox")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Density v centre/vs sex/output GHFF_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Density v centre/vs sex/Random effects GHFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()

####---------- LRFF only ---------- ####
gam_nested_full <- gamm(prop.M ~ scaled_eud.distance + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_LRFF, family=gaussian(link="identity"))
gam_nested_full_seasonal <- gamm(prop.M ~ scaled_eud.distance + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_LRFF, family=gaussian(link="identity")) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change
## note that if species are run as a fixed factor, AR term would need to be different as not all species are present in each session|site/subplot breakdown (so model will not run)

## Check and save model fits:
jpeg("Output/Model outputs/Density v centre/vs sex/Fit LRFF_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs sex/Fit LRFF nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Density v centre/vs sex/Fit LRFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Density v centre/vs sex/Fit LRFF nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Distance from centroid","Nested General Additive Model with a-seasonal sessional term - Little red flying-fox")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Density v centre/vs sex/output LRFF_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Distance from centroid","Nested General Additive Model with seasonal sessional term - Little red flying-fox")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Density v centre/vs sex/output LRFF_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Density v centre/vs sex/Random effects LRFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()


####------------------------------------------------------------------------------------------##
##--------------------------- Roosting height vs count per species ---------------------------##
####------------------------------------------------------------------------------------------##

## Prep data for models (format for measure of interest)
data_spp <- heights.subset %>%
  filter(species=="BFF"|species=="GHFF"|species=="LRFF")  %>%
  mutate(session = as.numeric(as.character(session))) %>%
  mutate(subplot = as.factor(subplot)) %>%
  mutate(species = as.factor(species)) %>%
  mutate(tree.accession = as.factor(tree.accession)) %>%
  mutate(max = as.numeric(as.character(max))) %>%
  mutate(site.code = as.factor(site.code))  %>%
  filter(!is.na(max)) %>%
  filter(!is.na(count))
summary(data_spp)

data_BFF <- data_spp %>%
  filter(species == "BFF")
data_GHFF <- data_spp %>%
  filter(species == "GHFF")
data_LRFF <- data_spp %>%
  filter(species == "LRFF")

####---------- Model fitting ---------- ####
## View data:
jpeg("Output/Model outputs/Height v count/data.jpg", width=1200, height=800)
p3 <- ggplot(data_spp, aes(session, max)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,30)) 
p1 <- ggplot(data_spp, aes(session, max)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,30)) + facet_grid(.~species) 
p2 <- ggplot(data_spp, aes(session, max)) + geom_point(col="grey") + geom_smooth(col="black") + theme_bw() + background_grid(major="x", colour.major = "grey95") + coord_cartesian(y=c(0,30)) + facet_grid(site.code~species) 
ggarrange(p3, p1)
dev.off()
jpeg("Output/Model outputs/Height v count/data_by site.jpg", width=1200, height=800)
plot(p2)
dev.off()

####---------- Combined species ---------- ####
gam_nested_full <- gamm(max ~ count + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_spp, family=Gamma(link="log"))
gam_nested_full_seasonal <- gamm(max ~ count + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_spp, family=Gamma(link="log")) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Check and save model fits:
jpeg("Output/Model outputs/Height v count/Fit_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Height v count/Fit nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Height v count/Fit_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Height v count/Fit nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Number of bats in tree", "Nested General Additive Model with a-seasonal sessional term")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Height v count/output_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Number of bats in tree","Nested General Additive Model with seasonal sessional term")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Height v count/output_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Height v count/Random effects_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()


####---------- BFF ---------- ####
gam_nested_full <- gamm(max ~ count + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_BFF, family=Gamma(link="log"))
gam_nested_full_seasonal <- gamm(max ~ count + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_BFF, family=Gamma(link="log")) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Check and save model fits:
jpeg("Output/Model outputs/Height v count/Fit BFF_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Height v count/Fit BFF nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Height v count/Fit BFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Height v count/Fit BFF nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Number of bats in tree", "Nested General Additive Model with a-seasonal sessional term - Black flying-fox")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Height v count/output BFF_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Number of bats in tree","Nested General Additive Model with seasonal sessional term - Black flying-fox")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Height v count/output BFF_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Height v count/Random effects BFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()

####---------- GHFF ---------- ####
gam_nested_full <- gamm(max ~ count + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_GHFF, family=Gamma(link="log"))
gam_nested_full_seasonal <- gamm(max ~ count + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_GHFF, family=Gamma(link="log")) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Check and save model fits:
jpeg("Output/Model outputs/Height v count/Fit GHFF_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Height v count/Fit GHFF nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Height v count/Fit GHFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Height v count/Fit GHFF nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Number of bats in tree", "Nested General Additive Model with a-seasonal sessional term - Grey-headed flying-fox")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Height v count/output GHFF_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Number of bats in tree","Nested General Additive Model with seasonal sessional term - Grey-headed flying-fox")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Height v count/output GHFF_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Height v count/Random effects GHFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()

####---------- LRFF ---------- ####
gam_nested_full <- gamm(max ~ count + s(session, bs="re") + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_LRFF, family=Gamma(link="log"))
gam_nested_full_seasonal <- gamm(max ~ count + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot + session, p=1), method = "REML", data=data_LRFF, family=Gamma(link="log")) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Check and save model fits:
jpeg("Output/Model outputs/Height v count/Fit LRFF_gam_nested_full.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full$gam)
dev.off()
#intervals(gam_nested_full$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Height v count/Fit LRFF nested terms_gam_nested_full.jpg", width=1200, height=800)
pacf(resid(gam_nested_full$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

jpeg("Output/Model outputs/Height v count/Fit LRFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_nested_full_seasonal$gam)
dev.off()
#intervals(gam_nested_full_seasonal$lme,which= "var-cov")$corStruct ##Estimated coefficient of AR(1) process. High values (>0.5) indicates a strong dependency on previous values of errors
jpeg("Output/Model outputs/Height v count/Fit LRFF nested terms_gam_nested_full_seasonal.jpg", width=1200, height=800)
pacf(resid(gam_nested_full_seasonal$lme,type="normalized"),lag.max=48,main="Nested") #Plot of values of a partial autocorrelation function applied to normalized residuals. It gives the partial correlation of residuals with its own lagged values. Optimal values of pACF should be within dashed blue lines
dev.off()

## Compare models:
#anova(gam_nested_full$lme, gam_nested_full_seasonal$lme)

## Evaluate output:
#summary(gam_nested_full$gam)
#summary(gam_nested_full_seasonal$gam)
output_gam_nested_full <- save.output.gamm_2lev(gam_nested_full, "Number of bats in tree", "Nested General Additive Model with a-seasonal sessional term - Little red flying-fox")
saveRDS(output_gam_nested_full, file = "Output/Model outputs/Height v count/output LRFF_gam_nested_full.Rds")
output_gam_nested_full_seasonal <- save.output.gamm_2lev(gam_nested_full_seasonal, "Number of bats in tree","Nested General Additive Model with seasonal sessional term - Little red flying-fox")
saveRDS(output_gam_nested_full_seasonal, file = "Output/Model outputs/Height v count/output LRFF_gam_nested_full_seasonal.Rds")

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Height v count/Random effects LRFF_gam_nested_full_seasonal.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_nested_full_seasonal$gam)
dev.off()

