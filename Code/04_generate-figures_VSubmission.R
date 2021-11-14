## Title: R script for analysing roosting structure of flying-fox roosts in SE QLD and NE NSW
## Manuscript: Conventional wisdom on roosting behavior of Australian flying-foxes—A critical review, and evaluation using new data <https://doi.org/10.1002/ece3.8079>
## Author: Tamika Lunn, Griffith University
## Version: VSubmission, created 13 November 2021

## V1-3 - manuscript preparation
## VSubmission - code copied from 'Figures_FF-1_V3.R'

rm(list=ls())

##############################################################################################
##---------------------------------Load data & set functions--------------------------------##
##############################################################################################

source ("Code/00_functions_VSubmission.R") ## Read from relative path
library(binom)

treebat <-readRDS("Data/Processed/treebat.RDS")
heights.subset <- readRDS("Data/Processed/heights.subset.RDS")
index.TREE.wide_plot <- readRDS("Data/Processed/index.TREE.wide_plot.RDS")
roostbat <- readRDS("Data/Processed/roostbat.RDS")
centroids <- read.csv("Data/Raw/Centroids.csv") %>%
  mutate(roost.centroid.N = ifelse(site.accession=="DLIS010",NA,roost.centroid.N)) %>% ## Remove May 2019 Lismore area point (DLIS010)
  mutate(roost.centroid.E = ifelse(site.accession=="DLIS010",NA,roost.centroid.E)) 
centroids$session <- as.factor(centroids$session)


##############################################################################################
##------------------------------------ Overview of data: -----------------------------------##
##############################################################################################

## See data_README

################################################################################################
##-----------------------------------------Start code-----------------------------------------##
################################################################################################

##-------------------------------------- Quick summaries -------------------------------------##

###################################
#### Check sex ratio in roosts ####
###################################

#### All species, all roost sites ####
treebat_summary <- treebat %>%
  mutate(all.M = BFF.M+GHFF.M + LRFF.M) %>%
  mutate(all.F = BFF.F+GHFF.F + LRFF.F) %>%
  filter(is.na(all.M)==FALSE) %>%
  filter(is.na(all.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","BFF.M", "GHFF.M", "LRFF.M", "all.M", "BFF.F", "GHFF.F", "LRFF.F", "all.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.M = sum(all.M),
        total.F = sum(all.F),
        total.M.F.ratio = total.M/total.F) %>%
  mutate(blank.category = "blah")  %>%
  filter(is.na(total.M.F.ratio)==FALSE) %>%
  #View(treebat_summary[,c("urban","site.code", "session", "total.M", "total.F", "total.M.F.ratio")])
  ## Calculate average ratio across all roosts*sessions
  ddply(c("blank.category"), summarise,
        mean.M.F.ratio = mean(total.M.F.ratio))
## For every 1 female there is 0.76 males

#### All species, by urban/non-urban area ####
treebat_summary <- treebat %>%
  mutate(all.M = BFF.M+GHFF.M + LRFF.M) %>%
  mutate(all.F = BFF.F+GHFF.F + LRFF.F) %>%
  filter(is.na(all.M)==FALSE) %>%
  filter(is.na(all.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","BFF.M", "GHFF.M", "LRFF.M", "all.M", "BFF.F", "GHFF.F", "LRFF.F", "all.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.M = sum(all.M),
        total.F = sum(all.F),
        total.M.F.ratio = total.M/total.F) %>%
  mutate(urban = ifelse(site.code=="DCAN"|site.code=="DCLU", "Non-urban", ## Roosts >10km from the urban edge
                        ifelse(site.code=="DLIS", "Neither", ## Lismore roost was only 0.4km from the urban edge
                               "Urban")))  %>% 
  filter(is.na(total.M.F.ratio)==FALSE) %>%
  #View(treebat_summary[,c("urban","site.code", "session", "total.M", "total.F", "total.M.F.ratio")])
  ddply(c("urban"), summarise,
        ## Calculate average ratio across all urban and non-urban roost types
        mean.M.F.ratio = mean(total.M.F.ratio))
## For every 1 female there is 0.84 males (Non-urban) and 0.73 males (Urban)

#### All species, by contemporary/non-contemporary roosts: ####
treebat_summary <- treebat %>%
  mutate(all.M = BFF.M+GHFF.M + LRFF.M) %>%
  mutate(all.F = BFF.F+GHFF.F + LRFF.F) %>%
  filter(is.na(all.M)==FALSE) %>%
  filter(is.na(all.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","BFF.M", "GHFF.M", "LRFF.M", "all.M", "BFF.F", "GHFF.F", "LRFF.F", "all.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.M = sum(all.M),
        total.F = sum(all.F),
        total.M.F.ratio = total.M/total.F) %>%
  mutate(contemp = ifelse(site.code=="DSUN"|site.code=="DAVO"|site.code=="DBUR", "contemporary", ## see Table 1
                          "non-contemporary"))  %>% 
  filter(is.na(total.M.F.ratio)==FALSE) %>%
  #View(treebat_summary[,c("contemp","site.code", "session", "total.M", "total.F", "total.M.F.ratio")])
  ddply(c("contemp"), summarise,
        ## Calculate average ratio across all contemporary and non-contemporary roosts
        mean.M.F.ratio = mean(total.M.F.ratio))
## For every 1 female there is 0.71 males (contemporary) and 0.79 males (non-contemporary)


#### BFF, all roost sites ####
treebat_summary <- treebat %>%
  filter(is.na(BFF.M)==FALSE) %>%
  filter(is.na(BFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","BFF.M", "BFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.M = sum(BFF.M),
        total.F = sum(BFF.F),
        total.M.F.ratio = total.M/total.F) %>%
  mutate(blank.category = "blah")  %>%
  filter(is.na(total.M.F.ratio)==FALSE) %>%
  #View(treebat_summary[,c("blank.category","site.code", "session", "total.M", "total.F", "total.M.F.ratio")])
  ## Calculate average ratio across all roosts*sessions
  ddply(c("blank.category"), summarise,
        mean.M.F.ratio = mean(total.M.F.ratio))
## For every 1 BFF female there is 0.92 BFF males

#### GHFF, all roost sites ####
treebat_summary <- treebat %>%
  filter(is.na(GHFF.M)==FALSE) %>%
  filter(is.na(GHFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","GHFF.M", "GHFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.M = sum(GHFF.M),
        total.F = sum(GHFF.F),
        total.M.F.ratio = total.M/total.F) %>%
  mutate(blank.category = "blah")  %>%
  filter(is.na(total.M.F.ratio)==FALSE) %>%
  filter(!total.M.F.ratio==Inf) %>%
  #View(treebat_summary[,c("blank.category","site.code", "session", "total.M", "total.F", "total.M.F.ratio")])
  ## Calculate average ratio across all roosts*sessions
  ddply(c("blank.category"), summarise,
        mean.M.F.ratio = mean(total.M.F.ratio))
## For every 1 GHFF female there is 0.64 GHFF males

#### LRFF, all roost sites ####
treebat_summary <- treebat %>%
  filter(is.na(LRFF.M)==FALSE) %>%
  filter(is.na(LRFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","LRFF.M", "LRFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.M = sum(LRFF.M),
        total.F = sum(LRFF.F),
        total.M.F.ratio = total.M/total.F) %>%
  mutate(blank.category = "blah")  %>%
  filter(is.na(total.M.F.ratio)==FALSE) %>%
  filter(!total.M.F.ratio==Inf) %>%
  #View(treebat_summary[,c("blank.category","site.code", "session", "total.M", "total.F", "total.M.F.ratio")])
  ## Calculate average ratio across all roosts*sessions
  ddply(c("blank.category"), summarise,
        mean.M.F.ratio = mean(total.M.F.ratio))
## For every 1 LRFF female there is 1.60 LRFF males


#### BFF, by urban/non-urban area ####
treebat_summary <- treebat %>%
  filter(is.na(BFF.M)==FALSE) %>%
  filter(is.na(BFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","BFF.M", "BFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.M = sum(BFF.M),
        total.F = sum(BFF.F),
        total.M.F.ratio = total.M/total.F) %>%
  mutate(urban = ifelse(site.code=="DCAN"|site.code=="DCLU", "Non-urban", ## Roosts >10km from the urban edge
                        ifelse(site.code=="DLIS", "Neither", ## Lismore roost was only 0.4km from the urban edge
                               "Urban")))  %>% 
  filter(is.na(total.M.F.ratio)==FALSE) %>%
  #View(treebat_summary[,c("urban","site.code", "session", "total.M", "total.F", "total.M.F.ratio")])
  ddply(c("urban"), summarise,
        ## Calculate average ratio across all urban and non-urban roost types
        mean.M.F.ratio = mean(total.M.F.ratio))
## For every 1 BFF female there is 0.77 BFF males (urban) and 1.3 BFF males (non-urban)

#### BFF, by contemporary/non-contemporary roosts ####
treebat_summary <- treebat %>%
  filter(is.na(BFF.M)==FALSE) %>%
  filter(is.na(BFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","BFF.M", "BFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.M = sum(BFF.M),
        total.F = sum(BFF.F),
        total.M.F.ratio = total.M/total.F) %>%
  mutate(contemp = ifelse(site.code=="DSUN"|site.code=="DAVO"|site.code=="DBUR", "contemporary", ## see Table 1
                          "non-contemporary"))  %>% 
  filter(is.na(total.M.F.ratio)==FALSE) %>%
  #View(treebat_summary[,c("contemp","site.code", "session", "total.M", "total.F", "total.M.F.ratio")])
  ddply(c("contemp"), summarise,
        ## Calculate average ratio across all contemporary and non-contemporary roosts
        mean.M.F.ratio = mean(total.M.F.ratio))
## For every 1 BFF female there is 0.74 BFF males (contemporary) and 1.02 BFF males (non-contemporary)

#### GHFF, by urban/non-urban area ####
treebat_summary <- treebat %>%
  filter(is.na(GHFF.M)==FALSE) %>%
  filter(is.na(GHFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","GHFF.M", "GHFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.M = sum(GHFF.M),
        total.F = sum(GHFF.F),
        total.M.F.ratio = total.M/total.F) %>%
  mutate(urban = ifelse(site.code=="DCAN"|site.code=="DCLU", "Non-urban", ## Roosts >10km from the urban edge
                        ifelse(site.code=="DLIS", "Neither", ## Lismore roost was only 0.4km from the urban edge
                               "Urban")))  %>% 
  filter(is.na(total.M.F.ratio)==FALSE) %>%
  filter(!total.M.F.ratio==Inf) %>%
  #View(treebat_summary[,c("urban","site.code", "session", "total.M", "total.F", "total.M.F.ratio")])
  ddply(c("urban"), summarise,
        ## Calculate average ratio across all urban and non-urban roost types
        mean.M.F.ratio = mean(total.M.F.ratio))
## For every 1 GHFF female there is 0.54 GHFF males (urban) and 0.65 GHFF males (non-urban)

#### GHFF, by contemporary/non-contemporary roosts ####
treebat_summary <- treebat %>%
  filter(is.na(GHFF.M)==FALSE) %>%
  filter(is.na(GHFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","GHFF.M", "GHFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.M = sum(GHFF.M),
        total.F = sum(GHFF.F),
        total.M.F.ratio = total.M/total.F) %>%
  mutate(contemp = ifelse(site.code=="DSUN"|site.code=="DAVO"|site.code=="DBUR", "contemporary", ## see Table 1
                          "non-contemporary"))  %>% 
  filter(is.na(total.M.F.ratio)==FALSE) %>%
  filter(!total.M.F.ratio==Inf) %>%
  #View(treebat_summary[,c("contemp","site.code", "session", "total.M", "total.F", "total.M.F.ratio")])
  ddply(c("contemp"), summarise,
        ## Calculate average ratio across all contemporary and non-contemporary roosts
        mean.M.F.ratio = mean(total.M.F.ratio))
## For every 1 GHFF female there is 0.45 GHFF males (contemporary) and 0.69 GHFF males (non-contemporary)


#### LRFF, by urban/non-urban area ####
treebat_summary <- treebat %>%
  filter(is.na(LRFF.M)==FALSE) %>%
  filter(is.na(LRFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","LRFF.M", "LRFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.M = sum(LRFF.M),
        total.F = sum(LRFF.F),
        total.M.F.ratio = total.M/total.F) %>%
  mutate(urban = ifelse(site.code=="DCAN"|site.code=="DCLU", "Non-urban", ## Roosts >10km from the urban edge
                        ifelse(site.code=="DLIS", "Neither", ## Lismore roost was only 0.4km from the urban edge
                               "Urban")))  %>% 
  filter(is.na(total.M.F.ratio)==FALSE) %>%
  filter(!total.M.F.ratio==Inf) %>%
  #View(treebat_summary[,c("urban","site.code", "session", "total.M", "total.F", "total.M.F.ratio")])
  ddply(c("urban"), summarise,
        ## Calculate average ratio across all urban and non-urban roost types
        mean.M.F.ratio = mean(total.M.F.ratio))
## LRFF were only observed in urban roosts. For every 1 LRFF female there is 1.60 BFF males


#### LRFF, by contemporary/non-contemporary roosts ####
treebat_summary <- treebat %>%
  filter(is.na(LRFF.M)==FALSE) %>%
  filter(is.na(LRFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","LRFF.M", "LRFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.M = sum(LRFF.M),
        total.F = sum(LRFF.F),
        total.M.F.ratio = total.M/total.F) %>%
  mutate(contemp = ifelse(site.code=="DSUN"|site.code=="DAVO"|site.code=="DBUR", "contemporary", ## see Table 1
                          "non-contemporary"))  %>% 
  filter(is.na(total.M.F.ratio)==FALSE) %>%
  filter(!total.M.F.ratio==Inf) %>%
  #View(treebat_summary[,c("contemp","site.code", "session", "total.M", "total.F", "total.M.F.ratio")])
  ddply(c("contemp"), summarise,
        ## Calculate average ratio across all contemporary and non-contemporary roosts
        mean.M.F.ratio = mean(total.M.F.ratio))
## For every 1 GHFF female there is 0.39 GHFF males (contemporary) and 1.90 GHFF males (non-contemporary)

##########################################################
#### Reproducing females, filtered to birthing season ####
#### (BFF & GHFF: October-December, LRFF: April-June) ####
##########################################################

#### BFF, October-December, all sites ####
treebat_summary <- treebat %>% 
  filter(session=="3"|session=="4"|session=="5") %>% ##October-December counts only
  filter(is.na(BFF.Fw)==FALSE) %>%
  filter(is.na(BFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","BFF.Fw", "BFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all females and reproducing females across each roost per session, and calculate Fw:F ratio for each site*session
        total.Fw = sum(BFF.Fw),
        total.F = sum(BFF.F),
        total.F.ratio = total.Fw/total.F) %>%
  mutate(blank.category = "blah")  %>%
  filter(is.na(total.F.ratio)==FALSE) %>%
  #View(treebat_summary[,c("blank.category","site.code", "session", "total.Fw", "total.F", "total.F.ratio")])
  ## Calculate average ratio across all roosts*sessions
  ddply(c("blank.category"), summarise,
        mean.F.ratio = mean(total.F.ratio))
## For every 1 BFF female there is 0.44 BFF females with young

#### GHFF, October-December, all sites ####
treebat_summary <- treebat %>%
  filter(session=="3"|session=="4"|session=="5") %>% ##October-December counts only
  filter(is.na(GHFF.Fw)==FALSE) %>%
  filter(is.na(GHFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","GHFF.Fw", "GHFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.Fw = sum(GHFF.Fw),
        total.F = sum(GHFF.F),
        total.Fw.F.ratio = total.Fw/total.F) %>%
  mutate(blank.category = "blah")  %>%
  filter(is.na(total.Fw.F.ratio)==FALSE) %>%
  filter(!total.Fw.F.ratio==Inf) %>%
  #View(treebat_summary[,c("blank.category","site.code", "session", "total.Fw", "total.F", "total.Fw.F.ratio")])
  ## Calculate average ratio across all roosts*sessions
  ddply(c("blank.category"), summarise,
        mean.Fw.F.ratio = mean(total.Fw.F.ratio))
## For every 1 GHFF female there is 0.44 females with young

#### LRFF, April-June, all sites ####
treebat_summary <- treebat %>%
  filter(session=="9"|session=="10"|session=="11") %>% ##October-December counts only
  filter(is.na(LRFF.Fw)==FALSE) %>%
  filter(is.na(LRFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","LRFF.Fw", "LRFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.Fw = sum(LRFF.Fw),
        total.F = sum(LRFF.F),
        total.Fw.F.ratio = total.Fw/total.F) %>%
  mutate(blank.category = "blah")  %>%
  filter(is.na(total.Fw.F.ratio)==FALSE) %>%
  filter(!total.Fw.F.ratio==Inf) %>%
  #View(treebat_summary[,c("blank.category","site.code", "session", "total.Fw", "total.F", "total.Fw.F.ratio")])
  ## Calculate average ratio across all roosts*sessions
  ddply(c("blank.category"), summarise,
        mean.Fw.F.ratio = mean(total.Fw.F.ratio))
## We did not observe LRFF females with pups

#### BFF, by urban/non-urban area, October-December ####
treebat_summary <- treebat %>%
  filter(session=="3"|session=="4"|session=="5") %>% ##October-December counts only
  filter(is.na(BFF.Fw)==FALSE) %>%
  filter(is.na(BFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","BFF.Fw", "BFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.Fw = sum(BFF.Fw),
        total.F = sum(BFF.F),
        total.Fw.F.ratio = total.Fw/total.F) %>%
  mutate(urban = ifelse(site.code=="DCAN"|site.code=="DCLU", "Non-urban", ## Roosts >10km from the urban edge
                        ifelse(site.code=="DLIS", "Neither", ## Lismore roost was only 0.4km from the urban edge
                               "Urban")))  %>% 
  filter(is.na(total.Fw.F.ratio)==FALSE) %>%
  #View(treebat_summary[,c("urban","site.code", "session", "total.Fw", "total.F", "total.Fw.F.ratio")])
  ddply(c("urban"), summarise,
        ## Calculate average ratio across all urban and non-urban roost types
        mean.Fw.F.ratio = mean(total.Fw.F.ratio))
## For every 1 BFF female there is 0.51 BFF females with pups (urban) and 0.33 BFF females with pups (non-urban)

#### BFF, by contemporary/non-contemporary roosts ####
treebat_summary <- treebat %>%
  filter(session=="3"|session=="4"|session=="5") %>% ##October-December counts only
  filter(is.na(BFF.Fw)==FALSE) %>%
  filter(is.na(BFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","BFF.Fw", "BFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.Fw = sum(BFF.Fw),
        total.F = sum(BFF.F),
        total.Fw.F.ratio = total.Fw/total.F) %>%
  mutate(contemp = ifelse(site.code=="DSUN"|site.code=="DAVO"|site.code=="DBUR", "contemporary", ## see Table 1
                          "non-contemporary"))  %>% 
  filter(is.na(total.Fw.F.ratio)==FALSE) %>%
  #View(treebat_summary[,c("contemp","site.code", "session", "total.Fw", "total.F", "total.Fw.F.ratio")])
  ddply(c("contemp"), summarise,
        ## Calculate average ratio across all contemporary and non-contemporary roosts
        mean.Fw.F.ratio = mean(total.Fw.F.ratio))
## For every 1 BFF female there is 0.51 BFF females with pups (contemporary) and 0.40 BFF females with pups (non-contemporary)


#### GHFF, by urban/non-urban area, October-December ####
treebat_summary <- treebat %>%
  filter(session=="3"|session=="4"|session=="5") %>% ##October-December counts only
  filter(is.na(GHFF.Fw)==FALSE) %>%
  filter(is.na(GHFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","GHFF.Fw", "GHFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.Fw = sum(GHFF.Fw),
        total.F = sum(GHFF.F),
        total.Fw.F.ratio = total.Fw/total.F) %>%
  mutate(urban = ifelse(site.code=="DCAN"|site.code=="DCLU", "Non-urban", ## Roosts >10km from the urban edge
                        ifelse(site.code=="DLIS", "Neither", ## Lismore roost was only 0.4km from the urban edge
                               "Urban")))  %>% 
  filter(is.na(total.Fw.F.ratio)==FALSE) %>%
  filter(!total.Fw.F.ratio==Inf) %>%
  #View(treebat_summary[,c("urban","site.code", "session", "total.Fw", "total.F", "total.Fw.F.ratio")])
  ddply(c("urban"), summarise,
        ## Calculate average ratio across all urban and non-urban roost types
        mean.Fw.F.ratio = mean(total.Fw.F.ratio))
## For every 1 GHFF female there is 0.41 GHFF females with pups (urban) and 0.46 GHFF females with pups (non-urban)


#### GHFF, by contemporary/non-contemporary roosts, October-December ####
treebat_summary <- treebat %>%
  filter(session=="3"|session=="4"|session=="5") %>% ##October-December counts only
  filter(is.na(GHFF.Fw)==FALSE) %>%
  filter(is.na(GHFF.F)==FALSE) %>%
  #View(treebat_summary[,c("site.accession","GHFF.Fw", "GHFF.F")])
  ddply(c("site.code", "session"), summarise,
        ## Sum all males and all females across each roost per session, and calculate M:F ratio for each site*session
        total.Fw = sum(GHFF.Fw),
        total.F = sum(GHFF.F),
        total.Fw.F.ratio = total.Fw/total.F) %>%
  mutate(contemp = ifelse(site.code=="DSUN"|site.code=="DAVO"|site.code=="DBUR", "contemporary", ## see Table 1
                          "non-contemporary"))  %>% 
  filter(is.na(total.Fw.F.ratio)==FALSE) %>%
  filter(!total.Fw.F.ratio==Inf) %>%
  #View(treebat_summary[,c("contemp","site.code", "session", "total.Fw", "total.F", "total.Fw.F.ratio")])
  ddply(c("contemp"), summarise,
        ## Calculate average ratio across all contemporary and non-contemporary roosts
        mean.Fw.F.ratio = mean(total.Fw.F.ratio))
## For every 1 GHFF female there is 0.47 GHFF males (contemporary) and 0.43 GHFF males (non-contemporary)

#####################################################
#### Average maximum roosting height per species ####
#####################################################

heights.subset %>%
  filter(!is.na(max)) %>%
  ddply(c("species"), summarise,
        max_height = mean(max),
        LIQR_max = Liqr(max),
        UIQR_max = Uiqr(max),
        min_height = mean(min),
        LIQR_min = Liqr(min),
        UIQR_min = Uiqr(min)
  )

########################################
#### Total number of trees occupied ####
########################################

tree_status <- treebat %>%
  filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>%
  ddply(c("site.code", "site.accession", "tree.accession"), summarise,
        occ = ifelse(BFF.index>0 | GHFF.index>0 | LRFF.index > 0, 1,0),
        non_occ = ifelse(BFF.index==0 & GHFF.index==0 & LRFF.index == 0, 1,0))
sum(tree_status$occ) #total number of trees with bats throughout the survey
nrow(tree_status) #total trees surveyed (should be somewhere around 32,786 = 2,522*13. But will be lower because of tree removal)

#################################
#### Co-occurrence in roosts ####
#################################

## Number of roosts (per survey) with more than one species:
multispecies_session <- treebat %>%
  filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>%
  ddply(c("site.code", "site.accession"), summarise,
        BFF = sum(BFF.index), 
        GHFF = sum(GHFF.index), 
        LRFF = sum(LRFF.index),
        BG = count(BFF>0&GHFF>0), #number of times they occurred in the same survey
        BR = count(BFF>0&LRFF>0), #number of times they occurred in the same survey
        GR = count(GHFF>0&LRFF>0), #number of times they occurred in the same survey
        #BFF_bin = ifelse(BFF==0,0,1),
        #GHFF_bin = ifelse(GHFF==0,0,1),
        #LRFF_bin = ifelse(LRFF==0,0,1),
        any_bin = ifelse(sum(BG, BR, GR)==0,0,1)
  ) 

## Above, summarized:
sum(multispecies_session$any_bin) #number of bat surveys where more than one species was present
sum(multispecies_session$any_bin)/nrow(multispecies_session) #proportion of bat surveys where more than one species was present

sum(multispecies_session$BG) #number of bat surveys where BFF and GHFF co-occurred
sum(multispecies_session$BG)/nrow(multispecies_session) #proportion of bat surveys where BFF and GHFF co-occurred

sum(multispecies_session$BR) #number of bat surveys where BFF and LRFF co-occurred
sum(multispecies_session$BR)/nrow(multispecies_session) #proportion of bat surveys where BFF and LRFF co-occurred

sum(multispecies_session$GR) #number of bat surveys where GHFF and LRFF co-occurred
sum(multispecies_session$GR)/nrow(multispecies_session) #proportion of bat surveys where GHFF and LRFF co-occurred

##--------------------------------------------------------------------------------------------##

##----------------------------------------- Figures ------------------------------------------##
#### Set labels for facets: ####

site.labs <- c("Avondale (C)", "Burleigh (C)", "Canungra", "Clunes", "Lismore", "Redcliffe", "Sunnybank (C)", "Toowoomba")
names(site.labs) <- c("DAVO", "DBUR", "DCAN", "DCLU", "DLIS", "DRED", "DSUN", "DTOW")

spp.labs <- c("All", "Black", "Grey-headed", "Little red")
names(spp.labs) <- c("all", "BFF", "GHFF", "LRFF")

sp.labs <- c("Black", "Grey-headed", "Little red")
names(sp.labs) <- c("BFF", "GHFF", "LRFF")

##ordered by site, and adjusted for species (i.e. site removed if species absent)
##All/BFF:
site.labs.ordered <- c("Avondale (C)", "Sunnybank (C)", "Burleigh (C)", "Redcliffe", "Toowoomba","Canungra", "Clunes", "Lismore")
names(site.labs.ordered) <- c("01-DAVO", "02-DSUN", "03-DBUR", "04-DRED", "05-DTOW", "06-DCAN", "07-DCLU", "08-DLIS")

##GHFF:
site.labs.ordered_GHFF <- c("Avondale (C)", "Sunnybank (C)", "Redcliffe", "Toowoomba","Canungra", "Clunes", "Lismore")
names(site.labs.ordered_GHFF) <- c("01-DAVO", "02-DSUN", "04-DRED", "05-DTOW", "06-DCAN", "07-DCLU", "08-DLIS")

##LRFF:
site.labs.ordered_LRFF <- c("Avondale (C)", "Sunnybank (C)", "Redcliffe", "Toowoomba")
names(site.labs.ordered_LRFF) <- c("01-DAVO", "02-DSUN", "04-DRED", "05-DTOW")

## Cheat/lazy facets
cheat.labs <- c("All roosts")
names(cheat.labs) <- c("all")


####------------------------------------------------------------------------------------------##
##----------------------------------------- Figure 2 -----------------------------------------##
##------------------------------ Subplot occupancy over time ---------------------------------##
####------------------------------------------------------------------------------------------##

##identify sessions there were ANY bats
bysite <- c("site.code", "session", "site.accession")
occ.all.output_site <- occ.all(treebat, bysite)
#occ.BFF.output_site <- occ.BFF(treebat, bysite)
#occ.GHFF.output_site <- occ.GHFF(treebat, bysite)
#occ.LRFF.output_site <- occ.LRFF(treebat, bysite)

##############################################
###### Total bats in subplot (boxplot) #######
##############################################

### Species combined
x <- "subplot"
y <- "N"
xlab <- "Subplot"
ylab <- "Total number of bats per subplot"
title <- "Subplot occupancy (total per subplot, approximated with weighted index value) when bats present at roost"
flab <- "Site"
clab <- "Site"
colour <- "Site"
data <- index.TREE.wide_plot %>%
  filter(species == "all") %>% #extract combined species measure only
  filter(site.accession %in% occ.all.output_site$site.accession) %>% #choose surveys when at least 1 bat was present, only
  mutate(plot.occ = ifelse(occ==0, 0,
                           ifelse(occ>0, 1,
                                  NA))) %>%
  mutate(subplot = as.numeric(as.character(subplot))) %>%
  create.Site()

## Facet by site, colour by site
p1_occup_site_fac_BOX <-  ggplot(data) + 
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
  geom_boxplot(aes(x=as.factor(.data[[x]]), y = .data[[y]], fill=Site)) +
  scale_fill_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore")) +
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered)) + # .~ facet into collumns and # ~. facet into rows
  scale_x_discrete(breaks=c(2,4,6,8, 10), labels=c("2","4","6","8", "10"))

##############################################
### Proportion of times occupied (barplot) ###
##############################################

### Species combined
x <- "subplot"
y <- "none"
xlab <- "Subplot"
ylab <- "Proportion of surveys subplot was occupied"
title <- "Subplot occupancy (proportion of times subplot was occupied) when bats present at roost"
flab <- "Site"
clab <- "Site"
colour <- "Site"
data <- index.TREE.wide_plot %>%
  filter(species == "all") %>% #extract combined species measure only
  filter(site.accession %in% occ.all.output_site$site.accession) %>% #choose surveys when at least 1 bat was present, only
  mutate(plot.occ = ifelse(occ==0, 0,
                           ifelse(occ>0, 1,
                                  NA))) %>%
  ddply(c("site.code", "subplot"), summarise,
        prop.occ = sum(plot.occ)/sum(!is.na(plot.occ))) %>%
  #mutate(subplot = as.numeric(as.character(subplot))) %>%
  create.Site()

## Facet by site, colour by site
p1_occup_site_fac_BAR <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  #geom_point(aes(x=data[[x]], y = data[[y]], colour=Site), size=1) +
  #scale_colour_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore")) +
  geom_bar(aes(x= .data[[x]], y = prop.occ, fill=Site), stat="identity") +
  scale_fill_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore")) +
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered)) +  # .~ facet into collumns and # ~. facet into rows
  scale_x_discrete(breaks=c(2,4,6,8, 10), labels=c("2","4","6","8", "10"))

### Arrange and save plots:
jpeg("Output/Figures/Figure 2.png", width=1400, height=1200)
ggarrange(p1_occup_site_fac_BOX+coord_cartesian(ylim=c(0, 750))+ggtitle("")+ theme(legend.position="none"), p1_occup_site_fac_BAR+ggtitle("")+ theme(legend.position="none"), nrow = 2, heights = c(3,2), labels = c("A", "B"), font.label = list(size = 24)) #26_plot occupancy when bats are present at site - total per subplot and proportion of times occupied_ZOOMED_1400 x 800
dev.off()


####------------------------------------------------------------------------------------------##
##----------------------------------------- Figure 3 -----------------------------------------##
##-------------------------- Density in core vs peripheral plots -----------------------------##
####------------------------------------------------------------------------------------------##

##-------------------------------- Produce histogram of values -------------------------------##
## Histogram of count of subplots occupied, when roosts occupied
## To derive threshold value of occupancy to define core and peripheral plot occupancy

##identify sessions there were ANY bats:
bysite <- c("site.code", "session", "site.accession")
occ.all.output_site <- occ.all(treebat, bysite)
#occ.BFF.output_site <- occ.BFF(treebat, bysite)
#occ.GHFF.output_site <- occ.GHFF(treebat, bysite)
#occ.LRFF.output_site <- occ.LRFF(treebat, bysite)

## Histogram of values:
p1_occup_site_fac_HIST <- index.TREE.wide_plot %>%
  filter(species == "all") %>% #extract combined species measure only
  filter(site.accession %in% occ.all.output_site$site.accession) %>% #choose surveys when at least 1 bat was present, only
  mutate(plot.occ = ifelse(occ==0, 0,
                           ifelse(occ>0, 1,
                                  NA))) %>%
  ddply(c("site.code", "subplot"), summarise,
        prop.occ = sum(plot.occ)/sum(!is.na(plot.occ))) %>%
  ggplot(aes(x=prop.occ)) + 
  geom_histogram(binwidth=0.05) +
  theme_bw() +
  background_grid("none")+
  labs(y="Count of subplots", x="Proportion of surveys subplot was occupied (occupancy when bats present at site)")+
  ggtitle("Histogram showing frequency of occupancy levels") +
  stat_bin(binwidth=0.05, geom="text", colour="white", size=3.5,
           aes(label=..count..), position=position_stack(vjust=0.5)) +
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  geom_vline(xintercept = 0.8, linetype="dashed", 
             color = "red", size=1) #27_frequency distribution of prop subplot occupancy_600 x 600

jpeg("Output/Figures/SI_frequency distribution of prop subplot occupancy.png", width=800, height=600)
plot(p1_occup_site_fac_HIST)
dev.off()

##-------------------------------- Make the comparison plots ---------------------------------##

### Identify subplots as core or peripheral:
threshold <- 0.8 #set threshold (from above)
core.plots <- index.TREE.wide_plot %>%
  filter(species == "all") %>% #extract combined species measure only
  filter(site.accession %in% occ.all.output_site$site.accession) %>% #choose surveys when at least 1 bat was present, only
  mutate(plot.occ = ifelse(occ==0, 0,
                           ifelse(occ>0, 1,
                                  NA))) %>%
  ddply(c("site.code", "subplot"), summarise,
        prop.occ = sum(plot.occ)/sum(!is.na(plot.occ))) %>% #calculate proportion of times each subplot was occupied. Use this to assign core vs peripheral rating
  mutate(occupancy.cat = ifelse(prop.occ>=threshold,"Core","Peripheral")) #set threshold for "core" occupancy

## Identify subplots there were ANY bats
byplot <- c("site.code", "session", "site.accession", "subplot", "rep")
occ.all.output_plot <- occ.all(treebat, byplot)
#occ.BFF.output_plot <- occ.BFF(treebat, byplot)
#occ.GHFF.output_plot <- occ.GHFF(treebat, byplot)
#occ.LRFF.output_plot <- occ.LRFF(treebat, byplot)

##############################################
#### VS total number of bats per subplot #####
##############################################
### Species combined
x <- "occupancy.cat"
y <- "N"
xlab <- "Subplot category"
ylab <- "Total number of bats per subplot"
title <- "Subplot occupancy (total per subplot, approximated with weighted index value) when bats present in subplot"
data <- index.TREE.wide_plot %>%
  full_join(core.plots, by = c("site.code", "subplot")) %>% #join with occupancy category
  filter(species == "all") %>% #extract combined species measure only
  filter(rep %in% occ.all.output_plot$rep) %>% #choose *plots* when at least 1 bat was present, only
  mutate(N = as.numeric(as.character(N))) %>%
  create.Site()

## Facet by site, colour by site
flab <- "Site"
clab <- "Site"
colour <- "Site"
p1_occup_site_fac_BOX <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  #geom_point(aes(x=.data[[x]], y = .data[[y]], colour=Site), size=1) +
  #scale_colour_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore")) +
  geom_boxplot(aes(x=.data[[x]], y = .data[[y]], fill=Site)) +
  scale_fill_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore"))+
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered)) # .~ facet into collumns and # ~. facet into rows

## No facet, colour by occupancy category
flab <- "Subplot category"
clab <- "Subplot category"
colour <- "Subplot category"
p1_occup_site_all_BOX <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  #geom_point(aes(x=.data[[x]], y = .data[[y]]), size=1) +
  geom_boxplot(aes(x=.data[[x]], y = .data[[y]], fill=occupancy.cat)) +
  scale_fill_manual(values=c("gray47", "gray87"), labels = c("Core", "Peripheral")) + 
  facet_grid(.~species, labeller = labeller(species = cheat.labs)) 

### Arrange and save plots:
#ggarrange(p1_occup_site_fac_BOX+theme(legend.position="none"), #28_Comparison of core & peri occupancy, with total number of bats per subplot_1000 x 800
#          p1_occup_site_all_BOX+ggtitle("") + theme(legend.position="none", axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()), 
#          nrow = 1, widths = c(3, 1), labels = c("A", "B") ) 
jpeg("Output/Figures/Figure 3.png", width=1400, height=800)
ggarrange(p1_occup_site_fac_BOX+theme(legend.position="none")+coord_cartesian(ylim=c(0, 750)), #28_Comparison of core & peri occupancy, with total number of bats per subplot_ZOOMED_1000 x 800
          p1_occup_site_all_BOX+coord_cartesian(ylim=c(0, 750))+ggtitle("") + theme(legend.position="none", axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()), 
          nrow = 1, widths = c(3, 0.35), labels = c("A", "B"), font.label = list(size = 24))
dev.off()


###########################################
##### VS proportion of occupied trees #####
###########################################
x <- "occupancy.cat"
y <- "none"
xlab <- "Subplot category"
ylab <- "Occupied trees (%)"
title <- "Subplot occupancy (% occupied trees per subplot)  when bats present in subplot"
data <- index.TREE.wide_plot %>%
  full_join(core.plots, by = c("site.code", "subplot")) %>% #join with occupancy category
  filter(species == "all") %>% #extract combined species measure only
  filter(occ >0) %>% #filter unoccupied subplots
  mutate(occ = as.numeric(as.character(occ))) %>%
  create.Site()

## Facet by site, colour by site
flab <- "Site"
clab <- "Site"
colour <- "Site"
p4_occup_site_fac_BOX <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  #geom_point(aes(x=.data[[x]], y = (occ/tree.count)*100, colour=Site), size=1) +
  #scale_colour_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore")) +
  geom_boxplot(aes(x=.data[[x]], y = (occ/tree.count)*100, fill=Site)) +
  scale_fill_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore"))+
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered)) # .~ facet into collumns and # ~. facet into rows

## No facet, colour by occupancy category
flab <- "Subplot category"
clab <- "Subplot category"
colour <- "Subplot category"
p4_occup_site_all_BOX <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  #geom_point(aes(x=.data[[x]], y = (occ/tree.count)*100, size=1) +
  geom_boxplot(aes(x=.data[[x]], y = (occ/tree.count)*100, fill=occupancy.cat)) +
  scale_fill_manual(values=c("gray47", "gray87"), labels = c("Core", "Peripheral")) + 
  facet_grid(.~species, labeller = labeller(species = cheat.labs)) 

### Arrange and save plots:
jpeg("Output/Figures/SI_prop occ_comparison of core & peri occupancy, with the prop of occupied trees in occupied subplots.png", width=1400, height=800)
ggarrange(p4_occup_site_fac_BOX+theme(legend.position="none"), #31_Comparison of core & peri occupancy, with the prop of occupied trees in occupied subplots_1000 x 800
          p4_occup_site_all_BOX+ggtitle("") + theme(legend.position="none", axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()), 
          nrow = 1, widths = c(3, 0.35), labels = c("A", "B"), font.label = list(size = 24))
dev.off()


####------------------------------------------------------------------------------------------##
##----------------------------------------- Supp Info ----------------------------------------##
##-------------------------------- Core vs peripheral trees -----------------------------------##
####------------------------------------------------------------------------------------------##

##-------------------------------- Produce histogram of values -------------------------------##
## Histogram of count of trees occupied, when roosts occupied
## To derive threshold value of occupancy to define core and peripheral tree occupancy

##identify sessions there were ANY bats
bysite <- c("site.code", "session", "site.accession")
occ.all.output_site <- occ.all(treebat, bysite)
#occ.BFF.output_site <- occ.BFF(treebat, bysite)
#occ.GHFF.output_site <- occ.GHFF(treebat, bysite)
#occ.LRFF.output_site <- occ.LRFF(treebat, bysite)

## Histogram of values (frequency of % occupied, per tree):
## Sites combined:
p1_occup_tree_site_all_HIST <- treebat %>%
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed. This is important to do so that proportion is calculated correctly
  mutate(all.index.weight = BFF.index.weight+GHFF.index.weight+LRFF.index.weight) %>%
  filter(site.accession %in% occ.all.output_site$site.accession) %>% #choose surveys when at least 1 bat was present, only (I think this is more appropriate for identifying core/peri trees, over occupied subplots
  mutate(tree.occ = ifelse(all.index.weight==0, 0, #Assign binary value to show if tree is occupied
                           ifelse(all.index.weight>0, 1,
                                  NA))) %>% #nrow of tree.occ should be 2522 trees (this is how many were tagged in total)
  ddply(c("site.code", "tree.accession"), summarise, #for each tree calculate the number of times it was occupied, across surveys where there was one bat present at roost
        tree.occ = sum(tree.occ)/sum(!is.na(tree.occ))) %>%
  create.Site() %>%
  ggplot(aes(x=tree.occ)) + 
  geom_histogram(binwidth=0.1) +
  theme_bw() +
  background_grid("none")+
  labs(y="Count of trees", x="Proportion of surveys tree was occupied (occupancy when bats present at site)")+
  ggtitle("Histogram showing frequency of tree occupancy") +
  stat_bin(binwidth=0.1, geom="text", colour="white", size=3.5,
           aes(label=..count..), position=position_stack(vjust=0.5)) +
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  geom_vline(xintercept = 0.8, linetype="dashed", 
             color = "red", size=1)

## Check value range of bins:
str(p1_occup_tree_site_all_HIST) #2482 obs. of  3 variables - so 2482 trees plotted in total. This should be correct (2522 total tagged - 40 felled in Canungra prior to first survey)
check_hist <- treebat %>%
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed. This is important to do so that proportion is calculated correctly
  mutate(all.index.weight = BFF.index.weight+GHFF.index.weight+LRFF.index.weight) %>%
  filter(site.accession %in% occ.all.output_site$site.accession) %>% #choose surveys when at least 1 bat was present, only (I think this is more appropriate for identifying core/peri trees, over occupied subplots
  mutate(tree.occ = ifelse(all.index.weight==0, 0, #Assign binary value to show if tree is occupied
                           ifelse(all.index.weight>0, 1,
                                  NA))) %>% #nrow of tree.occ should be 2522 trees (this is how many were tagged in total)
  ddply(c("site.code", "tree.accession"), summarise, #for each tree calculate the number of times it was occupied, across surveys where there was one bat present at roost
        tree.occ = sum(tree.occ)/sum(!is.na(tree.occ))) %>%
  create.Site()
nrow(check_hist[check_hist$tree.occ==1.00,]) #70
nrow(check_hist[check_hist$tree.occ<1 & check_hist$tree.occ>=0.90,]) 
nrow(check_hist[check_hist$tree.occ<0.90 & check_hist$tree.occ>=0.80,])
nrow(check_hist[check_hist$tree.occ<0.80 & check_hist$tree.occ>=0.75,]) 
nrow(check_hist[check_hist$tree.occ<0.70 & check_hist$tree.occ>=0.65,])
nrow(check_hist[check_hist$tree.occ<=0.3 & check_hist$tree.occ>0.15,]) #
nrow(check_hist[check_hist$tree.occ<=0.15 & check_hist$tree.occ>0,]) #295
nrow(check_hist[check_hist$tree.occ==0,]) #794

## Split by site:
p1_occup_tree_site_fac_HIST <- treebat %>%
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed. This is important to do so that proportion is calculated correctly
  mutate(all.index.weight = BFF.index.weight+GHFF.index.weight+LRFF.index.weight) %>%
  filter(site.accession %in% occ.all.output_site$site.accession) %>% #choose surveys when at least 1 bat was present, only (I think this is more appropriate for identifying core/peri trees, over occupied subplots
  mutate(tree.occ = ifelse(all.index.weight==0, 0, #Assign binary value to show if tree is occupied
                           ifelse(all.index.weight>0, 1,
                                  NA))) %>% #nrow of tree.occ should be 2522 trees (this is how many were tagged in total)
  ddply(c("site.code", "tree.accession"), summarise, #for each tree calculate the number of times it was occupied, across surveys where there was one bat present at roost
        tree.occ = sum(tree.occ)/sum(!is.na(tree.occ))) %>%
  create.Site() %>%
  ggplot(aes(x=tree.occ)) + 
  geom_histogram(binwidth=0.1, aes(fill=Site)) +
  theme_bw() +
  background_grid("none")+
  scale_fill_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore"))+
  labs(y="Count of trees", x="Proportion of surveys tree was occupied (occupancy when bats present at site)", fill="Site")+
  ggtitle("Histogram showing frequency of tree occupancy") +
  stat_bin(binwidth=0.1, geom="text", colour="white", size=3, aes(label=..count..), position=position_stack(vjust=0.5)) +
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  geom_vline(xintercept = 0.8, linetype="dashed", 
             color = "red", size=1) +
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered))  # .~ facet into collumns and # ~. facet into rows
#scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), labels=c("0", ">0-0.1", ">0.1-0.2", ">0.2-0.3", ">0.3-0.4", ">0.4-0.5", ">0.5-0.6", ">0.6-0.7", ">0.7-0.8", ">0.8-0.9", ">0.9-1.0")) 
# Note - can't seem to remove the labels without removing the colour of the bars
# Note - bins aren't 0.1, but not sure how to determine what they really are

### Arrange and save plots:
jpeg("Output/Figures/SI_frequency distribution of prop tree occupancy.png", width=1400, height=1000)
ggarrange(p1_occup_tree_site_fac_HIST+ theme(legend.position="none"), p1_occup_tree_site_all_HIST+ggtitle("")+ theme(legend.position="none"), nrow = 2, heights = c(3,2), labels = c("A", "B"), font.label = list(size = 24)) #38_frequency distribution of prop tree occupancy_1400 x 800
dev.off()


####------------------------------------------------------------------------------------------##
##----------------------------------------- Figure 4 -----------------------------------------##
##---------------------------- Occupancy vs distance from center -----------------------------##
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

##############################################
#### VS total number of bats per subplot #####
##############################################
### All species combined ###
x <- "scaled_eud.distance_mean"
y <- "N"
ylab <- "Total number of bats per subplot"
xlab <- "Mean scaled tree distance from roost center (scaled by maximum distance value) per subplot"
flab <- "Site"
title <- "Total bats per subplot, approximated with weighted index value"

## Calculate mean scaled distance (note - mean of scaled distances, not mean distance then scaled)
scaled_treebat <- treebat %>% 
  ddply(c("site.code", "session","site.accession", "subplot", "rep"), summarise,
        scaled_eud.distance_mean = mean(scaled_eud.distance)) #average per tree across subplot
index.TREE.wide_plot <- full_join(scaled_treebat, index.TREE.wide_plot, by = c("site.code", "session", "site.accession", "subplot")) 

### Species combined ###
clab <- "Site"
colour <- "Site"
data <- index.TREE.wide_plot %>% #Note N is estimated from weighted index value per tree, not true count per tree
  filter(species == "all")  %>% #extract combined species measure only
  #filter(site.accession %in% occ.all.output_plot$site.accession) %>% #don't want to filter unnocupied subplots in this case
  create.Site()

## No facet, no colour
p5_centroid_site_all <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  stat_smooth(aes(y = .data[[y]]), method = "loess", se = TRUE, size=1, colour="black") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) 

### Split by species ###
## Facet by site, colour by species 
dash <- "species"
dlab <- "Species"
flab <- "Species"
data_BFF <- index.TREE.wide_plot %>% #Note that this is estimated from weighted index value per tree, not true count per tree
  filter(species == "BFF")  %>% 
  #filter(site.accession %in% occ.BFF.output_plot$site.accession) %>% #don't want to filter unnocupied subplots in this case
  create.Site()
data_GHFF <- index.TREE.wide_plot %>% #Note that this is estimated from weighted index value per tree, not true count per tree
  filter(species == "GHFF")  %>% 
  #filter(site.accession %in% occ.GHFF.output_plot$site.accession) %>% #don't want to filter unnocupied subplots in this case
  create.Site()
data_LRFF <- index.TREE.wide_plot %>% #Note that this is estimated from weighted index value per tree, not true count per tree
  filter(species == "LRFF")  %>%
  #filter(site.accession %in% occ.LRFF.output_plot$site.accession) %>% #don't want to filter unnocupied subplots in this case
  create.Site()
combined_data <- rbind(data_BFF, data_GHFF, data_LRFF)

p5_centroid_site_spp <- plot_dash_smooth(combined_data, x, y, dash, ylab, xlab, title, dlab, flab) +
  scale_linetype_manual(values=c("solid", "longdash", "dotted"), labels = c("Black", "Grey-headed", "Little red"))+
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered), scales="free_x") 

### Arrange and save plots:
jpeg("Output/Figures/Figure 4.png", width=1400, height=1200)
ggarrange(p5_centroid_site_spp+coord_cartesian(ylim=c(0, 750))+ggtitle(""), p5_centroid_site_all+coord_cartesian(ylim=c(0, 750))+ggtitle(""), nrow = 2, heights = c(3,1.5), labels = c("A", "B"), font.label = list(size = 24), common.legend = TRUE, legend = "bottom") #44-02_mean scaled distance from roost center vs total bats per subplot_ZOOMED_1400 x 1200
dev.off()

###########################################
##### VS proportion of occupied trees #####
###########################################
### All species combined ###
x <- "scaled_eud.distance_mean"
y <- "none"
ylab <- "Occupied trees (%)"
xlab <- "Mean scaled tree distance from roost center (scaled by maximum distance value) per subplot"
flab <- "Site"
title <- "Proportion of occupied trees (% per subplot, un-occupied subplots removed)"

## Calculate mean scaled distance (note - mean of scaled distances, not mean distance then scaled)
scaled_treebat <- treebat %>% 
  ddply(c("site.code", "session","site.accession", "subplot", "rep"), summarise,
        scaled_eud.distance_mean = mean(scaled_eud.distance)) #average per tree across subplot
index.TREE.wide_plot <- full_join(scaled_treebat, index.TREE.wide_plot, by = c("site.code", "session", "site.accession", "subplot")) 

### Species combined ###
clab <- "Site"
colour <- "Site"
data <- index.TREE.wide_plot %>% #Note N is estimated from weighted index value per tree, not true count per tree
  filter(species == "all")  %>% #extract combined species measure only
  #filter(site.accession %in% occ.all.output_plot$site.accession) %>% #don't want to filter unnocupied subplots in this case
  create.Site()

## No facet, no colour 
p4_centroid_site_all <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  stat_smooth(aes(y = (occ/tree.count)*100), method = "loess", se = TRUE, size=1, colour="black")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) 

## Facet by site, colour by species 
colour <- "species"
clab <- "Species"
flab <- "Species"
data_BFF <- index.TREE.wide_plot %>% #Note that this is estimated from weighted index value per tree, not true count per tree
  filter(species == "BFF")  %>% 
  #filter(site.accession %in% occ.BFF.output_plot$site.accession) %>% #don't want to filter unnocupied subplots in this case
  create.Site()
data_GHFF <- index.TREE.wide_plot %>% #Note that this is estimated from weighted index value per tree, not true count per tree
  filter(species == "GHFF")  %>% 
  #filter(site.accession %in% occ.GHFF.output_plot$site.accession) %>% #don't want to filter unnocupied subplots in this case
  create.Site()
data_LRFF <- index.TREE.wide_plot %>% #Note that this is estimated from weighted index value per tree, not true count per tree
  filter(species == "LRFF")  %>%
  #filter(site.accession %in% occ.LRFF.output_plot$site.accession) %>% #don't want to filter unnocupied subplots in this case
  create.Site()
combined_data <- rbind(data_BFF, data_GHFF, data_LRFF)

p4_centroid_site_spp <- plot_blank(combined_data, x, y, ylab, xlab, title, clab, flab) +
  stat_smooth(aes(y = (occ/tree.count)*100, linetype=.data[[colour]]), method = "loess", se = TRUE, size=1, colour="black") +
  scale_linetype_manual(values=c("solid", "longdash", "dotted"), labels = c("Black", "Grey-headed", "Little red"))+
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered), scales="free_x") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) 

### Arrange and save plots:
jpeg("Output/Figures/SI_roost center vs prop trees occupied.png", width=1400, height=1200)
ggarrange(p4_centroid_site_spp+coord_cartesian(ylim=c(0, 100))+ggtitle(""), p4_centroid_site_all+ggtitle("")+coord_cartesian(ylim=c(0, 100)), nrow = 2, heights = c(3,1.5), labels = c("A", "B"), font.label = list(size = 24)) #45-02_mean scaled distance from roost center vs prop trees occupied_FIXED_1400 x 1200
dev.off()

####------------------------------------------------------------------------------------------##
##----------------------------------------- Figure 5 -----------------------------------------##
##------------------------------ Roost abundance vs roost area -------------------------------##
####------------------------------------------------------------------------------------------##

## Note - points are sites per session
data <- roostbat %>% filter(!is.na(site.code)) %>% create.Site()
x <- "index.abundance"
y <- "roost.area"
ylab <- "Total roost area"
xlab <- "Total roost abundance (index score)"
title <- "Global abundance vs roost area"
flab <- "none"
clab <- "Site"
colour <- "Site"

## Facet by site, colour by site
p4_GL_site_fac <- plot_black_smooth(data, x, y, colour, ylab, xlab, title, clab, flab) +
  geom_point(aes(y=.data[[y]], color=.data[[colour]]), size=4) +
  scale_color_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore"))+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7), labels=c("0","1-499","500-2,499","2,500 - 4,999","5,000 - 9,999","10,000 - 15,999","16,000 - 49,999","> 50,000")) +
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered), scales="free_x") 

## No facet, colour by site
p4_GL_site_all <- plot_black_smooth(data, x, y, colour, ylab, xlab, title, clab, flab) +
  geom_point(aes(y=.data[[y]], color=.data[[colour]]), size=4) +
  scale_color_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore"))+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7), labels=c("0","1-499","500-2,499","2,500 - 4,999","5,000 - 9,999","10,000 - 15,999","16,000 - 49,999","> 50,000")) 

## Arrange and save plots:
jpeg("Output/Figures/Figure 5.png", width=1400, height=1200)
ggarrange(p4_GL_site_fac+ggtitle("")+coord_cartesian(ylim=c(0, 90000)), p4_GL_site_all+coord_cartesian(ylim=c(0, 90000))+ggtitle(""), nrow = 2, labels = c("A", "B"),common.legend = TRUE, legend = "bottom", font.label = list(size = 24)) #12_Roost index abundance vs roost area_ZOOMED_1400 x 800
dev.off()

####------------------------------------------------------------------------------------------##
##----------------------------------------- Figure 6 -----------------------------------------##
##-------------------------------- Bar plots of co-occupation --------------------------------##
####------------------------------------------------------------------------------------------##

#########################################
#### Co-occurrence in occupied trees ####
#########################################

##identify session there were specific bat pairings:
bysite <- c("site.code", "session", "site.accession")
occ.BFF.output_site <- occ.BFF(treebat, bysite) 
occ.GHFF.output_site <- occ.GHFF(treebat, bysite) 
occ.LRFF.output_site <- occ.LRFF(treebat, bysite) 

## From sessions where both BFF and GHFF were present, how many trees did they co-occupy?
BFF_GHFF_tree <- treebat %>%
  filter(site.accession %in% occ.BFF.output_site$site.accession & site.accession %in% occ.GHFF.output_site$site.accession) %>% #identify sessions where the species pair was present
  filter(!is.na(BFF.index)&!is.na(GHFF.index)) %>%
  ddply(c("site.code", "site.accession"), summarise,
        BG = count(BFF.index>0&GHFF.index>0), #number of times they occurred in the same tree
        BG_occ = count(BFF.index>0|GHFF.index>0), #number of trees with either BFF or GHFF
        non_BG = BG_occ - BG
  ) 

BFF_LRFF_tree <- treebat %>%
  filter(site.accession %in% occ.BFF.output_site$site.accession & site.accession %in% occ.LRFF.output_site$site.accession) %>% #identify sessions where the species pair was present
  filter(!is.na(BFF.index)&!is.na(LRFF.index)) %>%
  ddply(c("site.code", "site.accession"), summarise,
        BR = count(BFF.index>0&LRFF.index>0), #number of times they occurred in the same tree
        BR_occ = count(BFF.index>0|LRFF.index>0), #number of trees with either BFF or LRFF
        non_BR = BR_occ - BR
  ) 

GHFF_LRFF_tree <- treebat %>%
  filter(site.accession %in% occ.GHFF.output_site$site.accession & site.accession %in% occ.LRFF.output_site$site.accession) %>% #identify sessions where the species pair was present
  filter(!is.na(GHFF.index)&!is.na(LRFF.index)) %>%
  ddply(c("site.code", "site.accession"), summarise,
        GR = count(GHFF.index>0&LRFF.index>0), #number of times they occurred in the same tree
        GR_occ = count(GHFF.index>0|LRFF.index>0), #number of trees with either GHFF or LRFF
        non_GR = GR_occ - GR
  ) 

temp <- rbind(BFF_GHFF_tree[,c(1:2)], GHFF_LRFF_tree[,c(1:2)], BFF_LRFF_tree[,c(1:2)])
length(unique(temp$site.accession)) ## Check that number of surveys with multiple species matches multispecies_session$any_bin (== 73)

sum(BFF_GHFF_tree$BG) #number of trees where BFF and GHFF co-occurred
sum(BFF_GHFF_tree$BG)/sum(BFF_GHFF_tree$BG_occ) #proportion of trees where BFF and GHFF co-occurred
BFF_GHFF <- binom.confint(x=sum(BFF_GHFF_tree$BG), n=sum(BFF_GHFF_tree$BG_occ), methods = 'wilson')
BFF_GHFF_nonocc <- binom.confint(x=sum(BFF_GHFF_tree$non_BG), n=sum(BFF_GHFF_tree$BG_occ), methods = 'wilson')

sum(BFF_LRFF_tree$BR) #number of trees where BFF and GHFF co-occurred
sum(BFF_LRFF_tree$BR)/sum(BFF_LRFF_tree$BR_occ) #proportion of trees where BFF and GHFF co-occurred
BFF_LRFF <- binom.confint(x=sum(BFF_LRFF_tree$BR), n=sum(BFF_LRFF_tree$BR_occ), methods = 'wilson')
BFF_LRFF_nonocc <- binom.confint(x=sum(BFF_LRFF_tree$non_BR), n=sum(BFF_LRFF_tree$BR_occ), methods = 'wilson')

sum(GHFF_LRFF_tree$GR) #number of trees where BFF and GHFF co-occurred
sum(GHFF_LRFF_tree$GR)/sum(GHFF_LRFF_tree$GR_occ) #proportion of trees where BFF and GHFF co-occurred
GHFF_LRFF <- binom.confint(x=sum(GHFF_LRFF_tree$GR), n=sum(GHFF_LRFF_tree$GR_occ), methods = 'wilson')
GHFF_LRFF_nonocc <- binom.confint(x=sum(GHFF_LRFF_tree$non_GR), n=sum(GHFF_LRFF_tree$GR_occ), methods = 'wilson')

### All species:
## From sessions where all species were present, how many trees did they co-occupy?
BFF_GHFF_LRFF_tree <- treebat %>%
  filter(site.accession %in% occ.BFF.output_site$site.accession & site.accession %in% occ.GHFF.output_site$site.accession & site.accession %in% occ.LRFF.output_site$site.accession) %>% 
  filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>%
  ddply(c("site.code", "site.accession"), summarise,
        BGR = count(BFF.index>0&GHFF.index>0&LRFF.index>0), #number of times they occurred in the same tree
        BGR_occ = count(BFF.index>0|GHFF.index>0|LRFF.index>0), #number of trees any species
        non_BGR = BGR_occ - BGR
  ) 
sum(BFF_GHFF_LRFF_tree$BGR) #number of trees where BFF and GHFF co-occurred
sum(BFF_GHFF_LRFF_tree$BGR)/sum(BFF_GHFF_LRFF_tree$BGR_occ) #proportion of trees where BFF and GHFF co-occurred
BFF_GHFF_LRFF <- binom.confint(x=sum(BFF_GHFF_LRFF_tree$BGR), n=sum(BFF_GHFF_LRFF_tree$BGR_occ), methods = 'wilson')
BFF_GHFF_LRFF_nonocc <- binom.confint(x=sum(BFF_GHFF_LRFF_tree$non_BGR), n=sum(BFF_GHFF_LRFF_tree$BGR_occ), methods = 'wilson')


#### Bar plot of tree co-occupation:
output <- rbind(BFF_GHFF, BFF_GHFF_nonocc, BFF_LRFF, BFF_LRFF_nonocc, GHFF_LRFF, GHFF_LRFF_nonocc, BFF_GHFF_LRFF, BFF_GHFF_LRFF_nonocc)
output$pair <- c("BFF_GHFF", "BFF_GHFF","BFF_LRFF", "BFF_LRFF","GHFF_LRFF","GHFF_LRFF", "X_BFF_GHFF_LRFF", "X_BFF_GHFF_LRFF")
output$occ <- c("Co-occur", "Separate","Co-occur", "Separate","Co-occur", "Separate", "Co-occur", "Separate")
output$label <- paste("Trees = ", output$n)
output[c(2,4,6,8),9] <- ""

pair.lab <- c("BFF & GHFF", "BFF & LRFF", "GHFF & LRFF", "BFF & GHFF & LRFF")
names(pair.lab) <- c("BFF_GHFF","BFF_LRFF","GHFF_LRFF", "X_BFF_GHFF_LRFF")

tree_cooc <- ggplot(output, aes(x=occ)) + 
  geom_bar(aes(x=occ, y = mean, fill=pair), stat="identity") +
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  #scale_x_discrete(breaks=c("BFF_GHFF","BFF_LRFF","GHFF_LRFF"),
  #                 labels=c("BFF & GHFF", "BFF & LRFF", "GHFF & LRFF")) +
  geom_text(aes(label=label), y=1, color="black", size=6)+
  theme_bw() +
  scale_fill_manual(values=c("gray41", "gray63", "gray79", "slategray2"), labels = c("BFF & GHFF", "BFF & LRFF", "GHFF & LRFF", "BFF & GHFF & LRFF")) +
  background_grid(major="x", colour.major = "grey95")+
  labs(y="Proportion of trees occupied", x="Species pair", colour="Species pair", fill="Species pair")+
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  theme(legend.position="none") +
  facet_grid(.~pair, labeller = labeller(pair = pair.lab))

jpeg("Output/Figures/Output/Figures/Figure 6_02.png", width=1000, height=500)
plot(tree_cooc)
dev.off()


#### Co-occurrence in occupied plots ####
##identify plots there were ANY bats:
byplot <- c("site.code", "session", "site.accession", "subplot")
occ.all.output_plot <- occ.all(treebat, byplot) %>%
  mutate(rep = paste(site.accession, subplot, sep="-"))
occ.BFF.output_plot <- occ.BFF(treebat, byplot) %>%
  mutate(rep = paste(site.accession, subplot, sep="-"))
occ.GHFF.output_plot <- occ.GHFF(treebat, byplot) %>%
  mutate(rep = paste(site.accession, subplot, sep="-"))
occ.LRFF.output_plot <- occ.LRFF(treebat, byplot) %>%
  mutate(rep = paste(site.accession, subplot, sep="-"))

multispecies_plot <- treebat %>%
  filter(rep %in% occ.all.output_plot$rep) %>%
  filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>%
  ddply(c("site.code", "site.accession","subplot"), summarise,
        BFF = sum(BFF.index), 
        GHFF = sum(GHFF.index), 
        LRFF = sum(LRFF.index), 
        BG = count(BFF>0&GHFF>0), #number of times they occurred in the same plot
        BR = count(BFF>0&LRFF>0), #number of times they occurred in the same plot
        GR = count(GHFF>0&LRFF>0), #number of times they occurred in the same plot
        #BFF_bin = ifelse(BFF==0,0,1),
        #GHFF_bin = ifelse(GHFF==0,0,1),
        #LRFF_bin = ifelse(LRFF==0,0,1),
        any_bin = ifelse(sum(BG, BR, GR)==0,0,1)
  ) 

sum(multispecies_plot$any_bin) #number of plots where more than one species was present
sum(multispecies_plot$any_bin)/nrow(multispecies_plot) #proportion of plots where more than one species was present
binom.confint(x=sum(multispecies_plot$any_bin), n=nrow(multispecies_plot), methods = 'wilson')


##identify session there were specific bat pairings:
bysite <- c("site.code", "session", "site.accession")
occ.BFF.output_site <- occ.BFF(treebat, bysite) 
occ.GHFF.output_site <- occ.GHFF(treebat, bysite) 
occ.LRFF.output_site <- occ.LRFF(treebat, bysite) 

## Identify sessions where both BFF & GHFF were present:
BFF_GHFF_plot <- treebat %>%
  filter(site.accession %in% occ.BFF.output_plot$site.accession & site.accession %in% occ.GHFF.output_plot$site.accession) %>% #identify sessions where the species pair was present
  filter(rep %in% occ.BFF.output_plot$rep | rep %in% occ.GHFF.output_plot$rep) %>% #identify plots occupied by either species of interest, only
  filter(!is.na(BFF.index)&!is.na(GHFF.index)) %>%
  ddply(c("site.code", "site.accession","subplot"), summarise,
        BFF = sum(BFF.index), 
        GHFF = sum(GHFF.index), 
        BG = count(BFF>0&GHFF>0) #number of times they occurred in the same plot
  ) 

## Identify sessions where both BFF & LRFF were present:
BFF_LRFF_plot <- treebat %>%
  filter(site.accession %in% occ.BFF.output_plot$site.accession & site.accession %in% occ.LRFF.output_plot$site.accession) %>% #identify sessions where the species pair was present
  filter(rep %in% occ.BFF.output_plot$rep | rep %in% occ.LRFF.output_plot$rep) %>% #identify occupied plots only
  filter(!is.na(BFF.index)&!is.na(LRFF.index)) %>%
  ddply(c("site.code", "site.accession","subplot"), summarise,
        BFF = sum(BFF.index), 
        LRFF = sum(LRFF.index), 
        BR = count(BFF>0&LRFF>0) #number of times they occurred in the same plot
  ) 

## Identify sessions where both GHFF & LRFF were present:
GHFF_LRFF_plot <- treebat %>%
  filter(site.accession %in% occ.GHFF.output_plot$site.accession & site.accession %in% occ.LRFF.output_plot$site.accession) %>% #identify sessions where the species pair was present
  filter(rep %in% occ.GHFF.output_plot$rep | rep %in% occ.LRFF.output_plot$rep) %>% #identify occupied plots only
  filter(!is.na(GHFF.index)&!is.na(LRFF.index)) %>%
  ddply(c("site.code", "site.accession","subplot"), summarise,
        GHFF = sum(GHFF.index), 
        LRFF = sum(LRFF.index), 
        GR = count(GHFF>0&LRFF>0) #number of times they occurred in the same plot
  ) 

sum(BFF_GHFF_plot$BG) #number of plots where BFF and GHFF co-occurred
nrow(BFF_GHFF_plot)-sum(BFF_GHFF_plot$BG) #number of plots they were separate
#sum(BFF_GHFF_plot$BG)/nrow(BFF_GHFF_plot) #proportion of plots where BFF and GHFF co-occurred
BFF_GHFF <- binom.confint(x=sum(BFF_GHFF_plot$BG), n=nrow(BFF_GHFF_plot), methods = 'wilson') #plots where BFF and GHFF co-occurred
BFF_GHFF_nonocc <- binom.confint(x=(nrow(BFF_GHFF_plot)-sum(BFF_GHFF_plot$BG)), n=nrow(BFF_GHFF_plot), methods = 'wilson') #plots they were separate


sum(BFF_LRFF_plot$BR) #number of plots where BFF and LRFF co-occurred
sum(BFF_LRFF_plot$BR)/nrow(BFF_LRFF_plot) #proportion of plots where BFF and LRFF co-occurred
BFF_LRFF <- binom.confint(x=sum(BFF_LRFF_plot$BR), n=nrow(BFF_LRFF_plot), methods = 'wilson')
BFF_LRFF_nonocc <- binom.confint(x=(nrow(BFF_LRFF_plot)-sum(BFF_LRFF_plot$BR)), n=nrow(BFF_LRFF_plot), methods = 'wilson')

sum(GHFF_LRFF_plot$GR) #number of plots where GHFF and LRFF co-occurred
sum(GHFF_LRFF_plot$GR)/nrow(GHFF_LRFF_plot) #proportion of plots where GHFF and LRFF co-occurred
GHFF_LRFF <- binom.confint(x=sum(GHFF_LRFF_plot$GR), n=nrow(GHFF_LRFF_plot), methods = 'wilson')
GHFF_LRFF_nonocc <- binom.confint(x=(nrow(GHFF_LRFF_plot)-sum(GHFF_LRFF_plot$GR)), n=nrow(GHFF_LRFF_plot), methods = 'wilson')

### All species:
## From sessions where all species were present, how many plots did they co-occupy?
BFF_GHFF_LRFF_plot <- treebat %>%
  filter(site.accession %in% occ.BFF.output_plot$site.accession & site.accession %in% occ.GHFF.output_plot$site.accession & site.accession %in% occ.LRFF.output_plot$site.accession) %>% #identify sessions where the species pair was present
  filter(rep %in% occ.GHFF.output_plot$rep | rep %in% occ.LRFF.output_plot$rep | rep %in% occ.BFF.output_plot$rep) %>% #identify occupied plots only
  filter(!is.na(GHFF.index)&!is.na(LRFF.index)&!is.na(BFF.index)) %>%
  ddply(c("site.code", "site.accession","subplot"), summarise,
        BFF = sum(BFF.index),
        GHFF = sum(GHFF.index), 
        LRFF = sum(LRFF.index), 
        BGR = count(BFF>0&GHFF>0&LRFF>0)
  ) 

sum(BFF_GHFF_LRFF_plot$BGR) #number of plots where GHFF and LRFF and BFF co-occurred
sum(BFF_GHFF_LRFF_plot$BGR)/nrow(BFF_GHFF_LRFF_plot) #proportion of plots where GHFF and LRFF and BFF co-occurred
BFF_GHFF_LRFF <- binom.confint(x=sum(BFF_GHFF_LRFF_plot$BGR), n=nrow(BFF_GHFF_LRFF_plot), methods = 'wilson')
BFF_GHFF_LRFF_nonocc <- binom.confint(x=(nrow(BFF_GHFF_LRFF_plot)-sum(BFF_GHFF_LRFF_plot$BGR)), n=nrow(BFF_GHFF_LRFF_plot), methods = 'wilson')


#### Bar plot of plot co-occupation:
output <- rbind(BFF_GHFF, BFF_GHFF_nonocc, BFF_LRFF, BFF_LRFF_nonocc, GHFF_LRFF, GHFF_LRFF_nonocc, BFF_GHFF_LRFF, BFF_GHFF_LRFF_nonocc)
output$pair <- c("BFF_GHFF", "BFF_GHFF","BFF_LRFF", "BFF_LRFF","GHFF_LRFF","GHFF_LRFF", "X_BFF_GHFF_LRFF", "X_BFF_GHFF_LRFF")
output$occ <- c("Co-occur", "Separate","Co-occur", "Separate","Co-occur", "Separate", "Co-occur", "Separate")
output$label <- paste("Subplots = ", output$n)
output[c(2,4,6,8),9] <- ""

pair.lab <- c("BFF & GHFF", "BFF & LRFF", "GHFF & LRFF", "BFF & GHFF & LRFF")
names(pair.lab) <- c("BFF_GHFF","BFF_LRFF","GHFF_LRFF", "X_BFF_GHFF_LRFF")

plot_cooc <- ggplot(output, aes(x=occ)) + 
  geom_bar(aes(x=occ, y = mean, fill=pair), stat="identity") +
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  #scale_x_discrete(breaks=c("BFF_GHFF","BFF_LRFF","GHFF_LRFF"),
  #                 labels=c("BFF & GHFF", "BFF & LRFF", "GHFF & LRFF")) +
  geom_text(aes(label=label), y=1, color="black", size=6)+
  theme_bw() +
  scale_fill_manual(values=c("gray41", "gray63", "gray79", "slategray2"), labels = c("BFF & GHFF", "BFF & LRFF", "GHFF & LRFF", "BFF & GHFF & LRFF")) +
  background_grid(major="x", colour.major = "grey95")+
  labs(y="Proportion of subplots occupied", x="Species pair", colour="Species pair", fill="Species pair")+
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  theme(legend.position="none") +
  facet_grid(.~pair, labeller = labeller(pair = pair.lab)) +
  coord_cartesian(ylim=c(0,1))

jpeg("Output/Figures/Output/Figures/Figure 6_01.png", width=1000, height=500)
plot(plot_cooc)
dev.off()

jpeg("Output/Figures/Output/Figures/Figure 6.png", width=1000, height=1000)
ggarrange(plot_cooc, tree_cooc, nrow = 2, labels = c("A", "B"), font.label = list(size = 24))
dev.off()

####------------------------------------------------------------------------------------------##
##----------------------------------------- Figure 7 -----------------------------------------##
##--------------------------------- Height range of species ----------------------------------##
####------------------------------------------------------------------------------------------##

tree_meta <- read.csv("Data/Raw/tree metadata.csv") %>%
  rename(c("site.code"= "site.accession",
           "subplot"= "Plot",
           "Height" = "Tree_Height")) %>%
  mutate(subplot = as.factor(subplot))
### Calculate the average canopy height per site ###
canopy_plots <- treebat %>%
  dplyr::left_join(tree_meta[, c('site.code', 'subplot', 'tree.accession', 'Height')], by=c('site.code', 'subplot', 'tree.accession'))  %>% #there is an error here, because merged data is slighly longer (32,798) than the original data (32,785)
  filter(crown=="C") %>%
  filter(!is.na(Height)) %>%
  filter(tree.accession != "DLIS10461")  %>% #error in height here
  ddply(c("site.code", "subplot"), summarise,
        canopy_height = mean(Height)
  )
canopy_site <- canopy_plots %>%
  mutate(species = ifelse(site.code=="DLIS" & subplot == 4,"BFF", 
                          ifelse(site.code=="DLIS" & subplot == 5,"BFF",
                                 ifelse(site.code=="DLIS" & subplot == 6,"BFF",
                                        ifelse(site.code=="DLIS" & subplot == 7,"BFF",
                                               ifelse(site.code=="DLIS" & subplot == 8,"BFF",
                                                      ifelse(site.code=="DLIS" & subplot == 10,"BFF",
                                                             ifelse(site.code=="DLIS" & subplot == 9,"GHFF",
                                                                    ifelse(site.code=="DLIS" & subplot == 2,"GHFF",
                                                                           ifelse(site.code=="DLIS" & subplot == 1,"GHFF",
                                                                                  ifelse(site.code=="DLIS" & subplot == 3,"GHFF",
                                                                                         ifelse(site.code=="DCLU" & subplot == 6,"GHFF",
                                                                                                ifelse(site.code=="DCLU" & subplot == 7,"GHFF",
                                                                                                       ifelse(site.code=="DCLU" & subplot == 8,"GHFF",
                                                                                                              ifelse(site.code=="DCLU" & subplot == 1,"GHFF",
                                                                                                                     ifelse(site.code=="DCLU" & subplot == 2,"BFF",
                                                                                                                            ifelse(site.code=="DCLU" & subplot == 3,"BFF",
                                                                                                                                   ifelse(site.code=="DCLU" & subplot == 4,"BFF",
                                                                                                                                          ifelse(site.code=="DCLU" & subplot == 5,"BFF",
                                                                                                                                                 ifelse(site.code=="DCLU" & subplot == 9,"BFF",
                                                                                                                                                        ifelse(site.code=="DCLU" & subplot == 10,"BFF",
                                                                                                                                                               
                                                                                                                                                               "all"))))))))))))))))))))) %>%
  ddply(c("site.code", "species"), summarise,
        mean_canopy = mean(canopy_height))

data <- heights.subset %>%
  filter(species=="BFF"|species=="GHFF"|species=="LRFF") %>%
  create.Site()
x <- "session"
y <- "none"
ylab <- "Roosting height (m)"
xlab <- "Survey month"
title <- "Height range (minimum to maximum roosting height) of species"
flab <- "Species"
clab <- "Species"
colour <- "species"

## Facet by site, colour by species
## Create temp plot to extract smoothed values:
plot <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  stat_smooth(aes(y = min, colour=species), method = "loess", se = FALSE) +
  stat_smooth(aes(y = max, colour=species), method = "loess", se = FALSE)+
  scale_color_manual(values=c("lightsteelblue3", "seagreen4", "coral3"), labels = c("Black", "Grey-headed", "Little red"))+
  scale_fill_manual(values=c("lightsteelblue3", "seagreen4", "coral3"), labels = c("Black", "Grey-headed", "Little red"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019"))+
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered)) # .~ facet into collumns and # ~. facet into rows

## Extract data from loess fit:
df <- extract_loess_fac(plot)

## Re-plot with fill between loess fits
hrange_site_fac <- plot_blank(df, x, y, ylab, xlab, title, clab, flab) +
  geom_ribbon(aes(ymin = min, ymax = max, fill=species), alpha=0.6)+ 
  scale_color_manual(values=c("coral3", "lightsteelblue3", "gray49"), labels = c("Black", "Grey-headed", "Little red"))+
  scale_fill_manual(values=c("coral3", "lightsteelblue3", "gray49"), labels = c("Black", "Grey-headed", "Little red"))+
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered))  # .~ facet into collumns and # ~. facet into rows

## No facet, colour by species
## Create temp plot to extract smoothed values:
plot <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  stat_smooth(aes(y = min, colour=species), method = "loess", se = FALSE) +
  stat_smooth(aes(y = max, colour=species), method = "loess", se = FALSE)+
  scale_color_manual(values=c("coral3", "lightsteelblue3", "gray49"), labels = c("Black", "Grey-headed", "Little red"))+
  scale_fill_manual(values=c("coral3", "lightsteelblue3", "gray49"), labels = c("Black", "Grey-headed", "Little red"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019"))

## Extract data from loess fit:
df <- extract_loess_all(plot)

## Re-plot with fill between loess fits
hrange_site_all <- plot_blank(df, x, y, ylab, xlab, title, clab, flab) +
  geom_ribbon(aes(ymin = min, ymax = max, fill=species), alpha=0.6)+ 
  scale_color_manual(values=c("coral3", "lightsteelblue3", "gray49"), labels = c("Black", "Grey-headed", "Little red"))+
  scale_fill_manual(values=c("coral3", "lightsteelblue3", "gray49"), labels = c("Black", "Grey-headed", "Little red"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019"))

### Arrange and save plots:
jpeg("Output/Figures/Figure 7.png", width=1400, height=1200)
ggarrange(hrange_site_fac+coord_cartesian(ylim=c(5, 25)), hrange_site_all+coord_cartesian(ylim=c(5, 25))+ggtitle(""), nrow = 2, heights = c(3,1.5), labels = c("A", "B"), font.label = list(size = 24), common.legend = TRUE, legend = "bottom") #17_Comparison of roosting height between species_ZOOMED_1400 x 800
dev.off()

####------------------------------------------------------------------------------------------##
##----------------------------------------- Figure 8 -----------------------------------------##
##---------------------------- Sex ratio vs distance from center -----------------------------##
####------------------------------------------------------------------------------------------##

## Refresh data
treebat <-readRDS("Data/Processed/treebat.RDS")
centroids <- read.csv("Data/Raw/Centroids.csv") %>%
  mutate(roost.centroid.N = ifelse(site.accession=="DLIS010",NA,roost.centroid.N)) %>% ## Remove May 2019 Lismore area point (DLIS010)
  mutate(roost.centroid.E = ifelse(site.accession=="DLIS010",NA,roost.centroid.E)) 
centroids$session <- as.factor(centroids$session)

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

## Note - have plotted the proportion of each sex per, not just all M/all F/mixed because we didn't look at every bat per tree
x <- "scaled_eud.distance"
y <- "prop.M"
ylab <- "Proportion of male bats per tree"
xlab <- "Distance from roost center scaled by maximum distance value"
title <- "Proportion of male bats per tree (per tree, count of male/female is a representative subset of total tree)"
flab <- "Site"
clab <- "Site"
colour <- "Site"

### All species combined ###
##identify occupied subplots:
byplot <- c("site.code", "session", "site.accession", "subplot", "rep")
#bysite <- c("site.code", "session", "site.accession")
occ.all.output_plot <- occ.all(treebat, byplot)
#occ.BFF.output_plot <- occ.BFF(treebat, byplot)
#occ.GHFF.output_plot <- occ.GHFF(treebat, byplot)
#occ.LRFF.output_plot <- occ.LRFF(treebat, byplot)

data_sex <- treebat %>% 
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed
  #filter(rep %in% occ.all.output_plot$rep) %>% #don't want to filter unnocupied subplots in this case
  mutate(all.index.weight = BFF.index.weight+GHFF.index.weight+LRFF.index.weight)  %>% 
  mutate(prop.M.BFF = (BFF.M/(BFF.M+BFF.F))*100)  %>%  
  mutate(prop.M.GHFF = (GHFF.M/(GHFF.M+GHFF.F))*100)  %>%  
  mutate(prop.M.LRFF = (LRFF.M/(LRFF.M+LRFF.F))*100)  %>% 
  mutate(prop.M.all = ((BFF.M+GHFF.M+LRFF.M)/(BFF.M+GHFF.M+LRFF.M+BFF.F+GHFF.F+LRFF.F))*100)  %>% 
  create.Site() %>%
  ## Reformat so that there is one column to identify species:
  melt(id.vars = c("Site","site.code", "session", "site.accession", "scaled_eud.distance", "tree.accession", "subplot", "rep"), measure.vars = c("prop.M.BFF", "prop.M.GHFF", "prop.M.LRFF", "prop.M.all"),
       variable.name = c("species.value"), value.name="value") %>%
  dplyr::mutate(species = str_extract(species.value, "GHFF|BFF|LRFF|all")) %>%
  dplyr::mutate(valuecat = str_extract(species.value, "prop.M")) %>%
  dplyr::select(-c(species.value)) %>%
  pivot_wider(names_from = valuecat, values_from = value)

## No facet, no colour
data <- data_sex %>%
  filter(species == "all")

p2_centroid_site_all <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  #geom_point(aes(x=.data[[x]], y=.data[[y]], color=.data[[colour]]), size=1) +
  #scale_color_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore")) +
  geom_smooth(aes(x=.data[[x]], y=.data[[y]]), method = "loess", se = TRUE, linetype = "solid", size=0.5, colour="black") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) 

## Facet by site, dash by species 
dash <- "species"
dlab <- "Species"
flab <- "Species"

data_BFF <- data_sex %>%
  filter(species == "BFF")  
data_GHFF <- data_sex %>%
  filter(species == "GHFF")  
data_LRFF <- data_sex %>%
  filter(species == "LRFF")  
combined_data <- rbind(data_BFF, data_GHFF, data_LRFF)

p2_centroid_site_spp <- plot_dash_smooth(combined_data, x, y, dash, ylab, xlab, title, dlab, flab) +
  scale_linetype_manual(values=c("solid", "longdash", "dotted"), labels = c("Black", "Grey-headed", "Little red"))+
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered), scales="free_x") 

### Arrange and save plots:
jpeg("Output/Figures/Figure 8.png", width=1400, height=1200)
ggarrange(p2_centroid_site_spp+coord_cartesian(ylim=c(0, 100))+ggtitle(""), p2_centroid_site_all+coord_cartesian(ylim=c(0, 100))+ggtitle(""), nrow = 2, heights = c(3,1.5), labels = c("A", "B"), font.label = list(size = 24), common.legend = TRUE, legend = "bottom") #22-02_scaled distance from roost center vs proportion of males per tree_ZOOMED_1400 x 1200
dev.off()

####------------------------------------------------------------------------------------------##
##----------------------------------------- Figure 9 -----------------------------------------##
##----------------------------- Temporal patterns in occupation ------------------------------##
####------------------------------------------------------------------------------------------##

##################################################
### Proportion of male bats per tree over time ###
##################################################

### Set rectangles to show approximate timing of parturition for grey-headed flying-foxes and black flying-foxes (blue) and little red flying-foxes (red)
parturition_blue <- dplyr::tibble(month_start = 9, 
                                  month_end = 10)
parturition_blue <- dplyr::mutate(parturition_blue, 
                                  birth = dplyr::if_else(month_start == 9, "parturition", "not parturition"),
                                  start_session = 2,
                                  end_session = 3)

parturition_red <- dplyr::tibble(month_start = 4, 
                                 month_end = 5)
parturition_red <- dplyr::mutate(parturition_blue, 
                                 birth = dplyr::if_else(month_start == 4, "parturition", "not parturition"),
                                 start_session = 9,
                                 end_session = 10)

## Note - have plotted the proportion of each sex per, not just all M/all F/mixed because we didn't look at every bat per tree
x <- "session"
y <- "prop.M"
ylab <- "Proportion of male bats per tree"
xlab <- "Survey month"
title <- "Proportion of male bats per tree (per tree, count of male/female is a representative subset of total tree)"

### All species combined ###
flab <- "Site"
clab <- "Site"
colour <- "Site"
##identify occupied subplots:
byplot <- c("site.code", "session", "site.accession", "subplot", "rep")
#bysite <- c("site.code", "session", "site.accession")
occ.all.output_plot <- occ.all(treebat, byplot)
#occ.BFF.output_plot <- occ.BFF(treebat, byplot)
#occ.GHFF.output_plot <- occ.GHFF(treebat, byplot)
#occ.LRFF.output_plot <- occ.LRFF(treebat, byplot)

data <- treebat %>% 
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed
  #filter(rep %in% occ.all.output_plot$rep) %>% #don't want to filter unnocupied subplots in this case
  mutate(all.index.weight = BFF.index.weight+GHFF.index.weight+LRFF.index.weight)  %>% 
  mutate(prop.M.BFF = (BFF.M/(BFF.M+BFF.F))*100)  %>%  
  mutate(prop.M.GHFF = (GHFF.M/(GHFF.M+GHFF.F))*100)  %>%  
  mutate(prop.M.LRFF = (LRFF.M/(LRFF.M+LRFF.F))*100)  %>% 
  mutate(prop.M.all = ((BFF.M+GHFF.M+LRFF.M)/(BFF.M+GHFF.M+LRFF.M+BFF.F+GHFF.F+LRFF.F))*100)  %>% 
  create.Site() %>%
  ## Reformat so that there is one column to identify species:
  melt(id.vars = c("Site","site.code", "session", "site.accession", "scaled_eud.distance", "tree.accession", "subplot", "rep", "month", "year"), measure.vars = c("prop.M.BFF", "prop.M.GHFF", "prop.M.LRFF", "prop.M.all"),
       variable.name = c("species.value"), value.name="value") %>%
  dplyr::mutate(species = str_extract(species.value, "GHFF|BFF|LRFF|all")) %>%
  dplyr::mutate(valuecat = str_extract(species.value, "prop.M")) %>%
  dplyr::select(-c(species.value)) %>%
  pivot_wider(names_from = valuecat, values_from = value)

## No facet, no colour
data <- data_sex %>%
  filter(species == "all")

p7_temporal_all_all <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  geom_smooth(aes(y=.data[[y]]), method = "loess", se = TRUE, linetype = "solid", size=1, colour="black") +
  geom_rect(data = parturition_blue, aes(x = NULL, xmin = start_session, xmax = end_session, fill = birth),
            ymin = -0.5, ymax = Inf, alpha = 0.2) +
  geom_rect(data = parturition_red, aes(x = NULL, xmin = start_session, xmax = end_session, fill = birth),
            ymin = -0.5, ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("red","blue", "white"),
                    labels = c("LRFF parturition", "BFF & GHFF parturition", "not parturition")) +  
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019"))

### Split by species ###
## Facet by site, colour by species 
dash <- "species"
dlab <- "Species"
flab <- "Species"

data_BFF <- data_sex %>%
  filter(species == "BFF")  
data_GHFF <- data_sex %>%
  filter(species == "GHFF")  
data_LRFF <- data_sex %>%
  filter(species == "LRFF")  
combined_data <- rbind(data_BFF, data_GHFF, data_LRFF)

p7_temporal_all_spp <- plot_blank(combined_data, x, y, ylab, xlab, title, clab, flab) +
  geom_smooth(aes(as.numeric(as.character(.data[[x]])), y=prop.M, linetype=.data[[dash]]), method = "loess", se = TRUE, size=1, colour="black") +
  scale_linetype_manual(values=c("solid", "longdash", "dotted"), labels = c("Black", "Grey-headed", "Little red"))+
  #scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019")) +
  geom_rect(data = parturition_blue, aes(x = NULL, xmin = start_session, xmax = end_session, fill = birth),
            ymin = -0.5, ymax = Inf, alpha = 0.2) +
  geom_rect(data = parturition_red, aes(x = NULL, xmin = start_session, xmax = end_session, fill = birth),
            ymin = -0.5, ymax = Inf, alpha = 0.2) +
  scale_fill_manual(values = c("red","blue", "white"),
                    labels = c("LRFF parturition", "BFF & GHFF parturition", "not parturition")) +  
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered)) + # .~ facet into collumns and # ~. facet into rows
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019"))

### Arrange and save plots:
jpeg("Output/Figures/SI_proportion of males per tree over time.png", width=1550, height=1200)
ggarrange(p7_temporal_all_spp+coord_cartesian(ylim=c(0, 100))+ggtitle(""), p7_temporal_all_all+coord_cartesian(ylim=c(0, 100))+ggtitle("")+theme(axis.text.y = element_blank(), axis.title.y = element_blank()), ncol = 2, widths = c(3,0.6), labels = c("A", "B"), font.label = list(size = 24), common.legend = TRUE, legend = "bottom") #43-02_proportion of males per tree over time_1400 x 1200
dev.off()

jpeg("Output/Figures/SI_proportion of males per tree over time_02.png", width=1400, height=1200)
ggarrange(p7_temporal_all_spp+coord_cartesian(ylim=c(0, 100))+ggtitle(""), legend = "bottom") #43-02_proportion of males per tree over time_1400 x 1200
dev.off()


##################################################
###### Index total abundance over sessions #######
##################################################
data <- roostbat %>%
  filter(!is.na(site.code)) %>%
  mutate(missing.abundance = replace_na(pop.estimate.total, 0.1)) %>% #create additional column to plot missing data. 0.1 so that it can be plotted, and the real zero values are retained
  mutate(missing.abundance = replace(missing.abundance, missing.abundance>0.1, NA)) %>% #remove any rows with data. Left with only missing rows
  create.Site()
x <- "session"
y <- "index.abundance"
colour <- "Site"
ylab <- "Total roost abundance (index)"
xlab <- "Survey month"
flab <- "none"
title <- "Total roost abundance (index)"

## Site facet, colour by site
clab <- "Site"
colour <- "Site"
p1_temporal_fac_02 <- plot_col_smooth(data, x, y, colour, ylab, xlab, title, clab, flab) +
  #geom_point(aes(y=missing.abundance), shape=4, color="red", size=1) + #plot missing rows with red crosses
  #geom_line(aes(y=data[[y]])) +
  scale_color_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore"))+
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered)) +
  coord_cartesian(ylim=c(0, 8)) #fix scale

## No facet, point colour by site
clab <- "Site"
colour <- "Site"
p1_temporal_all <- plot_black_smooth(data, x, y, colour, ylab, xlab, title, clab, flab) +
  #geom_point(aes(y=missing.abundance), shape=4, color="red", size=1) + #plot missing rows with red crosses
  #geom_line(aes(y=data[[y]])) +
  scale_color_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019")) +
  coord_cartesian(ylim=c(0, 8)) #fix scale

### Arrange and save plots
jpeg("Output/Figures/SI_roost abundance over time.png", width=1400, height=1200)
ggarrange(p1_temporal_fac_02+ggtitle(""), p1_temporal_all+ggtitle(""), nrow = 2, heights = c(3,2), labels = c("A", "B"), font.label = list(size = 24), common.legend = TRUE, legend = "bottom") #01_Total roost abundance (index) over time_1400 x 800
dev.off()


##################################################
############ Total area over sessions ############
##################################################
data <- roostbat %>% filter(!is.na(site.code)) %>% create.Site()
x <- "session"
y <- "roost.area"
colour <- "Site"
ylab <- "Total roost area"
xlab <- "Survey month"
title <- "Total roost area"
clab <- "Site"
flab <- "none"

## Facet by site, colour by site
p2_temporal_fac <- plot_col_smooth(data, x, y, colour, ylab, xlab, title, clab, flab) +
  #geom_line(aes(y=data[[y]])) +
  scale_color_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore"))+
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019"))  +
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered)) 

## No facet, point colour by site
p2_temporal_all <- plot_black_smooth(data, x, y, colour, ylab, xlab, title, clab, flab) +
  #geom_line(aes(y=data[[y]])) +
  scale_color_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore"))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019")) 

### Arrange and save plots
jpeg("Output/Figures/SI_roost area over time.png", width=1400, height=1200)
ggarrange(p2_temporal_fac+coord_cartesian(ylim=c(0, 95000))+ggtitle(""), p2_temporal_all+coord_cartesian(ylim=c(0, 95000))+ggtitle(""), nrow = 2, heights = c(3,2), labels = c("A", "B"), font.label = list(size = 24), common.legend = TRUE, legend = "bottom") 
dev.off()


###################################################
########### Tree occupancy over sessions ##########
###################################################
## Note - points are each subplot per session
x <- "session"
y <- "none"
ylab <- "Proportion of trees occupied" #"Occupied trees (%)"
xlab <- "Survey month"
title <- "Proportion of trees occupied" #"Occupied trees (% per subplot, un-occupied subplots removed)"
flab <- "none"

### All species combined ###
data <- index.TREE.wide_plot %>% 
  filter(species == "all") %>% #extract combined species measure only
  filter(occ >0) %>% #remove unoccupied subplots
  create.Site()

## Facet by site, point colour by site:
clab <- "Site"
p5_temporal_all_fac <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) + #blank plot function because y axis is formula, not single collumn to call to
  #geom_point(aes(y=(occ/tree.count)*100, color=Site), size=1) +
  geom_smooth(aes(y=(occ/tree.count)*100, color=Site), method = "loess", se = TRUE, linetype = "solid", size=1) +
  scale_color_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore"))+
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019"))  +
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered)) 

## No facet, point colour by site
clab <- "Site"
p5_temporal_all_all <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  #geom_point(aes(y=(occ/tree.count)*100, color=Site), size=1) + 
  geom_smooth(aes(y=(occ/tree.count)*100), color="black", method = "loess", se = TRUE, linetype = "solid", size=1) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019")) +
  scale_color_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore"))


### Split by species ###
data <- index.TREE.wide_plot %>% 
  filter(species == "BFF"|species == "GHFF"|species == "LRFF") %>% 
  #filter(species == "BFF"|species == "GHFF") %>% #facets with LRFF won't plot, so alternative is to plot BFF and GHFF only
  filter(occ >0) %>% #remove unoccupied subplots
  create.Site()

## Facet by site, dash by species:
clab <- "Species"
p5_temporal_spp_fac <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  #geom_point(aes(y=(occ/tree.count)*100, color=species), size=1) +
  geom_smooth(aes(y=(occ/tree.count)*100, linetype=species), method = "loess", se = TRUE, colour = "black", size=1) +
  scale_linetype_manual(values=c("solid", "longdash", "dotted"), labels = c("Black", "Grey-headed", "Little red"))+
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019"))  +
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered)) 

### Arrange and save plots
jpeg("Output/Figures/SI_prop occupied trees over time.png", width=1400, height=1200)
ggarrange(p5_temporal_spp_fac+ggtitle(""), p5_temporal_all_all+ggtitle(""), nrow = 2, heights = c(3,2), labels = c("A", "B"), font.label = list(size = 24), common.legend = TRUE, legend = "bottom") 
dev.off()


##############################################################################
## Number of bats per OCCUPIED subplot, as given by the index approximation ##
##############################################################################
x <- "session"
y <- "N"
ylab <- "Within-plot abundance"
xlab <- "Survey month"
title <- "Within-plot abundance over time"
flab <- "none"
colour <- "Site"

### All species combined ###
data <- index.TREE.wide_plot %>% 
  filter(species == "all") %>% #extract combined species measure only
  filter(occ >0) %>% #remove unoccupied subplots
  create.Site()

## No facet
clab <- "Site"
p7_temporal_all <- plot_black_smooth(data, x, y, colour, ylab, xlab, title, clab, flab) +
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019"))  +
  scale_color_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore"))

## Facet by site, point colour by site:
clab <- "Site"
p7_temporal_all_fac <- plot_col_smooth(data, x, y, colour, ylab, xlab, title, clab, flab) +
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019"))  +
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered)) +
  scale_color_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore"))

## Facet by site, dash by species:
data <- index.TREE.wide_plot %>% 
  filter(species == "BFF"|species == "GHFF"|species == "LRFF") %>% 
  #filter(species == "BFF"|species == "GHFF") %>% #facets with LRFF won't plot, so alternative is to plot BFF and GHFF only
  filter(occ >0) %>% #remove unoccupied subplots
  create.Site()

clab <- "Species"
p7_temporal_spp_fac <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  geom_smooth(aes(y=.data[[y]], linetype=species), method = "loess", se = TRUE, colour = "black", size=1) +
  scale_linetype_manual(values=c("solid", "longdash", "dotted"), labels = c("Black", "Grey-headed", "Little red"))+
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019"))  +
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered)) 

### Arrange and save subplots
jpeg("Output/Figures/SI_plot abundance over time.png", width=1400, height=1200)
ggarrange(p7_temporal_spp_fac+ggtitle("")+coord_cartesian(ylim=c(0, 1400)), p7_temporal_all+ggtitle(""), nrow = 2, heights = c(3,2), labels = c("A", "B"), font.label = list(size = 24), common.legend = TRUE, legend = "bottom") 
dev.off()


###########################################################################
## Number of bats per OCCUPIED tree, as given by the index approximation ##
###########################################################################
##identify occupied subplots:
byplot <- c("site.code", "session", "site.accession", "subplot", "rep")
#bysite <- c("site.code", "session", "site.accession")
occ.all.output_plot <- occ.all(treebat, byplot)
#occ.BFF.output_plot <- occ.BFF(treebat, byplot)
#occ.GHFF.output_plot <- occ.GHFF(treebat, byplot)
#occ.LRFF.output_plot <- occ.LRFF(treebat, byplot)

x <- "session"
y <- "all.index.weight"
ylab <- "Within-tree abundance" #"Weighted index value (#bats)"
xlab <- "Survey month"
title <-"Within-tree abundance"  #"Number of bats per occupied tree (per tree, approximated with weighted index value)"
flab <- "none"

### All species combined ###
data <- treebat %>%
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed
  mutate(all.index.weight = BFF.index.weight+GHFF.index.weight+LRFF.index.weight) %>% #sum WEIGHTS, not index values, to give indication of overall bat load
  filter(rep %in% occ.all.output_plot$rep) %>% #choose occupied subplots only
  filter(all.index.weight > 0) %>% #choose occupied trees only
  create.Site()

## Facet by site, point colour by site:
colour <- "Site"
clab <- "Site"
p6_temporal_all_fac <- plot_col_smooth(data, x, y, colour, ylab, xlab, title, clab, flab) +
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019"))  +
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered)) +
  scale_color_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore"))

## No facet, point colour by site:
colour <- "Site"
clab <- "Site"
p6_temporal_all_all <- plot_black_smooth(data, x, y, colour, ylab, xlab, title, clab, flab) +
  scale_color_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore")) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019")) 

### Split by species ###

## Facet by site, colour by species 
##create separate dfs for each species:
data <- treebat %>%
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed
  mutate(all.index.weight = BFF.index.weight+GHFF.index.weight+LRFF.index.weight) %>% #sum WEIGHTS, not index values, to give indication of overall bat load
  filter(rep %in% occ.all.output_plot$rep) %>% #choose occupied subplots only
  create.Site()
treebat_BFF <- data %>%
  filter(BFF.index.weight>0)
treebat_GHFF <- data %>%
  filter(GHFF.index.weight>0)
treebat_LRFF <- data %>%
  filter(LRFF.index.weight>0)
clab <- "Species"
flab <- "Species"

p6_temporal_spp_fac <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  geom_point(data=treebat_BFF, aes(x=as.numeric(as.character(session)), y=BFF.index.weight), size=1, colour="lightsteelblue3", label="Black") + 
  geom_smooth(data=treebat_BFF, aes(x=as.numeric(as.character(session)), y=BFF.index.weight), method = "loess", se = TRUE, linetype = "solid", size=1, colour="lightsteelblue3", label="Black") +
  geom_point(data=treebat_GHFF, aes(x=as.numeric(as.character(session)), y=GHFF.index.weight), size=1, colour="seagreen4", label="Grey-headed") + 
  geom_smooth(data=treebat_GHFF, aes(x=as.numeric(as.character(session)), y=GHFF.index.weight), method = "loess", se = TRUE, linetype = "solid", size=1, colour="seagreen4", label="Grey-headed") +
  geom_point(data=treebat_LRFF, aes(x=as.numeric(as.character(session)), y=LRFF.index.weight), size=1, colour="coral3", label="Red") + 
  geom_smooth(data=treebat_LRFF, aes(x=as.numeric(as.character(session)),y=LRFF.index.weight), method = "loess", se = TRUE, linetype = "solid", size=1, colour="coral3", label="Red") +
  scale_color_manual(values=c("lightsteelblue3", "seagreen4", "coral3"), labels = c("Black", "Grey-headed", "Little red"))+
  facet_grid(.~Site, labeller = labeller(Site = site.labs.ordered)) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019"))

## Facet by species, colour by site
data <- treebat %>%
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed
  mutate(all.index.weight = BFF.index.weight+GHFF.index.weight+LRFF.index.weight) %>% #sum WEIGHTS, not index values, to give indication of overall bat load
  filter(rep %in% occ.all.output_plot$rep) %>% #choose occupied subplots only
  create.Site()
treebat_BFF <- data %>%
  filter(BFF.index.weight>0) %>%
  mutate(species="BFF") #NOTE that this isn't actually correct. Only correct if the BFF data is then extracted for plotting, like in the below. Have done as a quick fix for the facetting
treebat_GHFF <- data %>%
  filter(GHFF.index.weight>0) %>%
  mutate(species="GHFF") #NOTE that this isn't actually correct. Only correct if the GHFF data is then extracted for plotting
treebat_LRFF <- data %>%
  filter(LRFF.index.weight>0) %>%
  mutate(species="LRFF") #NOTE that this isn't actually correct. Only correct if the LRFF data is then extracted for plotting

clab <- "Site"
flab <- "Site"
p6_temporal_spp_all <- plot_blank(data, x, y, ylab, xlab, title, clab, flab) +
  geom_point(data=treebat_BFF, aes(x=as.numeric(as.character(session)), y=BFF.index.weight, colour=Site), size=1) + 
  geom_smooth(data=treebat_BFF, aes(x=as.numeric(as.character(session)), y=BFF.index.weight, colour=Site), method = "loess", se = TRUE, linetype = "solid", size=1) +
  geom_point(data=treebat_GHFF, aes(x=as.numeric(as.character(session)), y=GHFF.index.weight, colour=Site), size=1) + 
  geom_smooth(data=treebat_GHFF, aes(x=as.numeric(as.character(session)), y=GHFF.index.weight, colour=Site), method = "loess", se = TRUE, linetype = "solid", size=1) +
  geom_point(data=treebat_LRFF, aes(x=as.numeric(as.character(session)), y=LRFF.index.weight, colour=Site), size=1) + 
  geom_smooth(data=treebat_LRFF, aes(x=as.numeric(as.character(session)),y=LRFF.index.weight, colour=Site), method = "loess", se = TRUE, linetype = "solid", size=1) +
  scale_color_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba", "Canungra", "Clunes", "Lismore")) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019")) +
  facet_grid(.~species, labeller = labeller(species = spp.labs)) # .~ facet into collumns and # ~. facet into rows


#######################
## Arranged together ##
#######################
jpeg("Output/Figures/Figure 9.png", width=1400, height=1800)
ggarrange(p1_temporal_fac_02+ggtitle("")+theme(axis.text.x = element_blank(), axis.title.x = element_blank()), p2_temporal_fac+ggtitle("")+ theme(axis.text.x = element_blank(), axis.title.x = element_blank()), p5_temporal_all_fac+ggtitle("")+ theme(axis.text.x = element_blank(), axis.title.x = element_blank()),p7_temporal_all_fac+ggtitle(""), 
          nrow = 4, labels = c("A", "B", "C", "D"), font.label = list(size = 24), 
          heights = c(1,1,1,1.25),
          common.legend = TRUE, legend = "bottom", 
          align="v") 
dev.off()


