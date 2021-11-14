# FF-roost-ecology
Data and code repository for Lunn et al (2021) Conventional wisdom on roosting behavior of Australian flying-foxes — A critical review, and evaluation using new data <https://doi.org/10.1002/ece3.8079>

To test commonly held understandings about the roosting ecology of Australian flying-foxes (and inform management practices of flying-fox roosts), we contribute a contemporary, fine-scale dataset on within-roost roosting structure, encompassing 13-monthly repeat measures from 2,522 spatially referenced roost trees across eight sites in southeastern Queensland and northeastern New South Wales. 

The main statistical comparisons tested with our empirical data were as follows: (1) whether density of occupation is greater for subplots in “core” areas of the roost compared with subplots in irregularly occupied “peripheral” areas (defined by occupation greater than or less than 80% of surveys respectively (Appendix S1); (2) whether bat occupation decreases with distance from the roost center (per species); (3) whether bat species segregate in vertical space; and (4) whether dominant individuals occupy the center of roosts, and subdominant individuals the outer area (per species). We also provide qualitative comparisons of (5) seasonal patterns of abundance and occupancy per species; and (6) whether bat species segregate in horizontal space. We utilized generalized additive mixed models (GAMMs) for all statistical comparisons.

The code is organised into an RProject containing six files:
- 00_functions_VSubmission.R - behind-the-scenes custom functions to help with data processing and visualisation
- 01_create-kernel-density-and-Voroni-polygons_VSubmission.R - spatial code to 1) create kernel density estimates from tree spatial point patterns weighted by bat index values in trees, and 2) calculate and visualise crown area of trees in subplots
- 02_format-raw-data_VSubmission.R - code to format the raw data into structures suitible for analyes, and to save processed data for later reference
- 03_fit-models_VSubmission.R - code to run general addative mixed models (GAMMs) for main statistical comparisons
- 04_generate-figures_VSubmission.R - code to generate figures presented in the published manuscript
- 05_bat-male-distribution-visuals_VSubmission.RMD - code to generate Appendix S4 of paper - Male composition per tree over space and time
- 06_bat-structure-cluster-visuals_VSubmission.RMD - code to generate Appendix S3 of paper - Species density over space and time
The rendered RMD files are available at:
- Appendix S3 ("Species density over space and time"): https://drive.google.com/file/d/1VNM5c91lPlf-aUKOcGhKU8TvtxQHIHiO/view?usp=sharing
- Appendix S4 ("Male composition per tree over space and time"): https://drive.google.com/file/d/1kfXRLUKsrlfNU_0I5SW2nutdAGqzt7KT/view?usp=sharing

** It is important to load the code as the RProject so that root folders in setwd are correct **

Data are organised into an RProject containing two folders
- /Raw : raw data input for code '01_create-kernel-density-and-Voroni-polygons_VSubmission.R' and '02_format-raw-data_VSubmission.R'
- /Processed : processed data input for code '03_fit-models_VSubmission.R' and '04_generate-figures_VSubmission.R'
More information on each file and data within files is given in the data_README file

Output is organised into an RProject containing:
- /Figures : figures presented in the main text or supplementary information in the published paper
- /Model output : GAMM model output for main statistical comparisons
More details on model structure are given in the published manuscript

