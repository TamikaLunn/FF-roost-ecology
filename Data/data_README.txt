OVERVIEW
This dataset presents data on three roosting structure of three species of flying-fox at eight sites in south-east Queensland and north-east New South Wales,Australia. 
Species included the grey-headed flying-fox (P. poliocephalus), black flying-fox (Pteropus alecto) and little red flying-fox (Pteropus scapulatus). 
All sites were previously documented as having a continuous population of grey-headed or black flying-foxes. 
Little red flying-foxes visited some roost sites intermittently, however no roost sites occurred within the distribution of spectacled flying-foxes.

RAW
We mapped the spatial arrangement of all overstory, canopy and midstory trees in a grid network of 10 stratified random subplots (20 x 20 meters each) per roost site. 
Subplots were stratified throughout perceived “core” (five subplots) and “peripheral” (five subplots) roosting areas, classed as areas observed to be frequently occupied (core) or infrequently (peripheral) by bats (Welbergen 2005). 
Core and peripheral areas were evaluated from regular observations made prior to roost tree mapping, though note that these categories were revised subsequently with the quantitative data. 
Trees were mapped and tagged using tree survey methods described in the “Ausplots Forest Monitoring Network, Large Tree Survey Protocol” (https://doi.org/10.1371/journal.pone.0137811).

To evaluate spatio-temporal patterns in flying-fox roosting, we revisited all tagged trees and scored the extent of species occupancy using the following tree abundance index: 0= zero bats; 1= 1-5 bats; 2=6-10 bats; 3=11-20 bats; 4=21-50 bats; 5=51-100 bats, 6=101-200 bats, 7= >200 bats. 
For a subset of trees (N=60 per site, consistent through time) absolute counts and minimum/maximum roosting heights of each species were taken. 
Overall roost perimeter (perimeter of area occupied) was mapped with GPS (accurate to 10 meters) immediately after the tree survey to estimate perimeter length and roost area. 
Total abundance at each roost was also estimated with a census count of bats where feasible (i.e., where total abundance was predicted to be <5,000 individuals), or by counting bats as they emerged in the evening from their roosts (“fly-out”), as per recommendations in Westcott et al. (2011). If these counts could not be conducted, population counts from local councils (conducted within ~a week of the bat surveys) were used, as the total abundance of roosts is generally stable over short timeframes (Nelson 1965b). 
Because roost estimates become more unreliable with increasing total abundance, and because our estimation methods were intrinsically linked with total abundance, we converted the total estimated abundance into an index estimate (where bin ranges increase with total abundance) for use in analyses, as per values used by the National Flying-Fox Monitoring Program (2017). Index categories were as follows: 1: 1-499 bats; 2: 500-2,499 bats; 3: 2,500-4,999 bats; 4: 5,000-9,999 bats; 5: 10,000-15,999 bats; 6: 16,000-49,999 bats; and 7: 50,000+ bats.
Roosting surveys were repeated once a month for 13 months (August 2018 - August 2019).
Methodological details are described in detail in the published paper 'Conventional wisdom on roosting behaviour of Australian flying foxes – a critical review, and evaluation using new data' (DOI https://doi.org/10.1002/ece3.8079). Raw data are available from this Dryad Dataset (https://doi.org/10.5061/dryad.g4f4qrfqv)

The files provided include:
tree survey data.csv - raw data on tree distributions measured within subplots:
     - tree.accession is unique tree ID. Name is structured as DAVO(site) 01(plot number) 001(tree number)
     - rand.tree indicates whether the tree was chosen for the subset measurements (Y) or not (N)
     - site.code is the roost (DAVO=Avondale, DBUR=Burleigh, DCAN=Canungra, DCLU=Clunes, DLIS=Lismore, DRED=Redcliffe, DSUN=Sunnybank, DTOW=Toowoomba)
     - subplot is the subplot number (1-10)
     - dead.alive is whether tree is dead (D) or alive (A). Dead trees are not assigned an "alive status" (below)
     - status is alive status of trees (A=alive, B=broken, C=leaning, D=fallen, E=fluted, F=hollow, G=rotten, H=multi-stemmed, I=few leaves, J=burnt, K=snapped, N=new, U=scar, V=dead top, X=buttress, Z=decline)
     - crown is crown class of tree (D=dominant, C=co-dominant, I=intermediate, CI=co-dominant but below canopy, S=suppressed, E=emergent, OG=open)
     - x is X distance within the SUBPLOT (0-20 meters only)
     - y is Y distance within the SUBPLOT (0-20 meters only)
     - Rx is X distance within the ROOST (0 ->20 meters, with (0,0) as SW corner)
     - Ry is X distance within the ROOST (0 ->20 meters, with (0,0) as SW corner)
     - E is UTM easting coordinate for tree. Calculated from Rx measure
     - N is UTM northing coordinate for tree. Calculated from Ry measure
spatial-bat-structure-data.csv - raw data on bat abundance per tree
     - tree.accession is unique tree ID. Name is structured as DAVO(site) 01(plot number) 001(tree number)
     - session is session number (1=August 2018, 2=September 2018, 3=October 2018, and so on)
     - site.accession is the unique roost*session ID number (DAVO=Avondale, DBUR=Burleigh, DCAN=Canungra, DCLU=Clunes, DLIS=Lismore, DRED=Redcliffe, DSUN=Sunnybank, DTOW=Toowoomba)
     - year/month/day of survey
     - BFF.index, GHFF.index, LRFF.index is extent of species occupancy scored with the tree abundance index (see above). BFF=black-flying-fox, GHFF=grey-headed flying-fox, LRFF=little ref flying-fox
     - BFF.count, GHFF.count, LRFF.count is exact count of species, measured for a subset of trees (see above)
     - BFF.min, BFF.max, GHFF.min, GHFF.max, LRFF.min, LRFF.max is minimum and maximum roosting hieghts of species, measured for a subset of trees (see above)
     - .M is number of male bats, .F is number of female bats, .Fw is number of female bats with young, .Fwo is number of female bats without young. Measured per species. Counts are for a portion of the tree only (hence totals do not add to count)
     - Additional columns are merged from the tree survey data (below). Tree survey data is repeated over sessions.
plot.data.csv -  data on trees per subplot
     - site.code is the roost (DAVO=Avondale, DBUR=Burleigh, DCAN=Canungra, DCLU=Clunes, DLIS=Lismore, DRED=Redcliffe, DSUN=Sunnybank, DTOW=Toowoomba)
     - subplot is the subplot number (1-10)
     - Ntrees is total number of trees measured in the subplot
     - E, D, C, CI, I, OG are the number of trees tallied per crown class in the subplot (D=dominant, C=co-dominant, I=intermediate, CI=co-dominant but below canopy, S=suppressed, E=emergent, OG=open)
     - core.peripheral is whether the plot was considered a core or peripheral area at the beginning of the tree survey, evaluated from regular observations made prior to roost tree mapping (though note that these categories were revised subsequently with the quantitative data)
     - mean.dist is the mean distance between trees in the subplot
     - plotID is a unique plot ID. Name is structured as DAVO(site) 001(plot number)
pixel-density-data.csv & pixel-density-data-nonzero.csv - data on fixed-bandwidth weighted kernel estimates calculated in spatstat from the tree distribition and bat occupancy data (see below)
     - v is the intensity/density value
     - Collumns give mean, min, max, median, lower interquartile range, upper interquartile range, and sd for values of v over the subplot, per species (BFF, GHFF, LRFF) or per species combined (all)
     - Summaries of v were calculated with either all pixel values in the subplot (pixel-density-data.csv) or with zero values excluded (pixel-density-data-nonzero.csv). Zero values represent empty space within the roost where there were no trees (and no bats)
Centroids.csv - the easting (roost.centroid.E) and northing (roost.centroid.N) of the spatial center of the roost, based on which trees were occupied at the time of the survey
ALL-roost-use-data.csv - roost-level information on bat occupancy
     - site.code is the roost (DAVO=Avondale, DBUR=Burleigh, DCAN=Canungra, DCLU=Clunes, DLIS=Lismore, DRED=Redcliffe, DSUN=Sunnybank, DTOW=Toowoomba)
     - session is session number (1=August 2018, 2=September 2018, 3=October 2018, and so on)
     - site.accession is the unique roost*session ID number
     - date/year/month/day of survey
     - pop.estimate.total is the total abundance at each roost was also estimated with a census count of bats where feasible (i.e., where total abundance was predicted to be <5,000 individuals), or by counting bats as they emerged in the evening from their roosts (“fly-out”), as per recommendations in Westcott et al. (2011). If these counts could not be conducted, population counts from local councils (conducted within ~a week of the bat surveys) were used, as the total abundance of roosts is generally stable over short timeframes (Nelson 1965b). 
     - index.abundance is the index estimate: Because roost estimates become more unreliable with increasing total abundance, and because our estimation methods were intrinsically linked with total abundance, we converted the total estimated abundance into an index estimate (where bin ranges increase with total abundance) for use in analyses, as per values used by the National Flying-Fox Monitoring Program (2017). Index categories were as follows: 1: 1-499 bats; 2: 500-2,499 bats; 3: 2,500-4,999 bats; 4: 5,000-9,999 bats; 5: 10,000-15,999 bats; 6: 16,000-49,999 bats; and 7: 50,000+ bats.
     - roost.area, roost.perimeter, pop.estimate.total, and index.abundance are measured for all species combined. 
     - pop.estimate.BFF, pop.estimate.GHFF, pop.estimate.LRFF are species specific population estimates (not always feasible to obtain).
     - pop.method is the method used to generate the estimate (see above)
/Voronoi Polygons/tree.tessellation.csv - the crown area estimated per tree, computed as the area of Dirichlet-Voronoi tessellations from tree distribution maps of canopy trees per subplot, calculated with the spatstat package in R (Baddeley 2010). Visuals of the tree distribution (circles) and tesselations (polygons) per plot are also given.
     - value_dirichlet is the calculated crown area
     - tree.accession is unique tree ID. Name is structured as DAVO(site) 01(plot number) 001(tree number)
     - site.code is the roost (DAVO=Avondale, DBUR=Burleigh, DCAN=Canungra, DCLU=Clunes, DLIS=Lismore, DRED=Redcliffe, DSUN=Sunnybank, DTOW=Toowoomba)
     - subplot is the subplot number (1-10)
     - Crown is crown class of tree (D=dominant, C=co-dominant, I=intermediate, CI=co-dominant but below canopy, S=suppressed, E=emergent, OG=open)
     - x/y is the X and Y distance within the SUBPLOT (0-20 meters only)


PROCESSED DATA
Information collected during the bat roosting surveys were used to calculate measures of bat density and abundance at three scales: roost-level, subplot-level and tree-level. 
For a visual summary of metrics see Figure 2 in the published paper ('Counterintuitive scaling between population abundance and local density: implications for modelling transmission of infectious diseases in bat populations'). 
Note that where index abundance scores were used in calculations, the middle value of the index range was taken. 

Roost-level density was calculated by dividing the total roost index abundance score by the total roost area (Figure 2A). 
Measures of subplot-level density were estimated with two methods: either as the tally of tree-level index abundance scores per subplot divided by subplot area (“subplot-level density”,  Figure 2B), or as the average of fixed-bandwidth weighted kernel estimates, estimated using the spatstat package in R (Diggle 1985) (“subplot-level kernel density”, Figure 2C). Kernel estimates are spatially explicit and give the density of a spatial pattern, estimated per pixel over a smoothed area (Baddeley 2010). Kernels were estimated from the spatial location of trees weighted by tree-level index abundance scores, with Gaussian kernel smoothing and a smoothing bandwidth of 0.6. Bandwidth was selected by comparing projected kernel density values to expected density values based on tree abundance and canopy area. Kernel averages were then calculated per subplot. To prevent dilution of the density estimates with unoccupied space, we included only occupied pixels in the subplot average (pixel size = 0.156 x 0.156 meters). This latter approach has the advantage of explicitly incorporating the spatial distribution of bats into the density estimate, and therefore gives a better representation of aggregations in occupied space. Note that neither roost nor subplot-based density measures consider the vertical distribution of bats.
Measures of tree-level density were estimated in either two-dimension (2-D; for comparison with other two-dimensional estimates) or three-dimension (3-D). Tree-level 2-D density was estimated from tree-level index abundance scores and canopy area (Figure 2D). Tree-level 3-D density was estimated for the tree subset, as the absolute count of bats divided by the volume of tree space occupied (i.e. per cubic metre rather than square metre, Figure 2E). Volume of tree space was calculated from the height range occupied (maximum height minus minimum height) and the approximate crown area of trees. To control for edge effects, and to prevent overestimation of crown area for overstory trees and trees outside of the canopy, we imposed a maximum crown area of 199 m2 (radius ~8 m). This value was selected based on mean values reported across species of eucalypts in New South Wales (Verma et al. 2014), eucalypts being broadly representative of trees in these roost sites (Brooks 2020). In total, 218 of the 2,522 tagged trees (8%) were imposed with the maximum crown area value. Crown area of midstory trees was assigned as the first quartile of canopy tree crown area (5.8 m2), to reflect observations that trees beneath the canopy were typically smaller than trees within the canopy. Mean calculated crown area was 30.4 m2 (crown radius ~ 3.1 m). To investigate whether the choice of maximum crown area impacted results, we also repeated analyses for additional values of maximum crown area (140 m2, 170 m2 and 230 m2) chosen to cover the range in smallest to largest mean values reported for individual eucalypt species in Verma et al. (2014).

The data provided here are the processed density and abundance values (response/dependent variables), and roost features (predictor/independent variables) calculated from raw data, used in the manuscript 'Counterintuitive scaling between population abundance and local density: implications for modelling transmission of infectious diseases in bat populations'
treebat.csv - is the re-formatted spatial-bat-structure-data.csv
index.TREE.wide_plot.csv - is the re-formatted plot.data.csv combined with the re-formatted pixel-density-data.csv and re-formatted pixel-density-data-nonzero.csv (with subplot-level kernel density estimates)
roostbat.csv - is the re-formatted ALL-roost-use-data.csv
kernel.TREE.wide_plot - is the re-formatted pixel-density-data.csv, with subplot-level kernel density estimates
kernelnz.TREE.wide_plot - is the re-formatted pixel-density-data-nonzero.csv, with subplot-level kernel density estimates
heights.subset - is the re-formatted spatial-bat-structure-data.csv, for the tree subset only, with tree-level 2D and 3D density estimates



