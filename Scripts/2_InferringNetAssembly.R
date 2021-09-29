###################################
## Code to replicate the analyses and figures shown in: Trait-based inference of ecological network assembly
## Ecological Monographs
# Authors: 
# Emma-Liina Marjakangas1,2,*, Gabriel Muñoz3,*, Shaun Turney3,*, Jörg Albrecht4, Eike Lena Neuschulz4, Matthias Schleuning4, Jean-Philippe Lessard3
# 
# 1 Centre for Biodiversity Dynamics, Department of Biology, Norwegian University of Science and Technology, N-7491 Trondheim, Norway
# 2 Finnish Museum of Natural History, University of Helsinki, P.O. Box 17, Helsinki FI-00014, Finland
# 3 Department of Biology, Faculty of Arts and Sciences, Concordia University, Montreal, Canada
# 4 Senckenberg Biodiversity and Climate Research Centre (SBiK-F), Senckenberganlage 25, 60325 Frankfurt am Main, Germany
# * These authors contributed equally to this study
###################################

#####
# Contact: Gabriel Muñoz - gabriel.munoz@concordia.ca
#####

## Load custom functions (important!)
source("Scripts/0_Functions.R") # path might change 



################
# Section 5: INFERRING NETWORK ASSEMBLY 
################




######
# a) Process-based species pool delineation
#####

# Slope and elevation are derived from a digital elevation model obtained with radar satellite remote sensing (NASA SRTM shuttle) (https://cgiarcsi.community/data/srtm-90m-digital-elevation-database-v4-1/)
# Temperature and precipitation correspond to the variables BIO01 and BIO12 of WORLDCLIM (https://developers.google.com/earth-engine/datasets/catalog/WORLDCLIM_V1_BIO#bands)

# Load libraries and corresponding data
library(gdistance)
library(ade4)
library(FD)
library(picante)
library(Rmisc)
library(spatstat)
library(RStoolbox)
library(MASS)

# load raw coordinates of collection data 
locCoord <- readxl::read_xlsx("data_package/data/CoordinatesEcuador.xlsx") ## path can change
# get centroids for each of collected sites
centroids <- data.frame("lon" = aggregate(locCoord$Lon, by = list(locCoord$Plot), mean)$x,
                        "lat" = aggregate(locCoord$Lat, by = list(locCoord$Plot), mean)$x,
                        "site" = unique(locCoord$Plot)[order(unique(locCoord$Plot))])

# Load environmental resistance layers (exported from earthengine, https://code.earthengine.google.com/14fe67bd4749c1f90a253f6cfaa05ff9)

# path can change 
SlopeLoja <- raster::raster("data_package/data/LojaSlope.tif")
ClimLoja <- raster::raster("data_package/data/LojaClim.tif")
PrecLoja <- raster::raster("data_package/data/LojaPrec.tif")
ElevLoja <- raster::raster("data_package/data/LojaElevation.tif")


# Let's apply the custom made function to calculate cost-path distances in
# environmental space and standardize the distance matrices based on their range maxima 

eDistSlope <- CustomCostDist(SlopeLoja, centroids)
eDistTemp <- CustomCostDist(ClimLoja, centroids)
eDistPrec <- CustomCostDist(PrecLoja, centroids)
eDistElev <- CustomCostDist(ElevLoja, centroids)

# normalize variables

slopeDist <- decostand(eDistSlope, "range")
ClimDist <- decostand(eDistTemp, "range")
eDistPrec <- decostand(eDistPrec, "range")
eDistElev <- decostand(eDistElev, "range")

## Geographical distance matrix 
LCPMat <- decostand((slopeDist+ClimDist+eDistPrec+eDistElev), "range")

## Generating environmental distance matrix 

# extract env data from site (added buffer)
EnvDat <- data.frame("Slope" = extract(SlopeLoja,round(centroids[,1:2], 3)),
                     "Clim" = extract(ClimLoja,round(centroids[,1:2], 3)),
                     "Prec" = extract(PrecLoja,round(centroids[,1:2], 3)),
                     "Elev" = extract(ElevLoja,round(centroids[,1:2], 3)))

# PCA to reduce dimensionality and extract scores of 1 axis to calculate euc.distances

SiteScor <- scores(rda(EnvDat))$sites[,1]
names(SiteScor) <- centroids$site
EnvMat <- dist(SiteScor)
# normalize
EnvMat <- decostand(EnvMat, "range")

# Get Probability matrix to sample species from sites by assigning equal weights to dispersal and establishment 

ProbPool <- ((1-EnvMat)*0.5 + (1-LCPMat)*0.5)

png("ElevLoja.png", 
    500, 500, pointsize = 20)
par(mar = c(2,2,2,2))
plot(ElevLoja, box = F, frame  = F, yaxt = "n",
     col = RColorBrewer::brewer.pal(4, "Greys")
)
#axis(2)
points(round(centroids[,1:2], 3))

dev.off()

png("SlopeLoja.png", 
    500, 500, pointsize = 20)
par(mar = c(2,2,2,2))
plot(SlopeLoja, box = F, frame  = F, yaxt = "n",
     col = RColorBrewer::brewer.pal(4, "Greys")
)
points(round(centroids[,1:2], 3))

dev.off()
png("ClimLoja.png", 
    500, 500, pointsize = 20)
par(mar = c(2,2,2,2))
plot(ClimLoja/10, box = F, frame  = F,yaxt = "n",
     col = RColorBrewer::brewer.pal(4, "Greys")
)
points(round(centroids[,1:2], 3))
dev.off()
png("PrecLoja.png", 
    500, 500, pointsize = 20)
par(mar = c(2,2,2,2))
plot(PrecLoja, box = F, frame  = F,yaxt = "n",
     col = RColorBrewer::brewer.pal(4, "Greys")
)
points(round(centroids[,1:2], 3))
dev.off()

image_append(c(
image_trim(image_read("ElevLoja.png")),
image_trim(image_read("SlopeLoja.png")),
image_trim(image_read("ClimLoja.png")),
image_trim(image_read("PrecLoja.png"))),
F)




png("EnvMat.png", 
    500, 500, pointsize = 20)
par(oma = c(3,3,3,3))
heatmap(1-EnvMat)
dev.off()

png("LCPMat.png", 
    500, 500, pointsize = 20)
par(oma = c(3,3,3,3))
heatmap(1-LCPMat)
dev.off()

png("ProbPool.png", 
    500, 500, pointsize = 20)
par(oma = c(3,3,3,3))
heatmap(ProbPool)
dev.off()

######
# a) Applying the null model approach to the empirical data 
#####


## RQL scaling 



###############################################################################
# load packages and data
################################################################################
library(bipartite)
library(ade4)
library(FD)
library(picante)
load("data_package/data/dataSet1.RData")
################################################################################
# check data
################################################################################
str(dataSet1)

data.frame(dataSet1$N[match(dataSet1$pT$plantCode, dataSet1$N$plantCode[-c(1:4)]),][c("site", "elevation")],
           dataSet1$N[match(dataSet1$pT$plantCode, dataSet1$N$plantCode[-c(1:4)]),][c("site", "elevation")])


dataSet1$pT[c("plantCode", "fruitLength")]
# pT - plant traits
# aT - animal traits
# N - plant-animal interactions
################################################################################
# get some first impression: plot interaction network on each elevation
# and count number of links and species
################################################################################
N <- xtabs(frequency ~ plantCode + animalCode + elevation, data = dataSet1$N)
par(mfrow = c(3, 1))
apply(N, 3, plotweb)
apply(N, 3, FUN = function(x) sum(x > 0)) # number of unique links on each elevation
apply(N, 3, FUN = function(x) sum(rowSums(x) > 0)) # number of plant species on each elevation
apply(N, 3, FUN = function(x) sum(colSums(x) > 0)) # number of animal species on each elevation
################################################################################
# setup data matrices for RLQ analysis based on the pooled binary metaweb
# (as in Albrecht et al. 2018, Nat. Comm. 9:3177)
################################################################################
pTRLQ <- sqrt(dataSet1$pT[, 2:5]) # sqrt-transform traits to improve normality
aTRLQ <- sqrt(dataSet1$aT[, 2:5]) # sqrt-transform traits to improve normality
nRLQ <- xtabs(frequency ~ plantCode + animalCode, data = dataSet1$N)
nRLQ <- nRLQ[rownames(pTRLQ), rownames(aTRLQ)]
nRLQ <- decostand(nRLQ, "pa")
################################################################################
# perform correspondence and principal components analyses
################################################################################
coa <- dudi.coa(as.data.frame(nRLQ), scannf = FALSE, nf = 4)
duP <- dudi.pca(as.data.frame(pTRLQ), scannf = FALSE, center = TRUE, scale = TRUE, nf = 4, row.w = coa$lw)
duA <- dudi.pca(as.data.frame(aTRLQ), scannf = FALSE, center = TRUE, scale = TRUE, nf = 4, row.w = coa$cw)
################################################################################
# perform RLQ analysis
################################################################################
RLQ <- rlq(duP, coa, duA, scannf = FALSE, nf = 4)
RLQ

par(oma = c(4,4,4,4))
heatmap(as.matrix(RLQ$tab))
################################################################################
# perform permutation tests
################################################################################
# Global test for association between trait spaces of plants and animals
################################################################################
trRLQ <- fourthcorner2(as.data.frame(pTRLQ),
                       as.data.frame(nRLQ),
                       as.data.frame(aTRLQ),
                       p.adjust.method.G = "fdr",
                       modeltype = 6, nrepet = 9999)$trRLQ
trRLQ
################################################################################
# Test for axis-specific correlations between plant and animal trait spaces
# NOTE: only tests of 1st vs. 1st and 2nd vs. 2nd axis are meaningful
################################################################################
fcRLQ <- fourthcorner.rlq(RLQ, type = "axes",
                          p.adjust.method.G = "fdr",
                          p.adjust.method.D = "fdr",
                          p.adjust.D = "levels",
                          modeltype = 6, nrepet = 9999)
fcRLQ
################################################################################
# Test for association of animal traits with first and second axis of animal trait space
################################################################################
qAxes <-  fourthcorner.rlq(RLQ, type = "Q.axes",
                           p.adjust.method.G = "fdr",
                           p.adjust.method.D = "fdr",
                           p.adjust.D = "levels",
                           modeltype = 6, nrepet = 9999)
qAxes
################################################################################
# Test for association of plant traits with first and second axis of plant trait space
################################################################################
rAxes <- fourthcorner.rlq(RLQ, type = "R.axes",
                          p.adjust.method.G = "fdr",
                          p.adjust.method.D = "fdr",
                          p.adjust.D = "levels",
                          modeltype = 6, nrepet = 9999)
rAxes
################################################################################
################################################################################
# Figure 1. Plot of RLQ-fourth-corner analysis
################################################################################
png("TraitMatching.png", width = 700, height = 500, pointsize = 14, res = 92)
m <- matrix(rep(c(1, 2), each = 4), 1, 8, byrow = TRUE)
m <- apply(m, 2, rep, each = 4)
m[16] <- 3
layout(m)
par(cex = 0.7, mar = c(0.5, 0.5, 0.5, 0.5) + 0.1, oma = c(2, 2, 2, 0))
p <- RLQ$mR
a <- RLQ$mQ
N <- nRLQ
p <- cbind(p, (rowSums(N > 0)/max(rowSums(N > 0)) + 0.25) * 2)
a <- cbind(a, (colSums(N > 0)/max(colSums(N > 0)) + 0.25) * 2)
p <- p[order(p[, 3], decreasing = TRUE), ]
a <- a[order(a[, 3], decreasing = TRUE), ]
Ntab <- as.data.frame.table(N > 0)
Ntab$pX <- p[as.character(Ntab$plantCode), 1]
Ntab$pY <- p[as.character(Ntab$plantCode), 2]
Ntab$aX <- a[as.character(Ntab$animalCode), 1]
Ntab$aY <- a[as.character(Ntab$animalCode), 2]
Ntab <- subset(Ntab, Freq > 0)
plot.new()
plot.window(xlim = c(-1, 1) * 1.1, ylim = c(-1, 1) * 1.1, asp = 1)
abline(h = 0, v = 0, col = "gray")
arrows(0, 0, RLQ$c1[, 1], RLQ$c1[, 2], col = "firebrick1",
       code = 2, angle = 30, length = 0.05, lwd = 1, lty = c(1, 1, 5, 2))
arrows(0, 0, RLQ$l1[, 1], RLQ$l1[, 2], col = "deepskyblue1",
       code = 2, angle = 30, length = 0.05, lwd = 1, lty = c(1, 1, 5, 2))
text(RLQ$c1 * 1.1, gsub("_", " ", rownames(RLQ$c1), fixed = TRUE),
     col = "firebrick1", font = 1, cex = 1, xpd = TRUE)
text(RLQ$l1 * 1.1, gsub("_", " ", rownames(RLQ$l1), fixed = TRUE),
     col = "deepskyblue1", font = 1, cex = 1, xpd = TRUE)
mtext(side = 3, line = 0, adj = -0.1, cex = 1, letters[1], font = 2)
legend("topright", bty = "n", lty = c(1, 5, 2), title = "Trait type",
       text.col = c("black", "black", "black"),
       legend = c("Matching", "Energy", "Foraging"))
plot.new()
plot.window(xlim = range(c(p[, 1], a[, 1])), ylim = range(c(p[, 2], a[, 2])), asp = 1)
with(Ntab, arrows(pX, pY, aX, aY, code = 0, lwd = 0.5, col = gray(0.75)))
points(p[, 1:2], pch = 21, cex = p[, 5], lwd = 0.5,
       col = "black", bg = "deepskyblue1")
points(a[, 1:2], pch = 21, cex = a[, 5], lwd = 0.5,
       col = "black", bg = "firebrick1")
mtext(side = 3, line = 0, adj = 0, cex = 1, letters[2], font = 2)
legend("topleft", bty = "n",
       text.col = c("deepskyblue1", "firebrick1", "gray"),
       legend = c("Plants", "Animals", "Interactions"))
par(mar = c(0.5, 0.5, 0.5, 0.5) + 0.1)
barplot(RLQ$eig/sum(RLQ$eig), ylim = c(0, 1), axes = FALSE)
mtext(side = 3, line = 0, adj = 0.5, cex = 0.5, "Eigenvalues")
axis(2, at = seq(0, 1, 1), las = 1, mgp = c(1.5, 0.75, 0))
axis(1, labels = c(1:4),at = seq(1, 4, 1), las = 1,tick = F, mgp = c(1.5, 0.75, 0))
dev.off()


################################################################################
################################################################################
# the scores from the RLQ analysis can be used to calculate functional diversity
# indices on each elevation and conduct null model analyses
################################################################################
# extract normed species scores of plants and animals for further analysis
# these scores represent the position of plants and animals in the
# respective trait spaces
################################################################################
# Animals - extract normed species scores in trait space & create distance matrix
################################################################################
RLQ$mQ
aDis <- dist(RLQ$mQ)
################################################################################
# Plants - extract normed species scores in trait space & create distance matrix
################################################################################
RLQ$mR
pDis <- dist(RLQ$mR)
################################################################################
## How interactions look into a trait space of RQL 
################################################################################
intRQL <- dplyr::bind_rows(RLQ$mQ,RLQ$mR)
rownames(intRQL) <- c(rownames(RLQ$mQ), rownames(RLQ$mR))
intRQL <- data.frame(dataSet1$N, "RLQan" = RLQ$mQ[match(dataSet1$N$animalCode,rownames(RLQ$mQ)),][,1])
intRQL <- data.frame(intRQL, "RLQpla" = RLQ$mQ[match(dataSet1$N$plantCode,rownames(RLQ$mR)),][,1])
################################################################################
pool <- data.frame(
  dataSet1$N,
  dataSet1$pT[match(dataSet1$N$plantCode, dataSet1$pT$plantCode),][,c(2:5)],
  dataSet1$aT[match(dataSet1$N$animalCode, dataSet1$aT$animalCode),][,c(2:5)],
  "RQLa" = intRQL$RLQan[match(dataSet1$N$animalCode, intRQL$animalCode)],
  "RQLp" = intRQL$RLQpla[match(dataSet1$N$plantCode, intRQL$plantCode)])
pool$intID <- paste0(pool$plantCode,pool$animalCode)


################################################################################
### c) Sensitivity analyses to different metrics to calculate SES 
##############################################################################


# Plants
# randomized pool -sd
myEF_r_sd_p <- test4EF(pool = pool,
                       traitName = "RQLp", 
                       nrep = 100, 
                       side = "P",
                       sdOrRange = "sd",
                       RoP = "Random",
                       ProbPool = ProbPool )
myEF2_r_sd_a <- test4EF(pool = pool, 
                        traitName = "RQLa",
                        nrep = 100, 
                        side = "A",
                        sdOrRange = "sd",
                        RoP = "Random", 
                        ProbPool = ProbPool )
# randomized pool -range
myEF_r_ran_p <- test4EF(pool = pool, 
                        traitName = "RQLp", 
                        nrep = 100, 
                        side = "P",
                        sdOrRange = "range",
                        RoP = "Random", 
                        ProbPool = ProbPool )



myEF2_r_ran_a <- test4EF(pool = pool, 
                         traitName = "RQLa",
                         nrep = 100,
                         side = "A",
                         sdOrRange = "range",
                         RoP = "Random", 
                         ProbPool = ProbPool )
# probabilistic pool -sd
myEF_pro_sd_p <- test4EF(pool = pool, 
                         traitName = "RQLp",
                         nrep = 100, 
                         side = "P", 
                         sdOrRange = "sd",
                         RoP = "Prob",
                         ProbPool = ProbPool )
myEF_pro_sd_p1 <- test4EF(pool = pool, 
                         traitName = "RQLp",
                         nrep = 100, 
                         side = "P", 
                         sdOrRange = "sd",
                         RoP = "Prob",
                         ProbPool = ProbPool )
myEF2_pro_sd_a <- test4EF(pool = pool, 
                          traitName = "RQLa", 
                          nrep = 100, 
                          side = "A", 
                          sdOrRange = "sd",
                          RoP = "Prob", 
                          ProbPool = ProbPool )
myEF2_pro_sd_a1 <- test4EF(pool = pool, 
                          traitName = "RQLa", 
                          nrep = 100, 
                          side = "A", 
                          sdOrRange = "sd",
                          RoP = "Prob", 
                          ProbPool = ProbPool )

# probabilistic pool -range
myEF_pro_ran_p <- test4EF(pool = pool, 
                          traitName = "RQLp", 
                          nrep = 100, 
                          side = "P", 
                          sdOrRange = "range",
                          RoP = "Prob", 
                          ProbPool = ProbPool )
myEF2_pro_ran_a <- test4EF(pool = pool, 
                           traitName = "RQLa",
                           nrep = 100, side = "A", 
                           sdOrRange = "range",
                           RoP = "Prob", 
                           ProbPool = ProbPool )
#####


EFSizeDF <- data.frame(myEF_r_sd_p,
                       myEF2_r_sd_a,
                       myEF_r_ran_p,
                       myEF2_r_ran_a,
                       myEF_pro_sd_p,
                       myEF2_pro_sd_a,
                       myEF_pro_ran_p,
                       myEF2_pro_ran_a)


########################
## Step 2. The influence of niche-based assembly processes on network structure
########################
## Start parallelization 
#####
library(snowfall)

sfInit(parallel=TRUE, cpus=10, type="SOCK")  # start cluster 
sfSource("Scripts/0_Functions.R") # Source functions to the cluster
# Export the general pool to sample species
sfExport("pool")
# Export the probabilistic matrix
sfExport("ProbPool")
## Run null networks in parallel
#####
####################
#CnRoIn
####################
CnRoIn <- sfLapply(1:length(unique(pool$site)), function(x) 
  replicate(10,NullNetMod(pool,
                          ProbPool, 
                          repPool = 1000,
                          NPlant =  F, 
                          NAnimal =  T, 
                          nulltype = "DD",
                          unique(pool$site)[x], 100)))


names(CnRoIn) <- unique(pool$site)
####################
# CoRoIn
####################
CoRoIn <- sfLapply(1:length(unique(pool$site)), function(x) 
  replicate(10,NullNetMod(pool,
                          ProbPool, 
                          repPool = 1000,
                          NPlant =  F, 
                          NAnimal =  F, 
                          nulltype = "DD",
                          unique(pool$site)[x], 100)))


names(CoRoIn) <- unique(pool$site)


####################
####################
#CnRnIn
####################
CnRnIn <- sfLapply(1:length(unique(pool$site)), function(x) 
  replicate(10,NullNetMod(pool,
                          ProbPool, 
                          repPool = 1000,
                          NPlant =  T, 
                          NAnimal =  T, 
                          nulltype = "DD",
                          unique(pool$site)[x], 100)))
names(CnRnIn) <- unique(pool$site)


####################

####################
#CoRnIn
####################
CoRnIn <- sfLapply(1:length(unique(pool$site)), function(x) 
  replicate(10,NullNetMod(pool,
                          ProbPool, 
                          repPool = 1000,
                          NPlant =  T, 
                          NAnimal =  F, 
                          nulltype = "DD",
                          unique(pool$site)[x], 100)))





# CoRnIn <- CoRnIn[c(1,3,5,6,8,10)] 
names(CoRnIn) <- unique(pool$site)


############
# Get the mean of replicates
######

CnRnIn_x <- sapply(1:6, function(x) apply(CnRnIn[[x]], 1, mean))
colnames(CnRnIn_x)<- names(CnRnIn)

CoRoIn_x <- sapply(1:6, function(x) apply(CoRoIn[[x]], 1, mean))
colnames(CoRoIn_x)<- names(CoRoIn)

CnRoIn_x <- sapply(1:6, function(x) apply(CnRoIn[[x]], 1, mean))
colnames(CnRoIn_x)<- names(CnRoIn)

CoRnIn_x <- sapply(1:6, function(x) apply(CoRnIn[[x]], 1, mean))
colnames(CoRnIn_x)<- names(CoRnIn)

## Calculating confidence intervals 

sapply(1:2, function(x) Rmisc::CI(CnRnIn_x[x,], 0.95 ))

############
# Get the sd of replicates
######

CnRnIn_sd <- sapply(1:6, function(x) apply(CnRnIn[[x]], 1, sd))
colnames(CnRnIn_x)<- names(CnRnIn)

CoRoIn_sd <- sapply(1:6, function(x) apply(CoRoIn[[x]], 1, sd))
colnames(CoRoIn_x)<- names(CoRoIn)

CnRoIn_sd <- sapply(1:6, function(x) apply(CnRoIn[[x]], 1, sd))
colnames(CnRoIn_x)<- names(CnRoIn)

CoRnIn_sd <- sapply(1:6, function(x) apply(CoRnIn[[x]], 1, sd))
colnames(CoRnIn_x)<- names(CoRnIn)




###########
# Process strength assymetry 
###########


PSdir <- data.frame("PSdir_random_range" = log(abs(myEF2_r_ran_a)/abs(myEF_r_ran_p)),
                    "PSdir_random_sd" = log(abs(myEF2_r_sd_a)/abs(myEF_r_sd_p)),
                    "PSdir_process_range" = log(abs(myEF2_pro_ran_a)/abs(myEF_pro_ran_p)),
                    "PSdir_process_sd" = log(abs(myEF2_pro_sd_a)/abs(myEF_pro_sd_p)))


stPSS <- PSdir$PSdir_process_sd

stPSS <-stPSS/max(abs(stPSS))

#####
# Calculate SESnet 
#####

# calculate observed modularity and nestedness 
netSites <- xtabs(frequency~plantCode + animalCode + site,pool)
nullDist <- apply(netSites, 3, NullModSen,100, "DD")
# make binary to calculate observed modularity
netSites <- apply(netSites, 3, makeBinar)
obsModular <- sapply(1:6, function(x)
  bipartite::computeModules(netSites[[x]])@likelihood)
obsNested <- sapply(1:6, function(x) vegan::nestednodf(netSites[[x]])$statistic[["NODF"]])
Zscores <- sapply(1:6, function(x) (obsModular[x] -nullDist[[x]][[1]]["NulMod.mean"]) / nullDist[[x]][[1]]["NulMod.sd"])
Zscores_nes <- sapply(1:6, function(x) (obsNested[x] -nullDist[[x]][[1]]["NulNest.mean"]) / nullDist[[x]][[1]]["NulNest.sd"])
plot(Zscores,Zscores_nes)
nullDist[[1]][[1]]["NulNest.mean"]
names(Zscores) <- unique(pool$site)
names(Zscores_nes) <- unique(pool$site)

# Calculate effect sizes from null networks 

# modularity 
CnRnIn_SES_mod <- (Zscores-CnRnIn_x[1,]) / CnRnIn_sd[1,]
CoRoIn_SES_mod <- (Zscores-CoRoIn_x[1,]) / CoRoIn_sd[1,]
CnRoIn_SES_mod <- (Zscores-CnRoIn_x[1,]) / CnRoIn_sd[1,]
CoRnIn_SES_mod <- (Zscores-CoRnIn_x[1,]) / CoRnIn_sd[1,]

# nestedness
CnRnIn_SES_nes <- (Zscores_nes-CnRnIn_x[2,]) / CnRnIn_sd[2,]
CoRoIn_SES_nes <- (Zscores_nes-CoRoIn_x[2,]) / CoRoIn_sd[2,]
CnRoIn_SES_nes <- (Zscores_nes-CnRoIn_x[2,]) / CnRoIn_sd[2,]
CoRnIn_SES_nes <- (Zscores_nes-CoRnIn_x[2,]) / CoRnIn_sd[2,]



## Partialying out separate effects of consumer, resources, and interactions 

# modularity 

INTef_SES_mod <- CnRnIn_SES_mod + CoRnIn_SES_mod -CnRoIn_SES_mod -  CoRoIn_SES_mod
INT_SES_mod <- CoRoIn_SES_mod
R_SES_mod <- CoRnIn_SES_mod- CoRoIn_SES_mod
C_SES_mod <- CnRoIn_SES_mod- CoRoIn_SES_mod

# nestedness 
INTef_SES_nes <- CnRnIn_SES_nes + CoRnIn_SES_nes -CnRoIn_SES_nes -  CoRoIn_SES_nes
INT_SES_nes <- CoRoIn_SES_nes
R_SES_nes <- CoRnIn_SES_nes- CoRoIn_SES_nes
C_SES_nes <- CnRoIn_SES_nes- CoRoIn_SES_nes

## Create a new dataframe with ef sizes

# modularity 
SESnet_mod <- data.frame(
  INTef_SES_mod, 
  INT_SES_mod,
  R_SES_mod,
  C_SES_mod)

#SESnet_mod <- SESnet_mod/max(abs(SESnet_mod))
colnames(SESnet_mod) <- c("InterEf", "Interactions", "Resources", "Consumers")


# nestednes 
SESnet_nes <- data.frame(
  INTef_SES_nes, 
  INT_SES_nes,
  R_SES_nes,
  C_SES_nes)

#SESnet_nes <- SESnet_nes/max(abs(SESnet_nes))
colnames(SESnet_nes) <- c("InterEf", "Interactions", "Resources", "Consumers")



# find partial sum of effects (no accounting for the interaction between effects)

# modularity 
sumAbsEf <- abs(SESnet_mod$Resources)+abs(SESnet_mod$Consumers)+abs(SESnet_mod$Interactions)
# modularity 
sumAbsEf_nes <- abs(SESnet_nes$Resources)+abs(SESnet_nes$Consumers)+abs(SESnet_nes$Interactions)

# Normalize data by the partial sum to find relative magnitudes

# modularity
SESnet_mod_norm <- SESnet_mod[,c(2,3,4)]/(sumAbsEf)

SESnet_mod_norm$elev <- dataSet1$N$elevation[match(rownames(SESnet_mod_norm),
                           dataSet1$N$site)]

# nestedness 

SESnet_nes_norm <- SESnet_nes[,c(2,3,4)]/(sumAbsEf_nes)

SESnet_nes_norm$elev <- dataSet1$N$elevation[match(rownames(SESnet_nes),
                                                   dataSet1$N$site)]



apply(SESnet_mod, 2, Rmisc::CI, 0.95)
apply(SESnet_nes, 2, Rmisc::CI, 0.95)


SESnet_mod <- SESnet_mod[match(rownames(PSdir),rownames(SESnet_mod)),]



SESnet_mod$elev <- dataSet1$N$elevation[match(rownames(SESnet_mod),
                                              dataSet1$N$site)]
SESnet_nes$elev <- dataSet1$N$elevation[match(rownames(SESnet_nes),
                                              dataSet1$N$site)]






siteElev <- c("Bombuscaro", "Copalinga", "ECSF", "Finca", "Bellavista", "Cajanuma")


par(mfrow = c(3,2), mar = c(2,0,2,0), oma = c(0,0,0,0))
for(i in rev(match( siteElev,dimnames(N2)$site))){
  graph2 <- N2[,,i]
  graph2 <- makeBinar(graph2)
  graph2.1 <- igraph::graph_from_incidence_matrix(graph2)

  dim(N2)[i]
  
  graph <- unique(intRQL[intRQL$site == dimnames(N2)$site[i],])
  a <- graph$RLQan
  b <- graph$RLQpla
  
  a_c <- aggregate(graph$RLQan, list(graph$animalCode), mean)$x
  b_c <- aggregate(graph$RLQpla, list(graph$plantCode), mean)$x
  
  graph_co <- unique(intRQL)
  a_co <- graph$RLQan
  b_co <- graph$RLQpla
  
  vertex_col <- f(c(a_co, b_co), 6, "BrBG", T)
  vertex_col <- vertex_col[match(c(a_c,b_c),c(a_co, b_co))]
  
  
  
  
  f1 <- kde2d(a, b, n = 100)
  f1$z <- ifelse(f1$z > 0.7, 0.7,f1$z)
  
  cor <-  data.frame(match(round(a),round(f1$x)),
                     match(round(b),round(f1$y)))
  
  intPro <- sapply(1:length(cor[,1]), function(x) f1$z[as.numeric(cor[x,][1]),
                                             as.numeric(cor[x,][2])])
  graph2.1 <-delete.vertices(simplify(graph2.1), degree(graph2.1)==0)

 

ord <- order(c(rep("#256B92",dim(graph2)[1] ),
                    rep("#B22222", dim(graph2)[2])),
             c(a_c,b_c))
  
  l <- layout.circle(graph2.1, order = ord)
  plot(graph2.1,
       vertex.size = 10,
       edge.curved=T,
       xlim = c(-3,3),
       edge.color = scales::alpha("grey", 0.9),
       main = dimnames(N2)$site[i],
       vertex.frame.color = scales::alpha(c(rep("#256B92",dim(graph2)[1] ),
                                            rep("#B22222", dim(graph2)[2])),0.5), 
       vertex.color =scales::alpha(c(rep("#256B92",dim(graph2)[1] ),rep("#B22222", dim(graph2)[2])),0.7),
       #vertex.color = f(c(a_c,b_c), 5,"PuOr", T),
       vertex.shape = c(rep("circle",dim(graph2)[1] ),rep("circle", dim(graph2)[2])),
       #edge.color = ifelse(intPro>0.5, "red", "grey50"),
       vertex.label = "",
       edge.width = 0.4,
       layout = l)
  par(new = T)
  igraph::plot.igraph(graph2.1,
       vertex.size = 10,
       edge.curved=T,
       xlim = c(-2.4,2.4),
       
       main = dimnames(N2)$site[i],
       vertex.frame.color = vertex_col, 
       #vertex.color =scales::alpha(c(rep("#256B92",dim(graph2)[1] ),rep("#B22222", dim(graph2)[2])),0.7),
       vertex.color = vertex_col,
       edge.width = 0,
       vertex.shape = c(rep("circle",dim(graph2)[1] ),rep("circle", dim(graph2)[2])),
       edge.color = "grey50",
       vertex.label = "",
       layout = l)
  
  

}






par(las = 1)
plot(c(-5:1)~rep(0,length(c(-5:1))),
     pch = 15,
     xlab= "",
     main = "Trait value",
     ylab = "",
     cex =10,
     xlim = c(-0.01,0.01),
     frame = F,
     axes = F,
     col =  f(c(-5:1), 6,"BrBG", F))
axis(side = 2, line = -16, gap.axis = 2)
  

par(las = 1)
plot(c(seq(0,1,0.1))~rep(0,
                        length(c(seq(0,1,0.1)))),
     pch = 15,
     xlab= "",
     main = "Trait value",
     ylab = "",
     cex =10,
     xlim = c(-0.01,0.01),
     frame = F,
     axes = F,
     ylim = c(0,1),
     col =  f(c(seq(0,1.2,0.1)), 7,"YlOrRd", T))
axis(side = 2, line = -16, gap.axis = 2)

  
  ########
  ## Figure 4:Effect of process strength asymmetry on plant-frugivore network structure along an elevational gradient
  #####


#### Figure 4.a) Strength of niche-based community assembly and assymetry processes for plants and frugivores across elevation
#############

png("FilteringGradient.png", width = 1000, height = 1000, pointsize = 20, res = 110)
par(mfrow = c(1,1), oma = c(2,2,2,2), mar = c(5,5,2,2))

plot(-myEF_pro_sd_p,elev_vec1, frame= F,
     xlim = c(-3,3),
     pch = 21,
     ylim  = c(0,4000),
     cex.lab = 1.5,
     cex.axis = 1,
     xlab = "Strength of community assembly",
     ylab = "Elevation",
     bg = scales::alpha("#B22222", 1),
     cex = 1.5)
segments(0,elev_vec1,-myEF_pro_sd_p,elev_vec1, col = "#B22222", lwd = 2)
abline(v=0)
abline(h = c(500,1500, 2500,3500), lty = 2)
points(-myEF2_pro_sd_a,c(elev_vec2),
       ylim = c(-3,8), 
       pch = 21,
       xlim  = c(0,4000),
       bg = scales::alpha("#256B92", 1),
       cex = 1.5)
segments(0,elev_vec2,-myEF2_pro_sd_a,elev_vec2, col = "#256B92" ,lwd = 2)

mtext("A",3, outer = T, adj = 0 , cex = 2)
legend(1.7,1600, "*", bty = "n", cex = 2)
mtext("Trait overdispersion",3, outer = F, adj = 0 , cex = 1)
mtext("Trait filtering",3, outer = F, adj = 1 , cex = 1)


dev.off()

png("PSA.png", width = 1000, height = 1000, pointsize = 20, res = 110)
par(mfrow = c(1,1), oma = c(2,2,2,2), mar = c(5,5,2,2))

plot(elev_vec2[order]~(-stPSS)[order],
     pch = 21,
     frame = F,
     xlab = "Process strength asymmetry",
     ylab = "Elevation",
     cex = 1.5,
     cex.lab = 1.5,
     xlim = c(-1,1),
     bg = ifelse(stPSS[order]>0, "#256B92", "#B22222"), 
     ylim = c(0,4000))

segments(0,elev_vec2[order],(-stPSS)[order],
         elev_vec2[order], 
         col = ifelse(stPSS[order]>0, "#256B92", "#B22222"), 
         lwd = 2)

abline(h = c(500,1500,2500,3500), v = 0,lty = 2)
mtext("B",3, outer = T, adj = 0 , cex = 2)
mtext("Consumer driven",3, outer = F, adj = 1 , cex = 1)
mtext("Resource driven",3, outer = F, adj = 0 , cex = 1)


dev.off()
png("StructureChange.png", width = 1000, height = 1000, pointsize = 20, res = 110)

par(mfrow = c(1,1), oma = c(2,2,2,2), mar = c(5,5,2,2))

elev_vec <- dataSet1$N$elevation[match(names(myEF_pro_sd_p1),dataSet1$N$site)]
elev_vec1 <- c(elev_vec-c(-150, -150, +250, +350, +250, -150))
elev_vec2 <- c(3250,1250, 2850, 750, 1850, 2250)


order <- match(names(Zscores),rownames(SESnet_mod))
plot(Zscores_nes,elev_vec2[order], 
     xlim = c(-50,50), 
     frame= F,
     pch = 1,
     ylim  = c(0,4000),
     cex.lab = 1.5,
     cex.axis = 1,
     xlab = "Network structure (Z-scores)",
     ylab = "Elevation",
     col = scales::alpha("#9C413D", 1),
     cex = 1.5)
segments(0,elev_vec2[order],Zscores_nes,elev_vec2[order], col = "#9C413D", lwd = )

points(Zscores,
       elev_vec1[match(names(Zscores),rownames(SESnet_mod))],
       pch = 16,       
       cex = 1.5,
       col = scales::alpha("#9C413D", 1))
abline(h = c(500,1500,2500,3500), v = 0,lty = 2)
segments(0,elev_vec1[order],Zscores,elev_vec1[order], col = "#9C413D", lwd = 2)

mtext("Random",3, outer = F, adj = 0.5 , cex = 1)

legend("bottomright", 
       pch=c(16,1),
       col = "#9C413D",
       bty = "n",
       legend = c("Modularity","Nestedness"))
mtext("C",3, outer = T, adj = 0 , cex = 2)


dev.off()








  elev_vec <- dataSet1$N$elevation[match(names(myEF_pro_sd_p1),dataSet1$N$site)]
  elev_vec1 <- c(elev_vec-c(-150, -150, +250, +350, +250, -150))
  elev_vec2 <- c(3250,1250, 2850, 750, 1850, 2250)
  
  
  
  
  
  
  png("AssemblyMode.png", width = 1000, height = 1000, pointsize = 20, res = 110)
  par(mfrow = c(1,1), oma = c(2,2,2,2), mar = c(5,5,2,2))
  
  assMode <- log(abs(SESnet_mod$Consumers)/abs(SESnet_mod$Resources))
  assMode2 <- log(abs(SESnet_nes$Consumers)/abs(SESnet_nes$Resources))
  
  plot(elev_vec2[order]~assMode[order],
       pch = 16,
       frame = F,
       xlab = "Network assembly mode",
       ylab = "Elevation",
       cex = 1.5,
       cex.lab = 1.5,
       ylim = c(0,4000),
       col = "#9C413D",
       xlim = c(-7,7))
  segments(0,elev_vec2[order],assMode[order],
           elev_vec2[order], 
           col = "#9C413D",
           lwd = 2)
  
  points(elev_vec1[order]~assMode2[order],
         col = "#9C413D",
         cex = 1.5,
         pch = 1)
  
  segments(0,elev_vec1[order],assMode2[order],
           elev_vec1[order], 
           col = "#9C413D", 
           lwd = 1)
  abline(h = c(500,1500,2500,3500), v = 0,lty = 2)
  mtext("D",3, outer = T, adj = 0 , cex = 2)
  mtext("Bottom-up",3, outer = F, adj = 0 , cex = 1)
  mtext("Top-down",3, outer = F, adj = 1 , cex = 1)
  
  legend("bottomright", 
         pch=c(16,1),
         col = "#9C413D",
         bty = "n",
         legend = c("Modularity","Nestedness"))
  dev.off()
 
  








  elev_vec <- dataSet1$N$elevation[match(names(myEF_pro_sd_p1),dataSet1$N$site)]
  elev_vec1 <- c(elev_vec-c(-150, -150, +250, +350, +250, -150))
  elev_vec2 <- c(3250,1250, 2850, 750, 1850, 2250)
  
  
  
  
  
  
  png("AssemblyMode.png", width = 1000, height = 1000, pointsize = 20, res = 110)
  par(mfrow = c(1,1), oma = c(2,2,2,2), mar = c(5,5,2,2))
  
  assMode <- log(abs(SESnet_mod$Consumers)/abs(SESnet_mod$Resources))
  assMode2 <- log(abs(SESnet_nes$Consumers)/abs(SESnet_nes$Resources))
  
  plot(elev_vec2[order]~assMode[order],
       pch = 16,
       frame = F,
       xlab = "Network assembly mode",
       ylab = "Elevation",
       cex = 1.5,
       cex.lab = 1.5,
       ylim = c(0,4000),
       col = "#9C413D",
       xlim = c(-7,7))
  segments(0,elev_vec2[order],assMode[order],
           elev_vec2[order], 
           col = "#9C413D",
           lwd = 2)
  
  points(elev_vec1[order]~assMode2[order],
         col = "#9C413D",
         cex = 1.5,
         pch = 1)
  
  segments(0,elev_vec1[order],assMode2[order],
           elev_vec1[order], 
           col = "#9C413D", 
           lwd = 1)
  abline(h = c(500,1500,2500,3500), v = 0,lty = 2)
  mtext("D",3, outer = T, adj = 0 , cex = 2)
  mtext("Bottom-up",3, outer = F, adj = 0 , cex = 1)
  mtext("Top-down",3, outer = F, adj = 1 , cex = 1)
  
  legend("bottomright", 
         pch=c(16,1),
         col = "#9C413D",
         bty = "n",
         legend = c("Modularity","Nestedness"))
  dev.off()
 
  
  
  
  # write the final image 
  
  magick::image_write(
    magick::image_append(c(
      magick::image_append(c(
        magick::image_append(c(
          magick::image_read("FilteringGradient.png"),
          magick::image_read("PSA.png"))))),
      magick::image_append(c(
        magick::image_append(c(
          magick::image_read("StructureChange.png"),
          magick::image_read("AssemblyMode.png")))))),
      stack = T),
    "Figure4Final.png")
  
  
  
  
  
  magick::image_write(
        magick::image_append(c(
          magick::image_read("FilteringGradient.png"),
          magick::image_read("StructureChange.png"),
          magick::image_read("AsymmetryEffects.png"))),
    "Figure4Final2.png")
  
  ####
  # Appendix images
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # 
  # png("AsymmetryEffects.png", width = 1000, height = 1000, pointsize = 20, res = 110)
  # par(mfrow = c(1,1), oma = c(2,2,2,2), mar = c(5,5,2,2))
  # plot(log(abs(SESnet_mod$Consumers)/abs(SESnet_mod$Resources))~
  #        stPSS,
  #      cex = 1.5,
  #      col = "#9C413D",
  #      pch = (SESnet_mod$elev/1000)+14,
  #      frame = F,
  #      cex.axis = 1, 
  #      xlim = c(-1,1),
  #      ylim = c(-4,8),
  #      cex.lab = 1.5,
  #      xlab = "Process strength asymmetry",
  #      ylab = "Network assembly mode")
  # 
  # 
  # 
  # points(log(abs(SESnet_mod$Consumers)/abs(SESnet_mod$Resources))~
  #          stPSS, cex= 1.5,
  #        pch = (SESnet_mod$elev/1000)-1)
  # 
  # 
  # points(log(abs(SESnet_nes$Consumers)/abs(SESnet_nes$Resources))~
  #          stPSS[match(rownames(SESnet_nes),rownames(SESnet_mod))], 
  #        pch = (SESnet_nes$elev/1000)-1,
  #        cex = 1.5, col = "#9C413D")
  # abline(h=0,v=0,lty = 2)
  # legend(0.5,8, 
  #        pch = c(15,0),
  #        c("Modularity","Nestedness"),
  #        cex = 0.7,
  #        col = c("#9C413D", "#9C413D"),
  #        bty = "n")
  # legend(0.40,8, 
  #        c("",""),
  #        cex = 0.7,
  #        pch = c(15,0)+1,
  #        col = c("#9C413D", "#9C413D"),
  #        bty = "n")
  # legend(0.30,8, 
  #        c("",""),
  #        cex = 0.7,
  #        pch = c(15,0)+2,
  #        col = c("#9C413D", "#9C413D"),
  #        bty = "n")
  # 
  # 
  # legend("bottomright",
  #        title = "Elevation",
  #        c("3000","2000","1000"),
  #        pch = c(17,16,15),
  #        bty = "n")
  # legend("bottomleft", 
  #        "Bottom-up",
  #        bty = "n" )
  # 
  # legend("topleft", 
  #        "Top-down",
  #        bty = "n" )
  # 
  # mtext("C",3, outer = T, adj = 0 , cex = 2)
  # 
  # 
  # dev.off()
  # 
  
  
  
  
  # Visualize sites in environmental space
  png("Maps.png", width = 500, height = 500, pointsize = 17)
  par(mfrow = c(2,2), mar = c(2,3,2,4), oma = c(2,2,2,2))
  plot(SlopeLoja, main = "Slope")
  points(centroids$lon,centroids$lat, pch = "+", cex = 2)
  plot(ClimLoja, main = "Temperature")
  points(centroids$lon,centroids$lat, pch = "+", cex = 2)
  plot(PrecLoja, main = "Precipitation")
  points(centroids$lon,centroids$lat, pch = "+", cex = 2)
  plot(ElevLoja, main = "Elevation")
  points(centroids$lon,centroids$lat, pch = "+", cex = 2)
  mtext("A",3, outer = T, adj = 0 , cex = 2)
  
  
  dev.off()
  
  # How many components are relevant for each variable? 
  
  png("PCAs.png", width = 500, height = 500, pointsize = 17)
  par(mfrow = c(2,2), mar = c(2,3,2,4), oma = c(2,2,2,2))
  plot(pcaSlope, main = "Slope")
  plot(pcaTemp, main = "Temperature")
  plot(pcaPrec, main = "Precipitation")
  plot(pcaElev, main = "Elevation")
  mtext("B",3, outer = T, adj = 0 , cex = 2)
  
  
  dev.off()
  
  # Let's visualize how our probabilisitic species pools are delineated among sites. 
  png("RDA.png", width = 500, height = 500, pointsize = 15)
  par(mfrow = c(1,1), mar = c(2,3,2,4), oma = c(2,2,2,2))
  biplot(EnvPCA, 
         xlim = c(-2,2),
         ylim = c(-1,1),
         frame  = F) # first axis = 59% of variation 
  mtext("C",3, outer = T, adj = 0 , cex = 2)
  
  
  dev.off()
  
  # Visualize the delination of species pool selection 
  png("Pool.png", width = 500, height = 500, pointsize = 15)
  par(mfrow = c(1,1), mar = c(2,3,2,4), oma = c(4,4,4,4))
  heatmap(ProbPool)
  
  mtext("D",3, outer = T, adj = 0 , cex = 2)
  
  dev.off()
  
  # write the final image 
  
  magick::image_write(
    magick::image_append(c(
      magick::image_append(c(
        magick::image_read("Maps.png"),
        magick::image_read("PCAs.png"))),
      magick::image_append(c(
        magick::image_read("RDA.png"),
        magick::image_read("Pool.png")))),
      stack = T),
    "AppendixCaseStudy.png")
  
  ##########

  
  
  
  
  png("AssemblyEffects.png", width = 1000, height = 1000, pointsize = 20, res = 110)
  par(mfrow = c(1,1), mar = c(4,5,2,2),oma = c(2,2,2,2), las = 2)
  boxplot(forBoxplot$EF~forBoxplot$treat+forBoxplot$elev,
          width = rep(0.5,12),
          frame = F,
          xaxt = "n",
          ylim = c(0.01,120),
          xlab = "Elevation",
          cex.lab = 1.5,
          cex.axis = 0.8,
          ylab = "Magnitude of effects",
          col = rep(RColorBrewer::brewer.pal(4, "Set1"),3),
          log = "y")
  abline(v = c(4.5, 4.5+4))
  axis(1, c(2.5,6.5,10.5), labels = c(1000,2000,3000), line = 0.3, las = 1)
  
  legend("bottomleft", fill =  RColorBrewer::brewer.pal(4, "Set1"), legend = levels(forBoxplot$treat),bty = "n")
  
  
  dev.off()
  
  ###########
  ## Figure appendix, change in network interaction probabilities
  
 
  
  
  
  
  png("Cajanuma.png", width = 500,500,pointsize = 12)
  plotCont("Cajanuma", intRQL)
  dev.off()
  
  png("Bombuscaro.png", width = 500,500,pointsize = 12)
  plotCont("Bombuscaro", intRQL)
  dev.off()
  
  png("Copalinga.png", width = 500,500,pointsize = 12)
  plotCont("Copalinga", intRQL)
  dev.off()
  
  png("ECSF.png", width = 500,500,pointsize = 12)
  plotCont("ECSF", intRQL)
  dev.off()
  
  png("Finca.png", width = 500,500,pointsize = 12)
  plotCont("Finca", intRQL)
  dev.off()
  
  png("Bellavista.png", width = 500,500,pointsize = 12)
  plotCont("Bellavista", intRQL)
  dev.off()
  
  
  
  
  
  
  
  
  # write the final image 
  
  magick::image_write(
    magick::image_append(c(
      magick::image_append(c(
        magick::image_append(c(
          magick::image_read("Bellavista.png"),
          magick::image_read("Cajanuma.png"))))),
      magick::image_append(c(
        magick::image_append(c(
          magick::image_read("ECSF.png"),
          magick::image_read("Finca.png"))))),
      magick::image_append(c(
        magick::image_append(c(
          magick::image_read("Bombuscaro.png"),
          magick::image_read("Copalinga.png")))))),
      stack = T),
    
    "DensityInteractions.png")
  
  
  
  
  
  #### Figure appendix) Relative contribution of community and interaction assembly processes to define network structure across the elevational gradient
  #############
  
  
  png("RelativeEffectsProcess.png", width = 1000, height = 1000, pointsize = 20, res = 110)
  par(mar = c(0,0,0,0), oma= c(2,2,2,2))
  Ternary::TernaryPlot(alab = " Interaction assembly (%)",
                       blab = " Resource assembly (%) ",
                       clab = " Consumer assembly (%)",
                       lab.col = c('darkgreen', '', '#B22222'),
                       point = 'up', 
                       lab.cex = 1.5, 
                       grid.minor.lines = 0,
                       grid.lty = 'solid',
                       col = "white",
                       grid.col = scales::alpha(c('darkgreen',"#256B92", "#B22222"),0.5), 
                       axis.col = rgb(0.6, 0.6, 0.6), 
                       ticks.col = rgb(0.6, 0.6, 0.6),
                       axis.rotate = FALSE,
                       padding = 0.08)
  
  Ternary::AddToTernary(points,
                        col = "orange",
                        abs(SESnet_mod_norm)[,c(1:3)] * 100, 
                        pch = as.numeric(as.factor(SESnet_mod_norm$elev))+14)
  
  Ternary::AddToTernary(points,
                        col = "purple",
                        abs(SESnet_nes_norm)[,c(1:3)] * 100, 
                        pch = as.numeric(as.factor(SESnet_nes_norm$elev))+14)
  
  Ternary::AddToTernary(points,
                        col = "black",
                        abs(SESnet_mod_norm)[,c(1:3)] * 100, 
                        pch = as.numeric(as.factor(SESnet_mod_norm$elev))-1)
  
  Ternary::AddToTernary(points,
                        col = "black",
                        abs(SESnet_nes_norm)[,c(1:3)] * 100, 
                        pch = as.numeric(as.factor(SESnet_nes_norm$elev))-1)
  
  legend("topleft",
         title = "Elevation",
         bty = "n",
         legend = c(1000,2000,3000),
         pch = c(15,16,17))
  
  
  legend("topright",
         bty = "n",
         legend = c("Modularity","Nestedness"),
         col = c("orange","purple"), 
         pch = 18)
  
  dev.off()
  
  
  
  
  png("AssDistr.png", width = 1000, height = 1000, pointsize = 20, res = 110)
  par(mfrow = c(1,1), oma = c(2,2,2,2), mar = c(4,5,2,2), las = 1)
  stPSS <- c(PSdir$PSdir_process_sd/max(abs(PSdir$PSdir_process_sd)))
  plot(abs(stPSS)~SESnet_mod$elev,
       ylim = c(0,1), 
       frame = F, 
       pch = 16,
       cex= 1.5,
       cex.lab = 1.5,
       ylab = "Process strength asymmetry",
       xlab = "Elevation",
       col = ifelse(sign(stPSS)<0,
                    "#B22222",
                    ""))
  
  
  points(abs(stPSS)~SESnet_mod$elev,
         cex = 1.5)
  mtext("C",3, outer = T, adj = 0 , cex = 2)
  legend("topleft", 
         cex = 0.7,
         bty = "n",
         c("Frugivore driven", "Plant driven"),
         fill = c("#B22222", ""))
  
  dev.off() 
  
  
  
  
  
  ###### Assymetry changes with metric
  
png("ProRan.png", width = 1000, height = 1000, pointsize = 20, res = 110)
  
  par(mfrow = c(1,1), mar = c(4,6,2,2), oma = c(1,2,2,1))
  plot(-myEF_pro_ran_p,elev_vec1, frame= F,
       xlim = c(-3,3), pch = c(18,16,18,16,17,17)-1,
       ylim  = c(0,4000),
       xaxt = "n",
       cex.lab = 1.5,
       cex.axis = 1,
       xlab = "Strength of community assembly",
       ylab = "Elevation",
       col = scales::alpha("#B22222", 0.7),
       cex = 1.5)
  axis(1, c(-3:8), c(-3:8), cex = 1.4)
  segments(0,elev_vec1,-myEF_pro_ran_p,elev_vec1, col = "#B22222", lwd = 2)
  points(-myEF_pro_ran_p,elev_vec1,
         cex =1.5, 
         pch = c(18,16,18,16,17,17)-16)
  abline(v=0)
  abline(h = c(500,1500, 2500,3500), lty = 2)
  points(-myEF2_pro_ran_a,c(elev_vec2),
         ylim = c(-3,8), 
         pch = c(18,16,18,16,17,17)-1,
         xlim  = c(0,4000),
         col = scales::alpha("#256B92", 0.7),
         cex = 1.5)
  segments(0,elev_vec2,-myEF2_pro_ran_a,elev_vec2, col = "#256B92" ,lwd = 2)
  points(-myEF2_pro_ran_a,c(elev_vec2),
         cex = 1.5, 
         pch = c(18,16,18,16,17,17)-16)
  mtext("A",3, outer = T, adj = 0 , cex = 2)

  dev.off()
  
  png("RanRan.png", width = 1000, height = 1000, pointsize = 20, res = 110)
  par(mfrow = c(1,1), mar = c(4,6,2,2), oma = c(1,2,2,1))
  
  plot(-myEF_r_ran_p,elev_vec1, frame= F,
       xlim = c(-3,3), pch = c(18,16,18,16,17,17)-1,
       ylim  = c(0,4000),
       xaxt = "n",
       cex.lab = 1.5,
       cex.axis = 1,
       xlab = "Strength of community assembly",
       ylab = "Elevation",
       col = scales::alpha("#B22222", 0.7),
       cex = 1.5)
  axis(1, c(-3:8), c(-3:8), cex = 1.4)
  segments(0,elev_vec1,-myEF_r_ran_p,elev_vec1, col = "#B22222", lwd = 2)
  points(-myEF_r_ran_p,elev_vec1,
         cex =1.5, 
         pch = c(18,16,18,16,17,17)-16)
  abline(v=0)
  abline(h = c(500,1500, 2500,3500), lty = 2)
  points(-myEF2_r_ran_a,c(elev_vec2),
         ylim = c(-3,8), 
         pch = c(18,16,18,16,17,17)-1,
         xlim  = c(0,4000),
         col = scales::alpha("#256B92", 0.7),
         cex = 1.5)
  segments(0,elev_vec2,-myEF2_r_ran_a,elev_vec2, col = "#256B92" ,lwd = 2)
  points(-myEF2_r_ran_a,c(elev_vec2),
         cex = 1.5, 
         pch = c(18,16,18,16,17,17)-16)
  mtext("B",3, outer = T, adj = 0 , cex = 2)
  dev.off()
  
  
  png("Ransd.png", width = 1000, height = 1000, pointsize = 20, res = 110)
  par(mfrow = c(1,1), mar = c(4,6,2,2), oma = c(1,2,2,1))
  
  
  plot(-myEF_r_sd_p,elev_vec1, frame= F,
       xlim = c(-3,3), pch = c(18,16,18,16,17,17)-1,
       ylim  = c(0,4000),
       xaxt = "n",
       cex.lab = 1.5,
       cex.axis = 1,
       xlab = "Strength of community assembly",
       ylab = "Elevation",
       col = scales::alpha("#B22222", 0.7),
       cex = 1.5)
  axis(1, c(-3:8), c(-3:8), cex = 1.4)
  segments(0,elev_vec1,-myEF_r_sd_p,elev_vec1, col = "#B22222", lwd = 2)
  points(-myEF_r_sd_p,elev_vec1,
         cex =1.5, 
         pch = c(18,16,18,16,17,17)-16)
  abline(v=0)
  abline(h = c(500,1500, 2500,3500), lty = 2)
  points(-myEF2_r_sd_a,c(elev_vec2),
         ylim = c(-3,8), 
         pch = c(18,16,18,16,17,17)-1,
         xlim  = c(0,4000),
         col = scales::alpha("#256B92", 0.7),
         cex = 1.5)
  segments(0,elev_vec2,-myEF2_r_sd_a,elev_vec2, col = "#256B92" ,lwd = 2)
  points(-myEF2_r_sd_a,c(elev_vec2),
         cex = 1.5, 
         pch = c(18,16,18,16,17,17)-16)
  mtext("C",3, outer = T, adj = 0 , cex = 2)
dev.off()  


magick::image_write(
  magick::image_append(c(
    magick::image_read("ProRan.png"),
    magick::image_read("RanRan.png"),
    magick::image_read("Ransd.png"))),
  path = "AssChangeProc.png")
  
