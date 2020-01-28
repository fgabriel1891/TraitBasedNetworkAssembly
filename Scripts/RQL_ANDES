################################################################################
## Case study networks andes ecuador 
# Jorg Albretch + M.Schleuning ---> G.Muñoz + E. L. Marjakangas 
################################################################################
# function to calculate functional dispersion and compare it to random communities 
################################################################################
nullFun <- function(d, a, nRan = 1e3, nIter = 1e3, method = "independentswap") {
  obs <- fdisp(d = d, a = a)$FDis
  ran <- replicate(nRan, fdisp(d = d, a = randomizeMatrix(a, null.model = method, 
                                                          iterations = nIter))$FDis)
  fd <- cbind(obs, ran)
  z <- apply(fd, 1, FUN = function(x) {
    (x[1] - mean(x[-1])) / sd(x[-1])
  })
  p <- apply(fd, 1, FUN = function(x) {
    ifelse(x[1] > mean(x[-1]), sum(x[1] <= x[-1]) / (nRan + 1), sum(x[1] >= x[-1]) / (nRan + 1))
  })
  out <- data.frame(elevation = names(obs), fdis = obs, z = z, p = p)
  return(out)
}
################################################################################
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
pdf("output/figure_01.pdf", width = 0.394 * 12, height = 0.394 * 6)
m <- matrix(rep(c(1, 2), each = 4), 1, 8, byrow = TRUE)
m <- apply(m, 2, rep, each = 4)
m[16] <- 3
layout(m)
par(cex = 0.7, mar = c(0.5, 0.5, 0.5, 0.5) + 0.1, oma = c(0, 1, 1, 0))
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
RLQ$mQ
aDis <- dist(RLQ$mQ)
# Plants - extract normed species scores in trait space & create distance matrix
RLQ$mR
pDis <- dist(RLQ$mR)
################################################################################
# create presence-absence matrices for plants and animals on each elevation
################################################################################
pComm <- as.matrix(xtabs(~ paste(elevation, site) + plantCode, data = dataSet1$N) > 0)
pComm <- pComm[, rownames(RLQ$mR)]
aComm <- as.matrix(xtabs(~ paste(elevation, site) + animalCode, data = dataSet1$N) > 0)
aComm <- aComm[, rownames(RLQ$mQ)]
aFD <- nullFun(aDis, aComm)
pFD <- nullFun(pDis, pComm)
# ALTERNATIVELY distance matrices based on the raw trait data could be used
#aDis <- vegdist(sqrt(aTRLQ), method = "mahalanobis")
#pDis <- vegdist(sqrt(pTRLQ), method = "mahalanobis")
#aFD <- nullFun(aDis, aComm)
#pFD <- nullFun(pDis, pComm)
################################################################################
################################################################################
# Figure 2. Plot of plant and animal FD against elevation
################################################################################
aFD$at <- rep(1:3, each = 2) + c(0.05, 0.1)
aFD$trophic_level <- "animals"
pFD$at <- rep(1:3, each = 2) - c(0.05, 0.1)
pFD$trophic_level <- "plants"
fd <- rbind(aFD, pFD)
fd$elevation <- substr(fd$elevation, 1, 4)
fd$col <- rep(c("firebrick1", "deepskyblue1"), each = 6)
pdf("output/figure_02.pdf", width = 0.394 * 6, height = 0.394 * 6)
par(cex = 0.7, las = 1, mar = c(4.5, 4.5, 2.5, 0.5) + 0.1)
plot.new()
plot.window(xlim = range(fd$at), ylim = range(fd$z))
with(fd, points(at, z, pch = 21, col = col, bg = ifelse(p < 0.05, col, "white")))
abline(h = 0, lty = 3)
axis(1, at = 1:3, labels = unique(fd$elevation)); axis(2); box()
mtext(side = 1, line = 2.75, cex = 0.7, "Elevation (m a.s.l.)")
mtext(side = 2, line = 2.75, cex = 0.7, las = 0, "SES FDis")
legend("topleft", bty = "n", ncol = 1, inset = c(0, -0.3),
       text.col = unique(fd$col),
       legend = unique(fd$trophic_level), xpd = TRUE)
legend("topright", bty = "n", ncol = 1, inset = c(0, -0.3),
       pch = 21, pt.bg = c("white", "black"),
       legend = c(expression(italic(p) > 0.05), expression(italic(p) < 0.05)), xpd = TRUE)
dev.off()
################################################################################
# END OF SCRIPT
################################################################################
################################################################################
### Implementation of the null model approach 
## Gabriel Munoz
## Jan 2020 



## Detect wheter interactions are clearly separated into subpools 

dataSet1$N$intID <- paste0(dataSet1$N$plantCode,"/",dataSet1$N$animalCode)
dataSet1$N$siteID <- paste0(dataSet1$N$site,"/",dataSet1$N$elevation)

inMa <- as.data.frame.matrix(table(dataSet1$N$intID,dataSet1$N$siteID))
intMDS <- vegan::metaMDS(inMa, distance = "bray")

## Pools diverge from plants to animals to interactions 
heatmap(pComm*1)
heatmap(aComm*1)
heatmap(as.matrix(inMa))


## How interactions look into a trait space of RQL 
intRQL <- dplyr::bind_rows(RLQ$mQ,RLQ$mR)

rownames(intRQL) <- c(rownames(RLQ$mQ), rownames(RLQ$mR))

intRQL <- data.frame(dataSet1$N, "RLQan" = RLQ$mQ[match(dataSet1$N$animalCode,rownames(RLQ$mQ)),][,1])
intRQL <- data.frame(intRQL, "RLQpla" = RLQ$mQ[match(dataSet1$N$plantCode,rownames(RLQ$mR)),][,1])

###################### Figure 1. Trait degree plots #############################
layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE))

hist(intRQL$RLQan, main = "MetaNetwork",
     col=scales::alpha("firebrick1", 0.5), xlab = "RQL 1axis")
hist(intRQL$RLQpla, col = scales::alpha("deepskyblue1", 0.5), add =T)


plot(density(intRQL$RLQan))

hist(intRQL$RLQan[intRQL$elevation == 1000], main = "1000",
     col=scales::alpha("firebrick1", 0.5), xlab = "RQL 1axis")
hist(intRQL$RLQpla[intRQL$elevation == 1000], col = scales::alpha("deepskyblue1", 0.5), add =T)



hist(intRQL$RLQan[intRQL$elevation == 2000], main = "2000",
     col=scales::alpha("firebrick1", 0.5), xlab = "RQL 1axis")
hist(intRQL$RLQpla[intRQL$elevation == 2000], col = scales::alpha("deepskyblue1", 0.5), add =T)



hist(intRQL$RLQan[intRQL$elevation == 3000], main = "3000",
     col=scales::alpha("firebrick1", 0.5), xlab = "RQL 1axis")
hist(intRQL$RLQpla[intRQL$elevation == 3000], 
     col = scales::alpha("deepskyblue1", 0.5), add =T)

###################### Figure 1. Trait degree plots #############################

###################### Figure 2. Quantify variation in network structure #############################

## function to subset networks by plant species (from Munoz et al 2019 JBI)
subNetbySp = function(net1, net2, colToSample = colToSample){ 
  
  # Check richness of networks of species, and find the smaller one 
  toSample = min(length(unique(net1[,colToSample])),
                 length(unique(net2[,colToSample])))
  # subsample both networks 
  if (length(unique(net1[,colToSample])) == toSample){
    eqNetworks = list("net1" = droplevels(net1),
                      "net2" = droplevels(net2[net2[,colToSample] 
                                               %in% sample(unique(net2[,colToSample]),toSample),]))
    print("Network 2 was subsampled")
  } else {
    eqNetworks = list("net1" = droplevels(net1[net1[,colToSample] 
                                               %in% sample(unique(net1[,colToSample])
                                                           ,toSample),]),
                      "net2" = droplevels(net2))
    print("Network 1 was subsampled")
  }
  print(paste("both networks have now",toSample, 
              "species of", names(net1)[colToSample]))
  return(eqNetworks)
}


## function to calculate modularity and make null distributions 
makeNull <- function(matrix,nullSim){
  matrix <-matrix[as.logical(rowSums(matrix != 0)), as.logical(colSums(matrix != 0))]
  
  matrix <- as.matrix.data.frame(matrix)
  print(dim(matrix))
  mod <- bipartite::computeModules(matrix)
  vn <- vegan::nullmodel(matrix, "r1")
  vns <- simulate(vn, nullSim)
  if(dim(matrix)[1] < 3){
    print("Too small of a network")
    return(0)
    
  } else{ 
    null <- sapply(1:nullSim, function(x) bipartite::computeModules(vns[,,x])@likelihood)
    modES <- (mod@likelihood - mean(null))/sd(null)
    return(modES)
    
  }
}
  
## function to calculate nestedness and make null distributions 
  
makeNullNes <- function(matrix, nullSim){
  # make into binary network first
  matrix[matrix > 0] <- 1
  # remove 0s 
  matrix <-matrix[as.logical(rowSums(matrix != 0)), as.logical(colSums(matrix != 0))]
  matrix <- as.matrix.data.frame(matrix)
  nes <- vegan::nestednodf(matrix)$statistic["NODF"]
  vn <- vegan::nullmodel(matrix, "r1")
  vns <- simulate(vn, nullSim)
  null <- sapply(1:nullSim, function(x) vegan::nestednodf(vns[,,x])$statistic["NODF"])
  modES <- (nes - mean(null))/sd(null)
  return(modES)}





## 1) Quantify meta-network structure and the variation in network structure when subsetting for the same number of interactions/plants/birds

# meta network object
metaNetEl <- xtabs(frequency~plantCode + animalCode, dataSet1$N)
# regional network object
metaNetEl1 <- xtabs(frequency~plantCode + animalCode  + elevation, dataSet1$N)
# observed network object 
metaNetEl2 <- xtabs(frequency~plantCode + animalCode + site + elevation, dataSet1$N)


# Establish a probabilisitc network  


# meta network object
pmetaNetEl <- metaNetEl/max(metaNetEl)
pmetaEl <- reshape2::melt(pmetaNetEl)
pmetaEl <- pmetaEl[!pmetaEl$value == 0,]


# regional network object

lnet <- reshape2::melt(metaNetEl1[,,1]) # 1000 
lnet <- lnet[!lnet$value == 0,]
lnet$value <- lnet$value/max(lnet$value)

length(unique(lnet$plantCode))


# local network object

onet <-  reshape2::melt(metaNetEl2[,,,1])
onet <- onet[!onet$value == 0,]
onet$value <- onet$value/max(onet$value)

## 2) Quantify local-network structure and the variation in network structure based on the same size as the observed networks

# make subsets 
## meta to regional


subNet <- subNetbySp(pmetaEl, lnet,1)
## regional to local 
subNet2 <- subNetbySp(lnet, onet,1)


## create networks  meta-regional
netpm <- xtabs(value~plantCode + animalCode, subNet$net1)
netbl <- xtabs(value~plantCode + animalCode, subNet$net2)

## create networks regional-local
netpm1 <- xtabs(value~plantCode + animalCode, subNet2$net1)
netbl1 <- xtabs(value~plantCode + animalCode, subNet2$net2)


## perform null models meta-regional

modPM <- makeNull(netpm,3)
modPM2 <- makeNull(netbl,3)

## perform null models meta-regional regional-local

modPM.1 <- makeNull(netpm1,3)
modPM2.1 <- makeNull(netbl1,3)


modPM, modPM2




subNetbySp(pmetaEl, lnet,1)$net2 == subNetbySp(pmetaEl, lnet,1)$net2


unique(pmetaEl$plantCode)
