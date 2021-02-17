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

###############
## Section 4: SIMULATING NETWORK ASSEMBLY
###############
## Load libraries 

library(ecolottery)
library(vegan)
library(bipartite)
library(tidyverse)
library(MASS)
library(corrplot)
library(foreach) # apply functions sequencially 
library(MASS) # stats
library(vegan) # eco stats 
library(yhat) # variance partitioning 

############################################################
## Step 1) define paramter objects and create combinations of assembly scenarios to simulate
############################################################

#### With similar processes at both trophic levels 

## Scenarios of environmental filtering + combinations of interaction assembly

Poss2 <- expand.grid("EA" = c("SEF"), # community assembly at trophic level A
                     "EB" = c("SEF"), # community assembly at trophic level B
                     "intHyp" = c("NL","FL", "MM"), # interaction assembly between A and B
                     "TraitA" = c(0.2), # trait optima  at trophic level A
                     "TraitB"= c(0.2), # trait optima  at trophic level B
                     "sigmaA" = seq(0.1,1,0.2), # strenght of filtering at trophic level A
                     "sigmaB" = seq(0.1,1,0.2),# strenght of filtering at trophic level B
                     "JA" = 500, # individuals at local community at trophic level A
                     "JB" = 500,  # individuals at local community at trophic level B
                     "mA" = 0.5,  # migration parameter (for neutral scenarios only, not used for niche-based assembly but needed in the object)
                     "mB" = 0.5    # migration parameter (for neutral scenarios only)
)

## Scenarios of neutral assembly  + combinations of interaction assembly 

Poss21 <- expand.grid("EA" = c("ND"),
                      "EB" = c("ND"), 
                      "intHyp" = c("NL", "FL", "MM"),
                      "TraitA" = c(0.2),
                      "TraitB"= c(0.2), 
                      "sigmaA" = (0), ## No filtering strength in neutral scenarios
                      "sigmaB" = (0), ## Idem 
                      "JA" = 500,
                      "JB" = 500,
                      "mA" = 0.5,
                      "mB" = 0.5
)


#### With contrasting processes at both trophic levels 


PossNET <- expand.grid("EA" = c("SEF"),
                       "EB" = c("ND"), 
                       "intHyp" = c("NL", "FL", "MM"),
                       "TraitA" = c(0.2),
                       "TraitB"= c(0), 
                       "sigmaA" = seq(0.1,1,0.2), ## No filtering strength in neutral scenarios
                       "sigmaB" = 0, ## Idem 
                       "JA" = 500,
                       "JB" = 500,
                       "mA" = 0.5,
                       "mB" = 0.5
)


PossNET2 <- expand.grid("EA" = c("ND"),
                        "EB" = c("SEF"), 
                        "intHyp" = c("NL", "FL", "MM"),
                        "TraitB" = c(0.2),
                        "TraitA"= c(0), 
                        "sigmaB" = seq(0.1,1,0.2), ## No filtering strength in neutral scenarios
                        "sigmaA" = 0, ## Idem 
                        "JA" = 500,
                        "JB" = 500,
                        "mA" = 0.5,
                        "mB" = 0.5
)



## Bind neutral and e.filtering scenarios into a single object
Poss <- rbind(Poss2, Poss21)
dim(Poss) # check the dimensions of the object  
PossNETall <- rbind(PossNET, PossNET2)
##############


############################################################
## Step 2) Run the simulations with the scenario object(s) created in step 1. 
############################################################
## Important: This is a time consuming step. Therefore, I strongly recommend running simulations in a cluster.
# However functions simulations can be parallelized at will. (See CustomFunctions.R)

## 2.1) Create species pool 
ss <- CreateSpPool(Jpool = c(5000, 5000),
                   J = c(500,500), 
                   poolShape = "uniform")

## hard save the species pool object.
## saveRDS(ss, "pool.rds")

# Important note if running in a cluster:
### Parallelization starts from here.
### set the numbers of cpus accordingly before running this section. 

cpus <- 80 

## Start parallelization (Make sure the snowfall library has been loaded)
library(snowfall)

sfInit(parallel=TRUE, cpus=cpus, type="SOCK")  # start cluster 
sfSource("Scripts/0_Functions.R") # Source functions to the cluster
# Export the scenario object(s) to the cluster (i.e. objects created in step 1)
sfExport("ss")
sfExport("PossNETall")

# Send packages to cluster 

sfLibrary(ecolottery)
sfLibrary(vegan)
sfLibrary(bipartite)
sfLibrary(tidyverse)
sfLibrary(MASS)
sfLibrary(corrplot)

### Running the simulations (time consuming step)

## Environmental filtering + neutral scenarios 

simulx <- RunSimul(Pool = ss,# species pool  
                   Scenario = Poss, # assembly scenario combinations (# modify this line to change the assembly scenario accordingly )
                   replicates = 10, # number of replicates of each assembly scenario
                   repNull = 100, # number of replicates of reshuffles to calculate Z-scores for network metrics 
                   quantile = 5,# quantile to consider "realized interactions" 
                   nulltype = "DD", # null model to choose. 
                   runParal = T) # Will it run in parallel?  (Strongly recomended)


simul_0 <- RunSimul(Pool = ss,# species pool  
                    Scenario = PossNETall, # assembly scenario combinations (# modify this line to change the assembly scenario accordingly )
                    replicates = 10, # number of replicates of each assembly scenario
                    repNull = 100, # number of replicates of reshuffles to calculate Z-scores for network metrics 
                    nulltype = "DD", # null model to choose. 
                    quantile = 5, # quantile to consider "realized interactions" 
                    runParal = T) # Will it run in parallel?  (Strongly recomended)

#### Hard save the object  (To ensure replicability)
# saveRDS(simulx, "simulx.rds")
# saveRDS(simul_0, "simul_0.rds.")
# (End of simulations)

#### End cluster
snowfall::sfStop()
################


############################################################
## Step 3) Extract info from simulation output and isolate the niche based effects.
############################################################


## Create an interpretable object  

## 78 equal scenarios
RES_DD <- Poss
RES_DD$QZscorex <- apply(sapply(1:3,function(y) sapply(1:78, function(x)simulx[[x]][,y]$metrics[[5]]$metrics$Q.Z.score)), 1, FUN = mean)
RES_DD$QZscoresd <- apply(sapply(1:3,function(y) sapply(1:78, function(x)simulx[[x]][,y]$metrics[[5]]$metrics$Q.Z.score)), 1, FUN = sd)
RES_DD$NODFx <- apply(sapply(1:3,function(y) sapply(1:78, function(x)simulx[[x]][,y]$metrics[[5]]$metrics$NODF.Z.score)), 1, FUN = mean)
RES_DD$NODFsd <- apply(sapply(1:3,function(y) sapply(1:78, function(x)simulx[[x]][,y]$metrics[[5]]$metrics$NODF.Z.score)), 1, FUN = sd)
RES_DD$con_x <- apply(sapply(1:3,function(y) sapply(1:78, function(x)simulx[[x]][,y]$metrics[[5]]$metrics$con)), 1, FUN = mean)
RES_DD$con_sd <- apply(sapply(1:3,function(y) sapply(1:78, function(x)simulx[[x]][,y]$metrics[[5]]$metrics$con)), 1, FUN = sd)


# 30 contrasting scenarios
RES_DD_2 <- PossNETall
RES_DD_2$QZscorex <- apply(sapply(1:3,function(y) sapply(1:30, function(x)simul_0[[x]][,y]$metrics[[5]]$metrics$Q.Z.score)), 1, FUN = mean)
RES_DD_2$QZscoresd <- apply(sapply(1:3,function(y) sapply(1:30, function(x)simul_0[[x]][,y]$metrics[[5]]$metrics$Q.Z.score)), 1, FUN = sd)
RES_DD_2$NODFx <- apply(sapply(1:3,function(y) sapply(1:30, function(x)simul_0[[x]][,y]$metrics[[5]]$metrics$NODF.Z.score)), 1, FUN = mean)
RES_DD_2$NODFsd <- apply(sapply(1:3,function(y) sapply(1:30, function(x)simul_0[[x]][,y]$metrics[[5]]$metrics$NODF.Z.score)), 1, FUN = sd)
RES_DD_2$con_x <- apply(sapply(1:3,function(y) sapply(1:30, function(x)simul_0[[x]][,y]$metrics[[5]]$metrics$con)), 1, FUN = mean)
RES_DD_2$con_sd <- apply(sapply(1:3,function(y) sapply(1:30, function(x)simul_0[[x]][,y]$metrics[[5]]$metrics$con)), 1, FUN = sd)

## calculate process strength symmetry and add it as a variable 
RES_DD$PSSdir <- log((RES_DD$sigmaA/RES_DD$sigmaB)^2)
RES_DD$PSSmag <- log(sqrt((RES_DD$sigmaA-RES_DD$sigmaB)^2)+1)
RES_DD$PSSdir <- RES_DD$PSSdir /max(RES_DD$PSSdir, na.rm = T)
RES_DD$PSSmag <- RES_DD$PSSmag /max(RES_DD$PSSmag, na.rm = T)

RES_DD_2$PSSdir <- log((RES_DD_2$sigmaA/RES_DD_2$sigmaB)^2)
RES_DD_2$PSSmag <- log(sqrt((RES_DD_2$sigmaA-RES_DD_2$sigmaB)^2)+1)
RES_DD_2$PSSdir <- RES_DD_2$PSSdir /max(RES_DD_2$PSSdir, na.rm = T)
RES_DD_2$PSSmag <- RES_DD_2$PSSmag /max(RES_DD_2$PSSmag, na.rm = T)

#########
## Step 4) Variance partitioning 
#########

RES_DD$ass <- log(RES_DD$sigmaA/RES_DD$sigmaB)

lmm <- lme4::lmer(QZscorex~ass + (1+ass|intHyp), RES_DD[complete.cases(RES_DD),])
lmm2 <- lme4::lmer(NODFx~ass + (1+ass|intHyp), RES_DD[complete.cases(RES_DD),])

lmm3 <- lme4::lmer(QZscoresd~ass + (1+ass|intHyp), RES_DD[complete.cases(RES_DD),])
lmm4 <- lme4::lmer(NODFsd~ass + (1+ass|intHyp), RES_DD[complete.cases(RES_DD),])



slopes <- c(lme4::fixef(lmm)["ass"]+lme4::ranef(lmm)$intHyp$ass, lme4::fixef(lmm2)["ass"]+lme4::ranef(lmm2)$intHyp$ass)
slopes2 <- c(lme4::fixef(lmm3)["ass"]+lme4::ranef(lmm3)$intHyp$ass, lme4::fixef(lmm4)["ass"]+lme4::ranef(lmm4)$intHyp$ass)









############################################################
## Step 5) Calculate Bottom-Up and Top-Down effects  
############################################################

cr_nl_a <- sapply(seq(0.1,1, 0.2), function(x) con_res_efz(RES_DD_2, RES_DD, x, "NL","a"))
cr_fl_a <- sapply(seq(0.1,1, 0.2), function(x) con_res_efz(RES_DD_2, RES_DD, x, "FL", "a"))
cr_mm_a <- sapply(seq(0.1,1, 0.2), function(x) con_res_efz(RES_DD_2, RES_DD, x, "MM","a"))


cr_nl_b <- sapply(seq(0.1,1, 0.2), function(x) con_res_efz(RES_DD_2, RES_DD, x, "NL","b"))
cr_fl_b <- sapply(seq(0.1,1, 0.2), function(x) con_res_efz(RES_DD_2, RES_DD, x, "FL", "b"))
cr_mm_b <- sapply(seq(0.1,1, 0.2), function(x) con_res_efz(RES_DD_2, RES_DD, x, "MM","b"))

## Rearrange dataframe 

plot_nl_a <- toplot(cr_nl_a)
plot_fl_a <- toplot(cr_fl_a)
plot_mm_a <- toplot(cr_mm_a)

plot_nl_b <- toplot(cr_nl_b)
plot_fl_b <- toplot(cr_fl_b)
plot_mm_b <- toplot(cr_mm_b)




############################################################
## Step 5) Plot the results
############################################################


##############
## Figure 3. Effects of community and interaction assembly processes on the structure of simulated resource-consumer networks
##############

########
## Figure 3a: Two dimensional space defined by modularity and nestedness axis.
#####
#######
dev.off()
par( mfrow = c(1,2),mar=c(6,6,4,4), las = 1)
plot(NODFx~QZscorex,
     frame = F,
     xlim = c(0, 150),
     ylim = c(0,250),
     pch = ifelse( RES_DD$intHyp == "NL", 1, 
                   ifelse( RES_DD$intHyp == "FL", 0,
                           2)),
     col = ifelse( RES_DD$intHyp == "NL", "#FADE43", 
                   ifelse( RES_DD$intHyp == "FL", "#E0A738",
                           "#9C413D")),
     xlab = "Modularity (Q-Zscore)", 
     ylab = "Nestedness (NODF-Zscore)", 
     cex.lab = 1.5,
     cex.axis = 1,
     cex = 1, 
     data = RES_DD[!RES_DD$EA == "ND",])


Plot_ConvexHull( RES_DD$QZscorex[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL")],
                 RES_DD$NODFx[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL")], 
                 scales::alpha("#FADE43",0.5), lwd = 3, lty = 2)

Plot_ConvexHull( RES_DD$QZscorex[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM")],
                 RES_DD$NODFx[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM")], 
                 scales::alpha("#9C413D",1), lwd = 3, lty = 2)

Plot_ConvexHull( RES_DD$QZscorex[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL")],
                 RES_DD$NODFx[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL")], lcolor = 
                scales::alpha("#E0A738",0.5), lwd = 3, lty = 2)

#### 
# neutral arrows

segments("x0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$QZscorex
         -RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$QZscoresd,
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$QZscorex
         +RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$QZscoresd,
         "y0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$NODFx, 
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$NODFx,
         lwd = 2, lty = 1, col = scales::alpha("#FADE43",0.9))

segments("y0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$NODFx
         -RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("NL"),]$NODFsd,
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$NODFx
         +RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("NL"),]$NODFsd,
         "x0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$QZscorex, 
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$QZscorex, 
         lwd = 2, lty = 1,col = scales::alpha("#FADE43",0.9))

segments("x0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$QZscorex
         -RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$QZscoresd,
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$QZscorex
         +RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$QZscoresd,
         "y0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$NODFx, 
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$NODFx,
         lwd = 2, lty = 1, col = scales::alpha("#9C413D",0.9))

segments("y0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$NODFx
         -RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("MM"),]$NODFsd,
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$NODFx
         +RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("MM"),]$NODFsd,
         "x0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$QZscorex, 
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$QZscorex, 
         lwd = 2, lty = 1,col = scales::alpha("#9C413D",0.9))

segments("x0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$QZscorex
         -RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$QZscoresd,
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$QZscorex
         +RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$QZscoresd,
         "y0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$NODFx, 
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$NODFx,
         lwd = 2, lty = 1, col = scales::alpha("#E0A738",0.9))

segments("y0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$NODFx
         -RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("FL"),]$NODFsd,
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$NODFx
         +RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("FL"),]$NODFsd,
         "x0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$QZscorex, 
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$QZscorex, 
         lwd = 2, lty = 1,col = scales::alpha("#E0A738",0.9))
###################


legend("topleft", 
       title = "Interaction assembly",
       legend = c("Morphological matching", 
                  "Forbidden links",
                  "Stochastic interactions"),
       pch = c(2,0,1), cex = 0.7,
       col = c("#9C413D","#E0A738","#FADE43"))

legend("bottomleft", 
       title = "Community assembly",
       legend = c("Environmental filtering", 
                  "Neutral assembly"),
       lty = c(2,1),cex = 0.7,
       col = c("gray80","black"))











barplot(slopes[c(1,4,2,5,3,6)],
        ylim = c(-5,15), 
        col = c("#FADE43","#FADE43","#E0A738",
                "#E0A738","#9C413D","#9C413D"),
        ylab = "ß Process strength assymetry", 
        cex.lab = 1.5)
legend("topright", 
       legend = c("Morphological matching",
                  "Forbidden links",
                  "Stochastic interactions"),
       cex = 0.7,
       fill = c("#FADE43","#E0A738","#9C413D"))

#################




barplot(slopes2[c(1,4,2,5,3,6)],
        ylim = c(0,1), 
        col = c("#FADE43","#FADE43","#E0A738","#E0A738","#9C413D","#9C413D"))








#######
# Figure 3: Bottom-up and Top-down effects of resource and community assembly processes. 
########


plot(log(abs(plot_nl_a$nodf)/abs(plot_nl_b$nodf))~
       log(plot_nl_a$size1/plot_nl_b$size2),
     xlim = c(-3,3),
     xlab = "Process strenght assymetry", 
     ylab = "Effect size log-ratio \n (SES consumer/SES resource)",
     col = "#FADE43",
     ylim = c(-3,3))
abline(h = 0, v = 0, lty = 2)
abline(lm(log(abs(plot_nl_a$nodf)/abs(plot_nl_b$nodf))~
            log(plot_nl_a$size1/plot_nl_b$size2)),
       col = "#FADE43")

points(log(abs(plot_fl_a$nodf)/abs(plot_fl_b$nodf))~
       log(plot_fl_a$size1/plot_fl_b$size2), 
       pch = 2, col = "#E0A738")

abline(lm(log(abs(plot_fl_a$nodf)/abs(plot_fl_b$nodf))~
            log(plot_fl_a$size1/plot_fl_b$size2)),
       col =  "#E0A738")

points(log(abs(plot_mm_a$nodf)/abs(plot_mm_b$nodf))~
       log(plot_mm_a$size1/plot_mm_b$size2), 
       pch = 1, col = "#9C413D")

abline(lm(log(abs(plot_mm_a$nodf)/abs(plot_mm_b$nodf))~
            log(plot_mm_a$size1/plot_mm_b$size2)),
       pch = 1,  col = "#9C413D")


###

plot(log(abs(plot_nl_a$q)/abs(plot_nl_b$q))~
       log(plot_nl_a$size1/plot_nl_b$size2),
     xlim = c(-3,3),
     col = "#FADE43",
     xlab = "Process strenght assymetry", 
     ylab = "Effect size log-ratio \n (SES consumer/SES resource)",
     ylim = c(-3,3))
abline(h = 0, v = 0, lty = 2)
abline(lm(log(abs(plot_nl_a$q)/abs(plot_nl_b$q))~
            log(plot_nl_a$size1/plot_nl_b$size2)),
       col = "#FADE43")

points(log(abs(plot_fl_a$q)/abs(plot_fl_b$q))~
         log(plot_fl_a$size1/plot_fl_b$size2), 
       pch = 2, col = "#E0A738")

abline(lm(log(abs(plot_fl_a$q)/abs(plot_fl_b$q))~
            log(plot_fl_a$size1/plot_fl_b$size2)),
       col =  "#E0A738")

points(log(abs(plot_mm_a$q)/abs(plot_mm_b$q))~
         log(plot_mm_a$size1/plot_mm_b$size2), 
       pch = 1, col = "#9C413D")

abline(lm(log(abs(plot_mm_a$q)/abs(plot_mm_b$q))~
            log(plot_mm_a$size1/plot_mm_b$size2)),
       pch = 1,  col = "#9C413D")











