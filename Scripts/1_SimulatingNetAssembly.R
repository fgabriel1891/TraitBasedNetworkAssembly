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
library(yhat) # variance partitioning 
library(lme4)
library(sjPlot)
# library(magick) # to export figures 
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
#saveRDS(simulx, "simulx.rds")
#saveRDS(simul_0, "simul_0.rds.")
# (End of simulations)

#### End cluster
snowfall::sfStop()
################


############################################################
## Step 3) Extract info from simulation output and isolate the niche based effects.
############################################################

RES_DD$EA
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
RES_DD$ass <- log(c(1-RES_DD$sigmaB)/c(1-RES_DD$sigmaA))
RES_DD$ass <- RES_DD$ass/max(abs(RES_DD$ass))


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



a_eff <- rbind(plot_nl_a,plot_fl_a,plot_mm_a)
b_eff <- rbind(plot_nl_b,plot_fl_b,plot_mm_b)



butd_eff <- data.frame(
        "nodf_eff"= log(abs(b_eff$nodf)/ abs(a_eff$nodf)),
        "q_eff"= log(abs(b_eff$q)/ abs(a_eff$q)),
        "ass" = a_eff$pss,
        "intHyp" =  c(rep("NL",25), rep("FL",25), rep("MM",25))
        
)







RES_DD$SymTF <- ifelse(RES_DD$ass>0,1,0)
intAss_mod_Q <- glm(RES_DD$QZscorex~RES_DD$ass*RES_DD$intHyp)
intAss_mod_NODF <- glm(RES_DD$NODFx~RES_DD$ass*RES_DD$intHyp)

summary(intAss_mod_Q)
with(summary(intAss_mod_Q), 1 - deviance/null.deviance)
with(summary(intAss_mod_NODF), 1 - deviance/null.deviance)

library(rcompanion)

nagelkerke(intAss_mod_Q)
nagelkerke(intAss_mod_NODF)
sjPlot::tab_model(intAss_mod_Q)
sjPlot::tab_model(intAss_mod_NODF)



CCData_mod <- yhat::commonalityCoefficients(RES_DD, "QZscorex", list("ass", "EA", "intHyp"))
CCData_nes <- yhat::commonalityCoefficients(RES_DD, "NODFx", list("ass", "EA","intHyp"))






dim(RES_ass)/3
RES_ass <- RES_DD[!RES_DD$ass == 0,]
RES_ass0 <- RES_DD[RES_DD$ass == 0,]

dim(RES_ass0)/3








butd_eff$key <- paste(a_eff$size1,a_eff$size2, butd_eff$intHyp, sep = "_")
RES_ass$key <- paste(RES_ass$sigmaA, RES_ass$sigmaB, RES_ass$intHyp, sep = "_")

ass_null <- data.frame(aggregate(RES_ass0$QZscorex, list(RES_ass0$intHyp), mean),
           "sd" = aggregate(RES_ass0$QZscorex, list(RES_ass0$intHyp), sd)$x)


ass_null2 <- data.frame(aggregate(RES_ass0$QZscorex, list(RES_ass0$intHyp), mean),
                       "sd" = aggregate(RES_ass0$QZscorex, list(RES_ass0$intHyp), sd)$x)








assdev <- ifelse(RES_ass$intHyp == "NL",
                 (RES_ass$QZscorex-ass_null2$x[1])/ass_null2$sd[1],
                 ifelse(RES_ass$intHyp == "MM",
                        (RES_ass$QZscorex-ass_null2$x[2])/ass_null2$sd[2],
                        (RES_ass$QZscorex-ass_null2$x[3])/ass_null2$sd[3]))


assdev2 <- ifelse(RES_ass$intHyp == "NL",
                 (RES_ass$NODFx-ass_null$x[1])/ass_null$sd[1],
                 ifelse(RES_ass$intHyp == "MM",
                        (RES_ass$NODFx-ass_null$x[2])/ass_null$sd[2],
                        (RES_ass$NODFx-ass_null$x[3])/ass_null$sd[3]))









## Step 4) Variance partitioning 
#########



try<-RES_DD[complete.cases(RES_DD),]

plot(dist(scale(try$QZscoresd))~dist(try$ass))



# model_aeff_nes <- lme4::lmer((NODFx)~ass + (1+ass|intHyp), data = RES_ass)
# model_aeff_q <- lme4::lmer((QZscorex)~ass + (1+ass|intHyp), data = RES_ass)
# 
# 
# 
# slopes <- c(lme4::ranef(model_aeff_nes)$intHyp$ass, 
#             lme4::ranef(model_aeff_q)$intHyp$ass)
# 
# slopes2 <- c(lme4::fixef(lmm3)["ass"]+
#                      lme4::ranef(lmm3)$intHyp$ass, 
#              lme4::fixef(lmm4)["ass"]+
#                      lme4::ranef(lmm4)$intHyp$ass)
# 
# 
# sjPlot::tab_model(model_aeff_nes)
# 
# sjPlot::tab_model(model_aeff_q)
# 
# 
# 
# 
# lme4::ranef(model_aeff_nes)
# 
# lme4::ranef(model_aeff_q)
# 
# sjPlot::tab_model(model_aeff_nes)
# sjPlot::tab_model(model_aeff_q)



matt3 <- apply(abs(matt), 2, function(x) x/rowSums(abs(matt)))
matt3 <- matt3*100




matt2 <- as.matrix.data.frame(cbind(a_eff$q,b_eff$q, butd_eff$ass))
apply(apply(matt2, 2, abs), 2, max)
matt2[,1] <- matt2[,1]/apply(apply(matt2, 2, abs), 2, max)[1]
matt2[,2] <- matt2[,2]/apply(apply(matt2, 2, abs), 2, max)[2]

matt2[,1] <- (matt2[,1]*50)+50
matt2[,2] <- (matt2[,2]*50)+50
matt2[,3] <- (matt2[,3]*50)+50




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

png("SimulRes2.png", 
    width = 1000, 
    height = 1000, 
    pointsize = 20, 
    res = 110)
par(mfrow = c(1,1), oma = c(1,1,1,1), mar = c(4,5,3,2), las = 1)
plot(NODFx~QZscorex,
     frame = F,
     xlim = c(30, 150),
     ylim = c(0,200),
     pch = ifelse( RES_DD$intHyp == "NL", 23, 
                   ifelse( RES_DD$intHyp == "FL", 21,
                           22)),
     col = "grey80",
     bg = scales::alpha(ifelse( RES_DD$intHyp == "NL", "#7fc97f", 
                   ifelse( RES_DD$intHyp == "FL", "#beaed4",
                           "#fdc086")),0.7),
     xlab = "Modularity (Q Z-score)", 
     ylab = "Nestedness (NODF Z-score)", 
     cex.lab = 1.5,
     cex.axis = 1,
     cex = 1.5, 
     data = RES_DD[!RES_DD$EA == "ND",])


#### 
# neutral arrows

segments("x0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$QZscorex
         -RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$QZscoresd,
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$QZscorex
         +RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$QZscoresd,
         "y0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$NODFx, 
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$NODFx,
         lwd = 4, lty = 1, col = scales::alpha("#1b9e77",0.4))

segments("y0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$NODFx
         -RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("NL"),]$NODFsd,
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$NODFx
         +RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("NL"),]$NODFsd,
         "x0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$QZscorex, 
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$QZscorex, 
         lwd = 4, lty = 1,col = scales::alpha("#1b9e77",0.4))

segments("x0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$QZscorex
         -RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$QZscoresd,
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$QZscorex
         +RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$QZscoresd,
         "y0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$NODFx, 
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$NODFx,
         lwd = 4, lty = 1, col = scales::alpha("#d95f02",0.4))

segments("y0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$NODFx
         -RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("MM"),]$NODFsd,
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$NODFx
         +RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("MM"),]$NODFsd,
         "x0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$QZscorex, 
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$QZscorex, 
         lwd = 4, lty = 1,col = scales::alpha("#d95f02",0.8))

segments("x0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$QZscorex
         -RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$QZscoresd,
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$QZscorex
         +RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$QZscoresd,
         "y0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$NODFx, 
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$NODFx,
         lwd = 4, lty = 1, col = scales::alpha("#7570b3",0.8))

segments("y0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$NODFx
         -RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("FL"),]$NODFsd,
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$NODFx
         +RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("FL"),]$NODFsd,
         "x0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$QZscorex, 
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$QZscorex, 
         lwd = 4, lty = 1,col = scales::alpha("#7570b3",0.8))
###################

# 
# Plot_ConvexHull( RES_DD$QZscorex[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL")],
#                  RES_DD$NODFx[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL")], 
#                  scales::alpha("#7fc97f",0.5), lwd = 3, lty = 2)
# 
# Plot_ConvexHull( RES_DD$QZscorex[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM")],
#                  RES_DD$NODFx[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM")], 
#                  scales::alpha("#fdc086",1), lwd = 3, lty = 2)
# 
# Plot_ConvexHull( RES_DD$QZscorex[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL")],
#                  RES_DD$NODFx[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL")], lcolor = 
#                 scales::alpha("#beaed4",0.5), lwd = 3, lty = 2)


legend("bottomleft", 
       legend = c("Morphological matching", 
                  "Forbidden links",
                  "Stochastic interactions"),
       pch = c(22,21,23), 
       pt.bg = c("#d95f02","#7570b3","#1b9e77"),
       cex = 0.7,
       bty = "n",
       lwd = 0.5,
       col =  c("#d95f02","#7570b3","#1b9e77"))

# legend("bottomleft", 
#        title = "Community assembly",
#        legend = c("Environmental filtering", 
#                   "Neutral assembly"),
#       pch = c("+","+"),
#        cex = 0.7, 
#        bty = "n",
#        col = c("gray80","black"))


 dev.off()


#######
# Figure 3: Bottom-up and Top-down effects of resource and community assembly processes. 
########


png("BUTDNODF.png", width = 1000, height = 1000, pointsize = 20, res = 110)

 par(mfrow = c(1,1), oma = c(1,1,3,1), mar = c(4,5,3,2))
 plot((butd_eff$nodf_eff)[match(RES_ass$key,butd_eff$key)]~c(assdev2/max(abs(assdev2))), 
     col = scales::alpha(ifelse(butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "NL", "#FADE43", 
                                ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "FL", "#E0A738",
                                        "#9C413D")),0.5),
     xlab = "PSA effect on nestedness",
     ylab = "Network assembly mode",
     cex.lab = 1.5,
     frame =F,
     xlim = c(-1,1),
     ylim = c(-6,6),
     pch = ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "NL", 16, 
                   ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "FL", 15,
                           17)))
points(butd_eff$nodf_eff[match(RES_ass$key,butd_eff$key)]~c(assdev2/max(abs(assdev2))) , 
       col = scales::alpha(ifelse(butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "NL", "#FADE43", 
                                  ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "FL", "#E0A738",
                                          "#9C413D")),0.9),
       pch = ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "NL", 1, 
                     ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "FL", 0,
                             2)))
abline(h=0,v=0, lty = 2)

legend("topleft", "Top-down", bty = "n")
legend("bottomleft", "Bottom-up", bty = "n")

legend("bottomright", 
       title = "Interaction assembly type",
       legend = c("Morphological matching", 
                  "Forbidden links",
                  "Stochastic interactions"),
       pch = c(2,0,1)+15,
       cex = 0.7,
       bty = "n",
       col = c("#9C413D","#E0A738","#FADE43"))
mtext("C",3, outer = T, adj = 0 , cex = 2.5)

dev.off()


cor((butd_eff$q_eff)[match(RES_ass$key,butd_eff$key)],c(assdev/max(abs(assdev))),method = "spearman")
cor((butd_eff$nodf_eff)[match(RES_ass$key,butd_eff$key)],c(assdev2/max(abs(assdev2))), method ="spearman")
     

###
png("BUTD_MOD.png", width = 1000, height = 1000, pointsize = 20, res = 110)

par(mfrow = c(1,1), oma = c(1,1,3,1), mar = c(4,5,3,2))

plot((butd_eff$q_eff)[match(RES_ass$key,butd_eff$key)]~c(assdev/max(abs(assdev))), 
     col = scales::alpha(ifelse(butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "NL", "#FADE43", 
                                ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "FL", "#E0A738",
                                        "#9C413D")),0.5),
     xlab = "PSA effect on modularity",
     ylab = "Network assembly mode",
     cex.lab = 1.5,
     frame =F,
     xlim = c(-1,1),
     ylim = c(-6,6),
     pch = ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "NL", 16, 
                   ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "FL", 15,
                           17)))
points(butd_eff$q_eff[match(RES_ass$key,butd_eff$key)]~c(assdev/max(abs(assdev))) , 
       col = scales::alpha(ifelse(butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "NL", "#FADE43", 
                                  ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "FL", "#E0A738",
                                          "#9C413D")),0.9),
       pch = ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "NL", 1, 
                     ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "FL", 0,
                             2)))
abline(h=0,v=0, lty = 2)

legend("topleft", "Top-down", bty = "n")
legend("bottomleft", "Bottom-up", bty = "n")

legend("bottomright", 
       title = "Interaction assembly type",
       legend = c("Morphological matching", 
                  "Forbidden links",
                  "Stochastic interactions"),
       pch = c(2,0,1)+15,
       cex = 0.7,
       bty = "n",
       col = c("#9C413D","#E0A738","#FADE43"))

mtext("B",3, outer = T, adj = 0 , cex = 2.5)

dev.off()

#########




# write the final image 

magick::image_write(
        magick::image_append(c(
                                magick::image_read("SimulRes.png"),
                                magick::image_read("BUTD_MOD.png"),
                                magick::image_read("BUTDNODF.png")),
                             stack = F),
        "Figure3Final.png")


##########
# Appendix figure
######




par(mfrow = c(1,2))

plot(abs(a_eff$nodf)~abs(b_eff$nodf),
     xlim = c(0.001,9), 
     ylim = c(0.001,9), 
     log = "xy",
     main = "Nestedness (NODF)",
     xlab = as.expression(bquote('SES'['NET']*' '*.("Consumer"))),
     ylab =  as.expression(bquote('SES'['NET']*' '*.("Resource"))),
     cex = 1+butd_eff$ass,
     col = scales::alpha(ifelse( butd_eff$intHyp == "NL", "#FADE43", 
                                 ifelse( butd_eff$intHyp == "FL", "#E0A738",
                                         "#9C413D")),0.4),
     pch = ifelse( butd_eff$intHyp == "NL", 16, 
                   ifelse( butd_eff$intHyp == "FL", 15,
                           17)))


points(abs(a_eff$nodf)~abs(b_eff$nodf),
       cex = 1+butd_eff$ass,
       col = scales::alpha(ifelse(butd_eff$intHyp == "NL", "#FADE43", 
                                  ifelse( butd_eff$intHyp == "FL", "#E0A738",
                                          "#9C413D")),0.9),
       pch = ifelse( butd_eff$intHyp == "NL", 1, 
                     ifelse( butd_eff$intHyp == "FL", 0,
                             2)))
abline(a=0,b=1, lty = 2)








plot(abs(a_eff$q)~abs(b_eff$q),
     xlim = c(0.001,9), 
     ylim = c(0.001,9), 
     log = "xy",
     main = "Modularity (Q)",
     xlab = as.expression(bquote('SES'['NET']*' '*.("Consumer"))),
     ylab =  as.expression(bquote('SES'['NET']*' '*.("Resource"))),
     cex = 1+butd_eff$ass,
     col = scales::alpha(ifelse( butd_eff$intHyp == "NL", "#FADE43", 
                                 ifelse( butd_eff$intHyp == "FL", "#E0A738",
                                         "#9C413D")),0.4),
     pch = ifelse( butd_eff$intHyp == "NL", 16, 
                   ifelse( butd_eff$intHyp == "FL", 15,
                           17)))


points(abs(a_eff$q)~abs(b_eff$q),
       cex = 1+butd_eff$ass,
       col = scales::alpha(ifelse(butd_eff$intHyp == "NL", "#FADE43", 
                                  ifelse( butd_eff$intHyp == "FL", "#E0A738",
                                          "#9C413D")),0.9),
       pch = ifelse( butd_eff$intHyp == "NL", 1, 
                     ifelse( butd_eff$intHyp == "FL", 0,
                             2)))

abline(a=0,b=1, lty = 2)








par(mfrow = c(2,2), mar = c(4,4,2,2))

plot(RES_ass$QZscorex~RES_ass$ass, 
     col = scales::alpha(ifelse(RES_ass$intHyp == "NL", "#FADE43", 
                                ifelse( RES_ass$intHyp == "FL", "#E0A738",
                                        "#9C413D")),0.9),
     pch = 16,
     ylab = "Q (Z-score)", 
     xlab = "Process strength asymmetry",
     cex = 1.5,
     frame = F)
points(RES_ass0$QZscorex~rep(0,dim(RES_ass0)[1]), 
       col = scales::alpha(ifelse(RES_ass0$intHyp == "NL", "#FADE43", 
                                  ifelse( RES_ass0$intHyp == "FL", "#E0A738",
                                          "#9C413D")),0.9),
       cex=  1.5)
       
plot(RES_ass$NODFx~RES_ass$ass, 
     col = scales::alpha(ifelse(RES_ass$intHyp == "NL", "#FADE43", 
                                ifelse( RES_ass$intHyp == "FL", "#E0A738",
                                        "#9C413D")),0.9),
     pch = 16,
     ylab = "NODF (Z-score)", 
     xlab = "Process strength asymmetry",
     cex = 1.5,
     frame = F)
points(RES_ass0$NODFx~rep(0,dim(RES_ass0)[1]), 
       col = scales::alpha(ifelse(RES_ass0$intHyp == "NL", "#FADE43", 
                                  ifelse( RES_ass0$intHyp == "FL", "#E0A738",
                                          "#9C413D")),0.9),
       cex=  1.5)


plot(RES_ass$QZscoresd~RES_ass$ass, 
     col = scales::alpha(ifelse(RES_ass$intHyp == "NL", "#FADE43", 
                                ifelse( RES_ass$intHyp == "FL", "#E0A738",
                                        "#9C413D")),0.9),
     pch = 16,
     ylab = "Q (Z-score) SD", 
     xlab = "Process strength asymmetry",
     cex = 1.5,
     frame = F)
points(RES_ass0$QZscoresd~rep(0,dim(RES_ass0)[1]), 
       col = scales::alpha(ifelse(RES_ass0$intHyp == "NL", "#FADE43", 
                                  ifelse( RES_ass0$intHyp == "FL", "#E0A738",
                                          "#9C413D")),0.9),
       cex=  1.5)

plot(RES_ass$NODFsd~RES_ass$ass, 
     col = scales::alpha(ifelse(RES_ass$intHyp == "NL", "#FADE43", 
                                ifelse( RES_ass$intHyp == "FL", "#E0A738",
                                        "#9C413D")),0.9),
     pch = 16,
     ylab = "NODF (Z-score) SD", 
     xlab = "Process strength asymmetry",
     cex = 1.5,
     frame = F)
points(RES_ass0$NODFsd~rep(0,dim(RES_ass0)[1]), 
       col = scales::alpha(ifelse(RES_ass0$intHyp == "NL", "#FADE43", 
                                  ifelse( RES_ass0$intHyp == "FL", "#E0A738",
                                          "#9C413D")),0.9),
       cex=  1.5)





