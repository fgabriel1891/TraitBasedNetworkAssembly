### Synthetized script to run trait-based simulations of network structure 
## Gabriel Muñoz
## gabriel.munoz@concordia.ca
## Jan 2020 

## Load custom functions (important!)
source("scripts/CustomFunctions.R")


## Note: Go trough the CustomFunctions.R script for details on the functions and the parameters 

############################################################
############################################################
## 1) define paramter objects and create combinations of assembly scenarios to simulate
############################################################
############################################################

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



## Bind neutral and e.filtering scenarios into a single object
Poss <- rbind(Poss2, Poss21)
dim(Poss) # check the dimensions of the object  

## Important!: Uncomment the lines below if you want to implement variable trait optimas but similar values between trophic levels (however this wont affect drastically 
# the results in network metrics and deviations, so far we leave the trait optima as fixed to ease the interpretation of results)

# Poss$TraitA =  rnorm(78, 0.5,0.2)
# Poss$TraitB =  Poss$TraitA 
# 


## Scenarios of limiting similarity + combinations of interaction assembly


Poss23 <- expand.grid("EA" = c("LS"),
                      "EB" = c("LS"), 
                      "intHyp" = c("FL", "MM", "NL"),
                      "TraitA" = c(0.2),
                      "TraitB"= c(0.2), 
                      "sigmaA" = (0), ## No filtering strength in neutral scenarios
                      "sigmaB" = (0), ## Idem
                      "chunkA" = c(2:6), # Strenght of limiting similarity (see CustomFunctions.R)
                      "chunkB" = c(2:6),
                      "JA" = 500,
                      "JB" = 500,
                      "mA" = 0.5,
                      "mB" = 0.5
)

dim(Poss23) # Check the dimensions of the resulting object  

## Important!: Uncomment the lines below if you want to implement variable trait optimas but similar values between trophic levels (however this wont affect drastically 
# the results in network metrics and deviations, so far we leave the trait optima as fixed to ease the interpretation of results)
# Poss23$TraitA =  rnorm(75, 0.5,0.2)
# Poss23$TraitB =  Poss23$TraitA 




######
# Simulations with varying trait optima between trophic levels ( varying trait optima with differences between trophic levels)

#####

### Environmental filtering

SimVarTrait <- expand.grid("EA" = c("SEF"),
                           "EB" = c("SEF"),
                           "intHyp" = c("NL", "FL", "MM"),
                           "sigmaA" = seq(0.1,1,0.30),
                           "sigmaB" = seq(0.1,1,0.30),
                           "JA" = 500,
                           "TraitA" = 0.2,
                           "TraitB" = 0.2,
                           "JB" = 500,
                           "mA" = 0.5,
                           "mB" = 0.5
)

# SimVarTrait$TraitA =  rnorm(dim(SimVarTrait)[1], 0.5,0.2)
# SimVarTrait$TraitB =  rnorm(dim(SimVarTrait)[1], 0.5,0.2)
# SimVarTrait$traitDif = abs(SimVarTrait$TraitA-SimVarTrait$TraitA)

### Neutral assemblages

SimVarTrait_NL <- expand.grid("EA" = c("ND"),
                              "EB" = c("ND"),
                              "intHyp" = c("NL","FL", "MM"),
                              "sigmaA" = 0,
                              "sigmaB" = 0,
                              "TraitA" = 0.2,
                              "TraitB" = 0.2,
                              "JA" = 500,
                              "JB" = 500,
                              "mA" = 0.5,
                              "mB" = 0.5
)
dim(SimVarTrait_NL)

# SimVarTrait_NL$TraitA =  rnorm(dim(SimVarTrait_NL)[1], 0.5,0.2)
# SimVarTrait_NL$TraitB =  rnorm(dim(SimVarTrait_NL)[1], 0.5,0.2)
# SimVarTrait_NL$traitDif = abs(SimVarTrait_NL$TraitA-SimVarTrait_NL$TraitB)


### Limiting similarity

SimVarTrait_LS <- expand.grid("EA" = c("LS"),
                              "EB" = c("LS"),
                              "intHyp" = c("NL","FL", "MM"),
                              "chunkA" = c(2:4),
                              "chunkB" = c(2:4),
                              "TraitA" = 0.2,
                              "TraitB" = 0.2,
                              "JA" = 500,
                              "JB" = 500,
                              "mA" = 0.5,
                              "mB" = 0.5
)

#
# SimVarTrait_LS$TraitA =  rnorm(dim(SimVarTrait_LS)[1], 0.5,0.2)
# SimVarTrait_LS$TraitB =  rnorm(dim(SimVarTrait_LS)[1], 0.5,0.2)
# SimVarTrait_LS$traitDif = abs(SimVarTrait_LS$TraitA-SimVarTrait_LS$TraitB)
#

############################################################
############################################################
## 2) Run the simulations with the scenario object(s) created in step 1. 
############################################################
############################################################



## Important: This is a time consuming step. Therefore, I strongly recommend running simulations in a cluster.
# However functions simulations can be parallelized at will. (See CustomFunctions.R)

## 2.1) Create species pool 
ss = CreateSpPool(c(5000, 5000), c(500,500), "uniform")

## hard save the species pool object (for later reproducibility )
saveRDS(ss, "pool.rds")

# Important note if running in a cluster:
### Parallelization starts from here. Thus, set the numbers of cpus accordingly before running this section. 
## Start parallelization (Make sure the snowfall library has been loaded)

sfInit(parallel=TRUE, cpus=50, type="SOCK")  # start cluster 
sfSource("scripts/CustomFunctions.R") # Source functions to the cluster

# Export the scenario object(s) to the cluster (i.e. objects created in step 1)
sfExport("ss")

# Fixed trait optima
sfExport("Poss")
sfExport("Poss23")

# Varying trait optima 
sfExport("SimVarTrait")
sfExport("SimVarTrait_LS")
sfExport("SimVarTrait_NL")


### Running the simulations (time consuming step)

## Environmental filtering + neutral scenarios 

simul <- RunSimul(Pool = ss,# species pool
                  Scenario = Poss, # assembly scenario combinations (# modify this line to change the assembly scenario accordingly )
                  replicates = 10, # number of replicates of each assembly scenario
                  repNull = 100, # number of replicates of reshuffles to calculate Z-scores for network metrics 
                  quantile = 5, # quantile to consider "realized interactions" 
                  runParal = T) # Will it run in parallel?  (Strongly recomended)

## Limiting similarity scenarios 
limSim <- RunSimul(Pool = ss,
                   Scenario = Poss23,
                   replicates = 10,
                   repNull = 100,
                   quantile = 5,
                   runParal = T)



#### Hard save the object  (To ensure replicability)
saveRDS(simul, "EF_NL_simul_fixAug.rds")
saveRDS(limSim, "LS_simul_fixAug.rds.")
# (End of simulations)

################################################
# This section can simulate bipartite communtities but not calculate network metrics (#Useful to rapidly create bipartite networks based on varying assembly scenarios)
################################################
# Note we are not using this for the paper, but it could be useful for testing the null model approach 

## Environmental Filtering
traitNet <- GenBiNet(Pool= ss, 
                     SimVarTrait,
                     replicates = 10,
                     repNull = 1,
                     quantile = 5, 
                     runParal = T)
## Neutral Assemblages
traitNet_NL <- GenBiNet(Pool= ss,
                        SimVarTrait_NL,
                        replicates = 10,
                        repNull = 1,
                        quantile = 5, 
                        runParal = T)
## Limiting Similarity 
traitNet_LS <- GenBiNet(Pool= ss, 
                        SimVarTrait_LS,replicates = 10,
                        repNull = 1, 
                        quantile = 5,
                        runParal = T)
## Hard save the object 
# 
# saveRDS(traitNet, "traitNet.rds")
# saveRDS(traitNet_NL, "traitNet_NL.rds")
# saveRDS(traitNet_LS, "traitNet_LS.rds")

################################################
# End of section to simulate bipartite communtities but not calculate network metrics (#Useful to rapidly create bipartite networks based on varying assembly scenarios)
################################################

#### End cluster

snowfall::sfStop()



############################################################
############################################################
## 3) Interpret and the raw simulation output and isolate the niche based effects.
############################################################
############################################################

## dependencies 

library(foreach) # apply functions sequencially 
library(MASS) # stats
library(vegan) # eco stats 
library(yhat) # variance partitioning 


## Create an interpretable object  

## Environmental filtering + neutral assembly
RES1 <- sumarizeResults(Poss, 
                        simul, 
                        quantile = 5, 
                        replicate = 10)
## Limiting similarity 
RES2 <- sumarizeResults(Poss23, 
                        limSim, 
                        quantile = 5, 
                        replicate = 10)


## Isolate niche based effects 

#####
# Neutral trophic assemblages 
### 

# Neutral interactions
# NODF
NODF_NN_NL_di = unlist(RES1[RES1$EA == "ND" & RES1$intHyp =="NL",]["NODFx"]) # NODF when netrl assembled and int neutral
NODF_NN_NL_sd = unlist(RES1[RES1$EA == "ND" & RES1$intHyp =="NL",]["NODFsd"]) #  NODF Variation 
# Q
Q_NN_NL_di = unlist(RES1[RES1$EA == "ND" & RES1$intHyp =="NL",]["QZscorex"]) # Q when netrl assembled and int neutral
Q_NN_NL_sd = unlist(RES1[RES1$EA == "ND" & RES1$intHyp =="NL",]["QZscoresd"])  #  Q Variation 
### Forbidden Links
# NODF
NODF_NN_FL_di = unlist(RES1[RES1$EA == "ND" & RES1$intHyp =="FL",]["NODFx"]) # NODF when netrl assembled and int match
NODF_NN_FL_sd = unlist(RES1[RES1$EA == "ND" & RES1$intHyp =="FL",]["NODFsd"]) #  NODF Variation 
# Q
Q_NN_FL_di = unlist(RES1[RES1$EA == "ND" & RES1$intHyp =="FL",]["QZscorex"]) # Q when netrl assembled and int match
Q_NN_FL_sd = unlist(RES1[RES1$EA == "ND" & RES1$intHyp =="FL",]["QZscoresd"])  #  Q Variation 
######
### Morphological matching
# NODF
NODF_NN_MM_di = unlist(RES1[RES1$EA == "ND" & RES1$intHyp =="MM",]["NODFx"]) # NODF when netrl assembled and int match
NODF_NN_MM_sd = unlist(RES1[RES1$EA == "ND" & RES1$intHyp =="MM",]["NODFsd"]) #  NODF Variation 
# Q
Q_NN_MM_di = unlist(RES1[RES1$EA == "ND" & RES1$intHyp =="MM",]["QZscorex"]) # Q when netrl assembled and int match
Q_NN_MM_sd = unlist(RES1[RES1$EA == "ND" & RES1$intHyp =="MM",]["QZscoresd"])  #  Q Variation 
######




######
## Calculate effect sizes relative to neutral assemblages 
####


# SES variation from null models 

### Forbidden Links

NODF_SESNull_FL <- (NODF_NN_FL_di-NODF_NN_NL_di)/NODF_NN_NL_sd
Q_SESNull_FL <- (Q_NN_FL_di-Q_NN_NL_di)/Q_NN_NL_sd

### Morphological match

NODF_SESNull_MM <- (NODF_NN_MM_di-NODF_NN_NL_di)/NODF_NN_NL_sd
Q_SESNull_MM <- (Q_NN_MM_di-Q_NN_NL_di)/Q_NN_NL_sd



######
## niche based community assembly effects relative to neutral assemblages 
####
# Environmental filtering 
# Neutral interactions
# Nestedness
NODF_EF_NL_di = makeMesh(RES1, "NL" , "NODFx", scale = F) # NODF when EF and neutral interactions
NODF_EF_NL_w = (NODF_EF_NL_di-NODF_NN_NL_di)/NODF_NN_NL_sd # effect size in relation to the variation of expected under neutral assemblages and neutral interactions

# Modularity
Q_EF_NL_di = makeMesh(RES1, "NL" , "QZscorex", scale = F) 
Q_EF_NL_w <-  (Q_EF_NL_di-Q_NN_NL_di)/Q_NN_NL_sd # effect size in relation to the variation of expected under neutral assemblages and neutral interactions

# Forbidden Links
# Nestedness
NODF_EF_FL_di = makeMesh(RES1, "FL" , "NODFx", scale = F) # env. filtering and functional matching
NODF_EF_FL_w <-  (NODF_EF_FL_di-NODF_NN_FL_di)/NODF_NN_FL_sd # effect size in relation to the variation of expected under neutral assemblages and neutral interactions

# Modularity 
Q_EF_FL_di = makeMesh(RES1, "FL" , "QZscorex", scale = F) 
Q_EF_FL_w <-  (Q_EF_FL_di-Q_NN_FL_di)/Q_NN_FL_sd # effect size in relation to the variation of expected under neutral assemblages and neutral interactions


# Morphological match
# Nestedness
NODF_EF_MM_di = makeMesh(RES1, "MM" , "NODFx", scale = F) # env. filtering and functional matching
NODF_EF_MM_w <-  (NODF_EF_MM_di-NODF_NN_MM_di)/NODF_NN_MM_sd # effect size in relation to the variation of expected under neutral assemblages and neutral interactions

# Modularity 
Q_EF_MM_di = makeMesh(RES1, "MM" , "QZscorex", scale = F) 
Q_EF_MM_w <-  (Q_EF_MM_di-Q_NN_MM_di)/Q_NN_MM_sd # effect size in relation to the variation of expected under neutral assemblages and neutral interactions

# Limiting similarity  
# Neutral interactions
# Nestedness
NODF_LS_NL_di =  makeMesh2(RES2, "NL", "NODFx", scale = F)  # NODF when LS and neutral interactions
NODF_LS_NL_w = (NODF_LS_NL_di-NODF_NN_NL_di)/NODF_NN_NL_sd # effect size in relation to the variation of expected under neutral assemblages and neutral interactions

# Modularity
Q_LS_NL_di =  makeMesh2(RES2, "MM", "QZscorex", scale = F)
Q_LS_NL_w = (Q_LS_NL_di-Q_NN_NL_di)/Q_NN_NL_sd # effect size in relation to the variation of expected under neutral assemblages and neutral interactions

# Forbidden links 
# Nestedness
NODF_LS_FL_di = makeMesh2(RES2, "FL" , "NODFx", scale = F) # env. filtering and functional matching
NODF_LS_FL_w <-  (NODF_LS_FL_di-NODF_NN_FL_di)/NODF_NN_FL_sd # effect size in relation to the variation of expected under neutral assemblages and neutral interactions
# Modularity
Q_LS_FL_di = makeMesh2(RES2, "MM" , "QZscorex", scale = F) # env. filtering and functional matching
Q_LS_FL_w <-  (Q_LS_MM_di-Q_NN_MM_di)/Q_NN_MM_sd # effect size in relation to the variation of expected under neutral assemblages and neutral interactions

#Morphological Match
# Nestedness
NODF_LS_MM_di = makeMesh2(RES2, "MM" , "NODFx", scale = F) # env. filtering and functional matching
NODF_LS_MM_w <-  (NODF_LS_MM_di-NODF_NN_MM_di)/NODF_NN_MM_sd # effect size in relation to the variation of expected under neutral assemblages and neutral interactions
# Modularity
Q_LS_MM_di = makeMesh2(RES2, "MM" , "QZscorex", scale = F) # env. filtering and functional matching
Q_LS_MM_w <-  (Q_LS_MM_di-Q_NN_MM_di)/Q_NN_MM_sd # effect size in relation to the variation of expected under neutral assemblages and neutral interactions
#######

######
## niche based interaction assembly effects relative to neutral assemblages 
####


NODF_NN_b <-  (c(NODF_NN_MM_di)-c(NODF_NN_NL_di)/NODF_NN_NL_sd)
Q_NN_b <- (c(Q_NN_MM_di)-c(Q_NN_NL_di)/Q_NN_NL_sd)

NODF_NN_b2 <-  (c(NODF_NN_FL_di)-c(NODF_NN_NL_di)/NODF_NN_NL_sd)
Q_NN_b2 <- (c(Q_NN_FL_di)-c(Q_NN_NL_di)/Q_NN_NL_sd)


# Morphological match
# Environmental filtering 
# Nestedness

EF_NODF_MM_b <- c(NODF_EF_MM_di)-c(NODF_EF_NL_di)/sd(NODF_EF_NL_di)
EF_Q_MM_b <- c(Q_EF_MM_di)-c(Q_EF_NL_di)/sd(Q_EF_NL_di)


# Limiting similairty 
LS_NODF_MM_b <- c(NODF_LS_MM_di)-c(NODF_LS_NL_di)/sd(NODF_LS_NL_di)
LS_Q_MM_b <- c(Q_LS_MM_di)-c(Q_LS_NL_di)/sd(Q_LS_NL_di)

# Forbidden Links
# Environmental filtering 
# Nestedness

EF_NODF_FL_b <- c(NODF_EF_FL_di)-c(NODF_EF_NL_di)/sd(NODF_EF_NL_di)
EF_Q_FL_b <- c(Q_EF_FL_di)-c(Q_EF_NL_di)/sd(Q_EF_NL_di)


# Limiting similairty 
LS_NODF_FL_b <- c(NODF_LS_FL_di)-c(NODF_LS_NL_di)/sd(NODF_LS_NL_di)
LS_Q_FL_b <- c(Q_LS_FL_di)-c(Q_LS_NL_di)/sd(Q_LS_NL_di)


#################################################################################
### This section is to normalize the effects across assembly scenarios (To run if interested into comparing effects of symmetry independantly across different community
# assembly scenarios) # may be to complex for this paper I think we are leaving this out but I leave this here if needed) 


### Calculating niche-based community  assembly effects
### calculating distances from center


x1 <- matrix(c(Q_EF_MM_w, NODF_EF_MM_w),nrow = length(NODF_EF_MM_w) , ncol =2)
x1.1 <- matrix(c(Q_LS_MM_w, NODF_LS_MM_w),nrow = length(NODF_LS_MM_w) , ncol =2)
x1.2 <- matrix(c(Q_EF_FL_w, NODF_EF_FL_w),nrow = length(NODF_EF_FL_w) , ncol =2)
x1.3 <- matrix(c(Q_LS_FL_w, NODF_LS_FL_w),nrow = length(NODF_LS_FL_w) , ncol =2)


x2 <- matrix(Q_EF_NL_w,NODF_EF_NL_w, nrow = 1, ncol = 2)
x2.2 <- matrix(Q_LS_NL_w,NODF_LS_NL_w, nrow = 1, ncol = 2)


euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))




each <- foreach(i = 1:nrow(x1), .combine = c ) %do% euc.dist(x1[i,],c(0,0))
each2 <- foreach(i = 1:nrow(x1.1), .combine = c ) %do% euc.dist(x1.1[i,],c(0,0))
each3 <- foreach(i = 1:nrow(x1.2), .combine = c ) %do% euc.dist(x1.2[i,],c(0,0))
each4 <- foreach(i = 1:nrow(x1.3), .combine = c ) %do% euc.dist(x1.3[i,],c(0,0))

eachAll <- rbind(
  data.frame("each" = each/max(each),"sym" = 1-col$diff,"SigmaA" = col$Var1, "SigmaB" = col$Var2, "intHyp" = rep("MM", length(col$value)), "EA" = rep("EF", length(col$value))),
  data.frame("each" = each2/max(each2),"sym" = 1-col$diff,"SigmaA" = col$Var1, "SigmaB" = col$Var2, "intHyp" = rep("MM", length(col$value)),"EA" = rep("LS", length(col$value))),
  data.frame("each" = each3/max(each3),"sym" = 1-col$diff,"SigmaA" = col$Var1, "SigmaB" = col$Var2,"intHyp" = rep("FL", length(col$value)),"EA" = rep("EF", length(col$value))),
  data.frame("each" = each4/max(each4),"sym" = 1-col$diff,"SigmaA" = col$Var1, "SigmaB" = col$Var2, "intHyp" = rep("FL", length(col$value)),"EA" = rep("LS", length(col$value)))
)

eachAll2 <- rbind(
  data.frame("each" = each, "sym" = 1-col$diff,"SigmaA" = col$Var1, "SigmaB" = col$Var2, "intHyp" = rep("MM", length(col$value)), "EA" = rep("EF", length(col$value))),
  data.frame("each" = each2, "sym" = 1-col$diff,"SigmaA" = col$Var1, "SigmaB" = col$Var2, "intHyp" = rep("MM", length(col$value)),"EA" = rep("LS", length(col$value))),
  data.frame("each" = each3,"sym" = 1-col$diff,"SigmaA" = col$Var1, "SigmaB" = col$Var2,"intHyp" = rep("FL", length(col$value)),"EA" = rep("EF", length(col$value))),
  data.frame("each" = each4,"sym" = 1-col$diff,"SigmaA" = col$Var1, "SigmaB" = col$Var2, "intHyp" = rep("FL", length(col$value)),"EA" = rep("LS", length(col$value)))
)



### Calculating niche-based interaction assembly effects
boxplot(x1)
# Morphological matching 
y1 <- matrix(c(EF_Q_MM_b, EF_NODF_MM_b),nrow = length(NODF_EF_MM_w) , ncol =2) ## EF
y1.1 <- matrix(c(LS_Q_MM_b, LS_NODF_MM_b),nrow = length(NODF_LS_MM_w) , ncol =2) ## LS
##  Forbiden links 
y1.2 <- matrix(c(EF_Q_FL_b, EF_NODF_FL_b),nrow = length(NODF_EF_FL_w) , ncol =2) ## EF
y1.3 <- matrix(c(LS_Q_FL_b, LS_NODF_FL_b),nrow = length(NODF_LS_FL_w) , ncol =2) ## LS

# Calculate int. assembly variation in Neutral assemblages 
# Morphological matching 
y2 <- matrix(c(Q_SESNull_MM,NODF_SESNull_MM), nrow = 1, ncol = 2)
# Forbidden links
y2.2 <- matrix(c(Q_SESNull_FL,NODF_SESNull_FL), nrow = 1, ncol = 2)

### calculating distances from neutral center
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))


Xeach <- foreach(i = 1:nrow(x1), .combine = c ) %do% euc.dist(y1[i,],y2)
Xeach2 <- foreach(i = 1:nrow(x1.1), .combine = c ) %do% euc.dist(y1.1[i,],y2.2)
Xeach3 <- foreach(i = 1:nrow(x1.2), .combine = c ) %do% euc.dist(y1.2[i,],y2)
Xeach4 <- foreach(i = 1:nrow(x1.3), .combine = c ) %do% euc.dist(y1.3[i,],y2.2)

XeachAll <- rbind(
  data.frame("each" = Xeach/max(Xeach),
             "sym" = 1-col$diff,
             "SigmaA" = col$Var1, 
             "SigmaB" = col$Var2, 
             "intHyp" = rep("MM", length(col$value)), 
             "EA" = rep("EF", length(col$value))),
  data.frame("each" = Xeach2/max(Xeach2),
             "sym" = 1-col$diff,
             "SigmaA" = col$Var1, 
             "SigmaB" = col$Var2, 
             "intHyp" = rep("MM", length(col$value)),
             "EA" = rep("LS", length(col$value))),
  data.frame("each" = Xeach3/max(Xeach3),
             "sym" = 1-col$diff,"SigmaA" = col$Var1,
             "SigmaB" = col$Var2,
             "intHyp" = rep("FL", length(col$value)),
             "EA" = rep("EF", length(col$value))),
  data.frame("each" = Xeach4/max(Xeach4),
             "sym" = 1-col$diff,
             "SigmaA" = col$Var1, 
             "SigmaB" = col$Var2, 
             "intHyp" = rep("FL", length(col$value)),
             "EA" = rep("LS", length(col$value)))
)

XeachAll2 <- rbind(
  data.frame("each" = Xeach,
             "sym" = 1-col$diff,
             "SigmaA" = col$Var1,
             "SigmaB" = col$Var2, 
             "intHyp" = rep("MM", length(col$value)), 
             "EA" = rep("EF", length(col$value))),
  data.frame("each" = Xeach2,
             "sym" = 1-col$diff,
             "SigmaA" = col$Var1,
             "SigmaB" = col$Var2, 
             "intHyp" = rep("MM", length(col$value)),
             "EA" = rep("LS", length(col$value))),
  data.frame("each" = Xeach3,
             "sym" = 1-col$diff,
             "SigmaA" = col$Var1,
             "SigmaB" = col$Var2,
             "intHyp" = rep("FL", length(col$value)),
             "EA" = rep("EF", length(col$value))),
  data.frame("each" = Xeach4,
             "sym" = 1-col$diff,
             "SigmaA" = col$Var1, 
             "SigmaB" = col$Var2, 
             "intHyp" = rep("FL", length(col$value)),
             "EA" = rep("LS", length(col$value)))
)
#################################################################################



#################################################################################
##### Variance Partitioning (Currently with the normalized results across all assembly scenarios, most likely it will be modified to see the results of only Env.Filtering)
#################################################################################

### Relative contributions in the variance of raw metrics
## Fit a good linear models 


lm <- lm(QZscorex~sigmaDif+sigmaA+sigmaB+intHyp+EA, RES) ## full model
lm1 <- lm(QZscorex~intHyp+EA, RES) ## only int hyp and wt 
lm2 <- lm(QZscorex~sigmaDif+intHyp + EA, RES) # int hy and strenght

anova(lm, lm1, lm2)
anova(lm1)



lm <- lm(NODFx~sigmaDif+sigmaA+sigmaB+intHyp+EA, RES) ## full model
lm1 <- lm(NODFx~intHyp+EA, RES) ## only int hyp and wt 
lm2 <- lm(NODFx~sigmaDif+intHyp + EA, RES) # int hy and strenght

anova(lm, lm1, lm2)
anova(lm1)

apsOut=commonalityCoefficients(RES,"QZscorex",
                               list("intHyp", "EA",
                                    "sigmaDif","sigmaA","sigmaB"))
apsOut=commonalityCoefficients(RES,"NODFx",
                               list("intHyp", "EA",
                                    "sigmaDif","sigmaA","sigmaB"))


apsOut=commonalityCoefficients(RES,"QZscorex",
                               list("intHyp", "EA", 
                                    "sigmaDif","sigmaA","sigmaB"))
apsOut2=commonalityCoefficients(RES,"QZscoresd",
                                list("intHyp", "EA", 
                                     "sigmaDif","sigmaA","sigmaB"))

apsOut3=commonalityCoefficients(RES,"NODFx",
                                list("intHyp", "EA", 
                                     "sigmaDif","sigmaA","sigmaB"))
apsOut4=commonalityCoefficients(RES,"NODFsd",
                                list("intHyp", "EA", 
                                     "sigmaDif","sigmaA","sigmaB"))

lm <- lm(QZscorex~sigmaDif+sigmaA+sigmaB+intHyp+EA, RES)
summary(lm)
lm <- lm(QZscoresd~sigmaDif+sigmaA+sigmaB+intHyp+EA, RES)
anova(lm)
lm <- lm(NODFx~sigmaDif+sigmaA+sigmaB+intHyp+EA, RES)
anova(lm)
lm <- lm(NODFsd~sigmaDif+sigmaA+sigmaB+intHyp+EA, RES)
anova(lm)

write.csv(cbind(
  apsOut$CCTotalbyVar,
  apsOut2$CCTotalbyVar,
  apsOut3$CCTotalbyVar,
  apsOut4$CCTotalbyVar
),"commonalityMetrics.csv")

#### End of variance partitioning section 
########################################################################################################################


############################################################
############################################################
## 4) Plot the results.
############################################################
############################################################



# Plot1 and plot 2 
# Plot 1: Space is build with raw Modularity and nestedness values. Has both MM, FL and NL,  
## Neutral assembly is shown as cross 
## Title plot 1: 
## Relative influence of community and interaction assembly on network structure (Note we are only showing the results of environmental filtering)

## Plot 2:  Spaces is build with the effect size, sizes represent the strenght of the filtering 
## at each trophic level. 
## We want to show that the trophic level with the strongest filtering determines the overall variation in 
## network metrics ## to remember to fix the lenght of the axis of the cross based on the neutral variance. 
## Effect of (asymmetry) in the  strenght of community assembly on network structure 


png("PlotSimulationResults.png", height = 3000, width = 3000, pointsize = 100)
plot(NODFx~QZscorex,
     frame = F,
     xlim = c(30, 200),
     ylim = c(-50,130),
     pch = ifelse(intHyp == "MM",
                  15,
                  ifelse(intHyp == "FL", 16,
                         17)),
     col = scales::alpha("red", 0.5),
     ylab = "Nestedness", 
     xlab = "Modularity", 
     cex = c(1-sigmaA) * 2, 
     data = RES1[!RES1$EA == "ND" & RES1$intHyp %in% c("MM", "FL","NL") ,])
points(NODFx~QZscorex,
       pch = ifelse(intHyp == "MM",
                    0,
                    ifelse(intHyp == "FL", 1,
                           2)),
       cex = c(1-sigmaA) * 2, 
       data = RES1[!RES1$EA == "ND" & RES1$intHyp %in% c("MM", "FL","NL") ,])
points(NODFx~QZscorex,
       pch = ifelse(intHyp == "MM",15,
                    ifelse(intHyp == "FL", 16,17)),
       col = scales::alpha("blue", 0.5),
       cex = c(1-sigmaB) * 2, 
       data = RES1[!RES1$EA == "ND",])
points(NODFx~QZscorex,
       pch = ifelse(intHyp == "MM",0,
                    ifelse(intHyp == "FL", 1,2)),
       cex = c(1-sigmaB) * 2, 
       data = RES1[!RES1$EA == "ND",])
abline(h= 0, lty = 2, col = "grey60")

segments("x0" = RES1[RES1$EA == "ND",]$QZscorex-RES1[RES1$EA == "ND",]$QZscoresd,
         "x1" = RES1[RES1$EA == "ND",]$QZscorex+RES1[RES1$EA == "ND",]$QZscoresd,
         "y0" = RES1[RES1$EA == "ND",]$NODFx, 
         "y1" = RES1[RES1$EA == "ND",]$NODFx, lwd = 15, lty = 1)

segments("y0" = RES1[RES1$EA == "ND",]$NODFx-RES1[RES1$EA == "ND",]$NODFsd,
         "y1" = RES1[RES1$EA == "ND",]$NODFx+RES1[RES1$EA == "ND",]$NODFsd,
         "x0" = RES1[RES1$EA == "ND",]$QZscorex, 
         "x1" = RES1[RES1$EA == "ND",]$QZscorex, lwd = 15, lty = 1)
legend("right", c("Neutral interactions", "Morphological match", "Forbidden links"),
       bty ="n", cex = 1,pch = c(2,0,1))
legend("topright", c("Trophic level A", "Trophic level B"),
       bty ="n", cex = 1, fill = scales::alpha(c("red", "blue"), 0.5))
dev.off()



# Alternative plot zooming into MM
# 
# ## Plot 2 
# plot(QZscorex~NODFx,
#      pch = 16, 
#      ylim = c(80, 200),
#      xlab = "Nestedness", 
#      ylab = "Modularity",
#      xlim = c(-50, 30),
#      col = scales::alpha("#0F2080", 0.7),
#      cex = c(1-sigmaA) * 3, 
#      data = RES1[!RES1$EA == "ND" & RES1$intHyp %in% c("MM") ,])
# points(QZscorex~NODFx, pch = 3,
#        lwd = 1, lty =3, cex = 3,data = RES1[RES1$EA == "ND",])
# points(QZscorex~NODFx,
#        pch =16,
#        col = scales::alpha("#F5793A", 0.7),
#        cex = c(1-sigmaB) * 3, 
#        data = RES1[!RES1$EA == "ND" & RES1$intHyp %in% c("MM") ,])
# 
# segments("y0" = RES1[RES1$EA == "ND",]$QZscorex-RES1[RES1$EA == "ND",]$QZscoresd,
#          "y1" = RES1[RES1$EA == "ND",]$QZscorex+RES1[RES1$EA == "ND",]$QZscoresd,
#          "x0" = RES1[RES1$EA == "ND",]$NODFx, 
#          "x1" = RES1[RES1$EA == "ND",]$NODFx)
# 
# segments("x0" = RES1[RES1$EA == "ND",]$NODFx-RES1[RES1$EA == "ND",]$NODFsd,
#          "x1" = RES1[RES1$EA == "ND",]$NODFx+RES1[RES1$EA == "ND",]$NODFsd,
#          "y0" = RES1[RES1$EA == "ND",]$QZscorex, 
#          "y1" = RES1[RES1$EA == "ND",]$QZscorex)
# 
# 
# 
# points(QZscorex~NODFx,
#        cex = c(1-sigmaA) * 3, 
#        data = RES1[!RES1$EA == "ND" & RES1$intHyp %in% c("MM") ,])
# points(QZscorex~NODFx,
#        cex = c(1-sigmaB) * 3, 
#        data = RES1[!RES1$EA == "ND" & RES1$intHyp %in% c("MM") ,])
# legend("bottomleft", c("Filtering strengh consumer", "Filtering strengh resource"),
#        bty ="n", cex = 1,pch = 16,
#        col =c(scales::alpha("#0F2080", 0.9),scales::alpha("#F5793A", 0.7)))
# 
# dev.off()
# ###
# 

#### End of main plot(s) section 
########################################################################################################################


#### Appendix plots 


######
# Appendix plot nestedness and modularity effect sizes distribution 
pdf("Figs/AppendixGrids.pdf", height = 20, width = 10, pointsize = 5)
par(mfrow = c(4,2), oma = c(2,2,2,2))
makeGrid(NODF_EF_NL_w,  main = "NODF_EF_NL", limSim = F )
makeGrid(NODF_EF_MM_w,  main = "NODF_EF_MM", limSim = F)

makeGrid(NODF_LS_NL_w,  main = "NODF_LS_NL", limSim = T)
makeGrid(NODF_LS_MM_w,  main = "NODF_LS_MM", limSim = T)

#makeGrid(Q_EF_NL_w,  main = "Q_EF_NL", limSim = F)
makeGrid(Q_EF_MM_w,  main = "Q_EF_MM", limSim = F)

makeGrid(Q_LS_NL_w,  main = "Q_LS_NL", limSim = T)
makeGrid(Q_LS_MM_w,  main = "Q_LS_MM", limSim = T)
dev.off()

###########


#####
# Appendix plot showing the relative variation of within vs between 

par(mfrow = c(1,2))

plot(c(NODF_EF_MM_w)~EF_NODF_MM_b,
     pch = 16,
     col = "black", 
     ylim = c(-10,10), xlim = c(-50,10))
points(c(NODF_LS_MM_w)~LS_NODF_MM_b, pch = 16,col = "red")
abline(h=1.96, v = 1.96, lty = 2)
abline (h=-1.96 ,v = -1.96, lty = 2)
abline(h=0, v=0)

plot(c(Q_EF_MM_w)~EF_Q_MM_b,
     pch = 16,
     ylim = c(-5,5),
     xlim = c(0,250),
     col = "black")
points(c(Q_LS_MM_w)~LS_Q_MM_b, pch = 16,col = "red")
abline(h=1.96, v = 1.96, lty = 2)
abline (h=-1.96 ,v = -1.96, lty = 2)
abline(h=0, v=0)
##########




###########
# Boxplots showing the raw variation of modularity and nestedness across simulated scenarios



RES <- rbind(RES1,RES2)

RES$intHyp <- factor(RES$intHyp, levels = c("NL","FL","MM"))
RES$EA <- factor(RES$EA, levels = c("ND","SEF","LS"))

par(mfrow = c(1,4), las = 1)
boxplot(RES$QZscorex~RES$intHyp, notch = T, log= "y", xaxt = "n", xlab = "Interaction assembly", ylab = "Modularity" )
axis(1,c(1,2,3), c("Neutral", "Forbidden Links", "Morphological match"))
stripchart(QZscorex~intHyp, data = RES[RES$EA == "SEF",], col = scales::alpha("black", 0.5), vertical = T, add = T, pch  = 16, method = "jitter")
stripchart(QZscorex~intHyp, data = RES[RES$EA == "LS",], col = scales::alpha("red", 0.5), vertical = T, add = T, pch  = 16, method = "jitter")
stripchart(QZscorex~intHyp, data = RES[RES$EA == "ND",], col = scales::alpha("yellow", 1), vertical = T, add = T, pch  = 16, method = "jitter")
legend("topleft", pch = c(16,16,16),col = c("black", "red", "yellow"), 
       legend = c("Environmental filtering", "Limiting similarity", "Neutral"),
       title = "Community assembly")

boxplot(RES$NODFx~RES$intHyp, notch = T, xaxt = "n", xlab = "Interaction assembly", ylab = "Nestedness" )
axis(1,c(1,2,3), c("Neutral", "Forbidden Links", "Morphological match"))
stripchart(NODFx~intHyp, data = RES[RES$EA == "SEF",], col = scales::alpha("black", 0.5), vertical = T, add = T, pch  = 16, method = "jitter")
stripchart(NODFx~intHyp, data = RES[RES$EA == "LS",], col = scales::alpha("red", 0.5), vertical = T, add = T, pch  = 16, method = "jitter")
stripchart(NODFx~intHyp, data = RES[RES$EA == "ND",], col = scales::alpha("yellow", 1), vertical = T, add = T, pch  = 16, method = "jitter")
legend("topleft", pch = c(16,16,16),col = c("black", "red", "yellow"), legend = c("Environmental filtering", "Limiting similarity", "Neutral"),
       title = "Community assembly")

boxplot(RES$NODFx~RES$EA,  xaxt = "n",xlab = "Interaction assembly", ylab = "Nestedness" )
axis(1,c(1,2,3), c("Neutral", "Environmental filtering", "Limiting similarity"))
stripchart(NODFx~EA, data = RES[RES$intHyp == "MM",], col = scales::alpha("red", 0.5), vertical = T, add = T, pch  = 16, method = "jitter")
stripchart(NODFx~EA, data = RES[RES$intHyp == "FL",], col = scales::alpha("blue", 0.5), vertical = T, add = T, pch  = 16, method = "jitter")
stripchart(NODFx~EA, data = RES[RES$intHyp == "NL",], col = scales::alpha("black", 1), vertical = T, add = T, pch  = 16, method = "jitter")
legend("topleft", pch = c(16,16,16),col = c("red", "blue", "black"), legend = c("Morphological matching", "Forbidden links", "Neutral"),
       title = "Interaction assembly")

boxplot(RES$QZscorex~RES$EA,  xaxt = "n",xlab = "Interaction assembly", ylab = "Nestedness" )
axis(1,c(1,2,3), c("Neutral", "Environmental filtering", "Limiting similarity"))
stripchart(QZscorex~EA, data = RES[RES$intHyp == "MM",], col = scales::alpha("red", 0.5), vertical = T, add = T, pch  = 16, method = "jitter")
stripchart(QZscorex~EA, data = RES[RES$intHyp == "FL",], col = scales::alpha("blue", 0.5), vertical = T, add = T, pch  = 16, method = "jitter")
stripchart(QZscorex~EA, data = RES[RES$intHyp == "NL",], col = scales::alpha("black", 1), vertical = T, add = T, pch  = 16, method = "jitter")
legend("topleft", pch = c(16,16,16),col = c("red", "blue", "black"), legend = c("Morphological matching", "Forbidden links", "Neutral"),
       title = "Interaction assembly")






#################
# Grid plots showing the variance in interaction probabilities across simulated scenarios in a trait space 

par(mfrow = c(5,3), 
    mar = c(1,1,3,1),
    oma = c(0,0,2,0))

for(i in c(1:length(traitNet))){
  makeTraitPlot(traitNet[i])
  leg <-  paste(abs(diff(as.numeric(
    str_split(names(traitNet)[i],
              "_",
              simplify = T)[6:7]))), "-",
    abs(mean(as.numeric(
      str_split(names(traitNet)[i],
                "_",
                simplify = T)[6:7]))))
  legend("topleft",
         legend = leg)
  
}



#### End of appendix plot section 
########################################################################################################################

















