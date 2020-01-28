### Proposed implementation for the case-study

## Gabriel Muñoz 
# January 2020

source("scripts/CustomFunctions.R") # path might change 

############################################################
############################################################
## 1) define paramter objects and create combinations of assembly scenarios to simulate
############################################################
############################################################


## dimensions of metanetwork
length(unique(intRQL$plantCode))
length(unique(intRQL$animalCode))


## simulate a species pool with uniform distribution of traits 

simPool = CreateSpPool(c(1020, 1270), c(102,127), "uniform")

hist(simPool$poolA$tra1, col = "red")
points(density(simPool$poolA$trait), col = "blue", type = "l")

## correct the proper pool distribution of traits 

hist(RLQ$mQ$NorS1)
hist(RLQ$mR$NorS1)

length(RLQ$mR$NorS1)
## 1) How do we assign relative abundances? do we assume more generalist are also more abundant? 
## all species have the same relative abundances 

## considering no considerable interspecific variation 
simPool$poolA$trait <- sapply(simPool$poolA$sp, function(x) RLQ$mR$NorS1[x])
simPool$poolB$trait <- sapply(simPool$poolB$sp, function(x) RLQ$mQ$NorS1[x])

## Now trait distribution in the simulated pools match the distribution of traits in the metanetwork
hist(simPool$poolA$trait)
hist(simPool$poolB$trait)

## Now lets simulate similar sized networks for each elevation from the simulated pool. 

# Set parameters from know observations 
# trait optimas at each elevation 

lnet$anT <- RLQ$mQ$NorS1[match(lnet$animalCode,rownames(RLQ$mQ))]
lnet$plT <- RLQ$mR$NorS1[match(lnet$plantCode,rownames(RLQ$mR))]

mean(lnet$anT) # trait optima a
mean(lnet$plT) # trait optima b 

sd(lnet$anT) # strenght filtering a
sd(lnet$plT) # strenght filtering b 

## Example for network at 1000 elevation 

PossCS1000  <- c("EA" = "SEF", 
                 "EB" = "SEF", 
                 "intHyp" = "MM", 
                 "TraitA" = 0.141533,
                 "TraitB" = 0.2264987,
                 "sigmaA" =  0.8859347, # strenght of filtering at trophic level A
                 "sigmaB" = 0.6875488 ,# strenght of filtering at trophic level B
                 "JA" = 500, # individuals at local community at trophic level A
                 "JB" = 500,  # individuals at local community at trophic level B
                 "mA" = 0.5,  # migration parameter (for neutral scenarios only, not used for niche-based assembly but needed in the object)
                 "mB" = 0.5 )  # migration parameter (for neutral scenarios only))


CaseStudySim <- RunSimul(Pool = simPool,# species pool
                         data.frame(t(PossCS1000)), # assembly scenario combinations (# modify this line to change the assembly scenario accordingly )
                  replicates = 10, # number of replicates of each assembly scenario
                  repNull = 10, # number of replicates of reshuffles to calculate Z-scores for network metrics 
                  quantile = 5, # quantile to consider "realized interactions" 
                  runParal = F) # Will it run in parallel?  (Strongly recomended)











