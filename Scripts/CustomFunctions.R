### Functions to run the trait-based simulations of bipartite network structure
## Gabriel Muñoz 
## Jan 2020
### Custom functions to run simulations of different scenarios of assembly process and calculate outcomes in network structure


## Load libraries 
library(snowfall)
library(ecolottery)
library(vegan)
library(bipartite)
library(tidyverse)
library(MASS)
library(corrplot)


###########################
###########################

## Generate a species pool for two different trophic levels

## CreateSpPool function to create a regional pool of species given total abundance and richness values
#' @param Jpool Total number of individuals in the pool 
#' @param J Total number of species in the pool 
    
CreateSpPool <- function(Jpool, J, poolShape = c("uniform", "log-series")) {
  
  if (poolShape == "uniform") { 
   
    # With uniform trait values in the species pool
    poola <- data.frame(cbind(1:Jpool[1], rep(1:J[1], Jpool[1]/J[1])))
    poola <- poola[order(poola$X2),]
    poola$X3 <- c(sapply(runif(J[1],0,1), function(x) rnorm(10, x, 0.001)))
    colnames(poola) = c("ind"  ,"sp",   "trait")
    
    
    poolb <- data.frame(cbind(1:Jpool[2], rep(1:J[2], Jpool[2]/J[2])))
    poolb <- poolb[order(poolb$X2),]
    poolb$X3 <- c(sapply(runif(J[2],0,1), function(x) rnorm(10, x, 0.001)))
    colnames(poolb) = c("ind"  ,"sp",   "trait")
    
    
  }
  
  if(poolShape == "log-series"){
    # With log-series trait values in the speciespool
    
    poola <- coalesc(Jpool[1], m = 0.5, theta = 50)$pool
    poolb <- coalesc(Jpool[2], m = 0.5, theta = 50)$pool
    
  }

  
  
 
  return(list( "poolA" = poola, "poolB" = poolb))
  
}



###########################
###########################


###### Helper function to calculate null models for modularity and nestedness (From Munoz et al. 2019. JBI )
#' NullModSen: Helper Function to calculate mean and sd of null models of modularity and nestedness for a given matrix/network
#'
#' @param matrix A network organized as a matrix in with species A as columns  Species B as rows.
#' @param rep  integer specifying the number of replicates for the null model matrices, i.e. shuffles of the web, passes to bipartite::nullmodel
#' @return a list with the mean and sd of calculated modularity and nestedness for the null model.
#' @examples NullModSen(matrix,rep)

#' createMetaMatrix: Helper function to create a metaMatrix simulating interactions as probability of relative abundances of species of each trophic level 
#' @param Plants Regional pool of plant trophic level resulting from the coalescent model from the Big species pool
#' @param Pollinators Regional pool of Pollinators |or any other| trophic level resulting from the coalescent model from the Big species pool
#' @param intHyp Hypothesis to assembly interactions 
#' @return a MetaMatrix of interactions as an edgelist with variables describing pairwise interaction probabilities and traits 

###########
NullModSen =  function(matrix, rep){
  
  # Generate null models "r1 null model": non-sequential algorithm for binary matrices that preserves the site (row) frequencies, but uses column marginal frequencies as probabilities of selecting species.
  
  nnm = simulate(vegan::nullmodel(matrix, "r1"), rep)
  
  # Calculate nestedness 
  null.nest = sapply(1:rep, function(x) bipartite::nested(nnm[,,x], method = "NODF2"))
  
  # Calculate modularity 
  null.mod = sapply(1:rep, function(x) bipartite::LPA_wb_plus(nnm[,,x])$modularity)
  
  null = list(c("NulMod" = c("mean" = mean(null.mod), "sd" = sd(null.mod)),
                "NulNest" = c("mean" = mean(null.nest), "sd" = sd(null.nest))))
  return(null)
  
}
#############


## Define gaussian function to set the strenght of the environmental filtering 
##########
filt_gaussian <- function(t,x) exp(-(x-t)^2/(2*sigma^2))
#########



###########################
###########################


#'  Wrapper function to simulate community assembly mechanisms, simulate interactions and calculate network metrics, many parameters pass on different functions of to the `coalesc` function from the `eccolottery` package
#' @param sigma Degree of spacing from the species pool. Modulating it would mean the "strenght" of the environmental filtering process
#' @param repNull Number of replicates for the null model network metrics 
#' @param quantile Define threshold defining the spliting to transform the weighted networks to bipartite networks
#' @param TraitA Trait mean of A trophic level
#' @param TraitB Trait mean of B trophic level
#' @param mA Migration rate from A trophic level 
#' @param mB Migration rate from B trophic level 
#' @param JA Number of species in the regional pool for A
#' @param JB Number of species in  regional pool for B
#' @param EF.A  assembly process for A
#' @param EF.B  assembly process for B 
#' @param intHyp Interaction assembly hypothesis between A and B 

simulateEF <- function(TraitA, TraitB, mA, mB, JA, JB, poolA, poolB, EF.A = c("ND","SEF", "DEF", "LS"), 
                       EF.B = c("ND","SEF", "DEF", "LS"), sigmaA, sigmaB, repNull = 10, quantile = 10, 
                       chunkA = chunkA, chunkB = chunkB, 
                       intHyp = c("MM", "FL", "NL")) {
  
  
  comA <- SimulateSelCom(TraitA, 
                         mA, 
                         JA,
                         poolA,
                         EF = EF.A, 
                         sigma = sigmaA,
                         chunk = chunkA)
  comB <- SimulateSelCom(TraitB, 
                         mB, 
                         JB, 
                         poolB,
                         EF = EF.B, 
                         sigma = sigmaB, 
                         chunk = chunkB)
  
  # Simulate interactions 
  metaMatrix <- createMetaMatrix(comA$com , 
                                 comB$com, 
                                 intHyp = intHyp)
  
  #met <- tryCatch({calculateMetrix(metaMatrix,quantile = quantile, repNull, intHyp = intHyp)}, 
  #              error = function(e) { print("NA")})
  met <- calculateMetrix(metaMatrix,
                         quantile = quantile,
                         repNull = repNull)
  
  
  # return object 
  
  return(met)
  
  
}

###########################
###########################



###########################
###########################


### Slight variation in SimulateEF to create networks but not to calculate paramenters 

simulateEF2 <- function(TraitA, TraitB, mA, mB, JA, JB, poolA, poolB, EF.A = c("ND","SEF", "DEF", "LS"), 
                       EF.B = c("ND","SEF", "DEF", "LS"), sigmaA, sigmaB, repNull = 10, quantile = 10, 
                       chunkA = chunkA, chunkB = chunkB, 
                       intHyp = c("MM", "FL", "NL")) {
  
  
  comA <- SimulateSelCom(TraitA, 
                         mA,
                         JA,
                         poolA,
                         EF = EF.A,
                         sigma = sigmaA, 
                         chunk = chunkA)
  comB <- SimulateSelCom(TraitB,
                         mB, 
                         JB, 
                         poolB, 
                         EF = EF.B, 
                         sigma = sigmaB, 
                         chunk = chunkB)
  
  # Simulate interactions 
  metaMatrix <- createMetaMatrix(comA$com ,
                                 comB$com, 
                                 intHyp = intHyp)
  
  # Create the weighted meta-network 
  intData <- data.frame(reshape2::acast(
    data = metaMatrix, 
    Var1~Var2, 
    sum,
    value.var = "intProb"))
  
  print("Subsetting metaMatrix on different int.Prob threshold") # Harcode to only select as realized interactions those whitin the 5th highest percentile of pairwise interaction probability
  metaMatrix <- metaMatrix %>%
    mutate(quantile = ntile(metaMatrix$intProb, 5))
  # subset network based on interaction threshold 
  
  subSet <-  metaMatrix[metaMatrix$quantile >= 5,]
  
  # Make binary matrix
  
  intData2 <- table(subSet$Var1, subSet$Var2)
  attributes(intData2)$class <- "matrix"
  
  # return object 
  
  return(list("net" = intData2,
              "com" = na.omit(subSet)))
  
  
}


###########################
###########################


###########################
###########################


### Helper function to call ecolottery package functions to simulate assembly processes in ecological communities
#' All parameters come from wrapper function "simulateEF" --- see above --- 
#########

SimulateSelCom  <- function(Trait, m, J, pool, 
                            EF = c("ND","SEF", "DEF", "LS"), 
                            sigma = sigma, 
                            chunk = chunk){
  
  
  # Simulate environmental filtering - Stabilizing filtering
  
  if (EF == "SEF") { 
    
    sigma <- sigma
    filt_gaussian <- function(t,x) exp(-(x-t)^2/(2*sigma^2))
    
    comA <- coalesc(J, m, filt = function(x) filt_gaussian(Trait, x), pool = pool)
    
    
  }
  
  # Simulate environmental filtering - Directional filtering 
  
  if (EF == "DEF") {
    
    
    comA <- coalesc(J, m, filt = function(x) 1 - min(x,1), pool = pool)
    
    
  }
  
  # Simulate neutral dynamics 
  
  if (EF == "ND") { 
    
    comA <- coalesc(J, m, pool = pool) 
    
    
  }
  
  # Simulate limiting similarity 
  
  if (EF == "LS") {
    
    overlap <- seq(0.1, 0.001, length.out = 5)
    
    num <- density(limSim(J, 
                          partition = chunk, 
                          overlap = overlap[chunk-1])) 
    # because chunks should always start at 2 and max 6 (harcoded, but to improve if really necessary)
           
    comA <- coalesc(J, m, filt = function(t) num$y[which.min(abs(num$x - t))] , pool = pool)
    
  }
  
  return(comA)
  
}

###########################
###########################



###########################
###########################

### Function to simulate effects of limiting similarity 

limSim <- function(J,partition,overlap){
  lim <-c(sapply(seq(0,1,length.out = partition),
           function(x) rnorm(J, x, overlap)))
  lim <- ifelse(lim < 0,0,
                ifelse(lim > 1, 1, lim))
  return(lim)}


###########################
###########################

###########################
###########################
#' calculateMetrix: function calculate a meta-matrix of interactions after the community assembly
#' @param Plants: Resource communtity
#' @param Pollinators: Consumer communtity
#' @param intHyp: Interaction hypothesis (MM = Morphological Match, FL = Forbidden Links, NL = Neutral interactions (i.e. abundance based))
#' @return a matrix of pairwise interactions in long format

createMetaMatrix <- function(Plants, Pollinators, intHyp = c("MM", "FL", "NL")){ 
  
  ## Coercing column names
  names(Plants) <- c("id", "sp", "trait")
  names(Pollinators) <- c("id", "sp", "trait")
  
  ### Simulating interactions 
  print("Simulating interactions")
  metaMatrix <- expand.grid(unique(Plants$id), unique(Pollinators$id))
  metaMatrix$relAbPl <- prop.table(table(Plants$id))[match(metaMatrix$Var1, names(prop.table(table(Plants$id))))]
  metaMatrix$trait1 <- Plants$trait[match(metaMatrix$Var1, Plants$id)]
  metaMatrix$sp1 <- Plants$sp[match(metaMatrix$Var1, Plants$id)]
  
  metaMatrix$relAbPol <- prop.table(table(Pollinators$id))[match(metaMatrix$Var2, names(prop.table(table(Pollinators$id))))]
  metaMatrix$trait2 <- Pollinators$trait[match(metaMatrix$Var2, Pollinators$id)]
  metaMatrix$sp2 <- Pollinators$sp[match(metaMatrix$Var2, Pollinators$id)]
  
  metaMatrix <- na.omit(metaMatrix) # Erase non interacting species
  
  metaMatrix$intProb <- metaMatrix$relAbPl*metaMatrix$relAbPol ### As a probability of abundance (i.e. int. probability equals the product of both relative abundances)
  
  #metaMatrix$intProb <- metaMatrix$relAbPl^2/(metaMatrix$relAbPl^2+metaMatrix$relAbPol^2) # disregard this part and not uncomment. Was a tryout 

  # Change interaction probabilities based on interaction assembly hypothesis 
  
    if (intHyp == "FL"){
      ## Forbidden links 1-trait <= 0 --> int prob = 0 
      metaMatrix$intProb <- ifelse(metaMatrix$trait1-metaMatrix$trait2 >= 0, metaMatrix$intProb * 1, metaMatrix$intProb * 0)
    } 
    
    if (intHyp == "MM") {
      ## Morphological matching, interaction depends on the frequency of differences between traits 
      metaMatrix$intProb <- 1- abs(metaMatrix$trait1-metaMatrix$trait2)/max(abs(metaMatrix$trait1-metaMatrix$trait2))
      
    }
  
  if (intHyp == "NL") {
    # Neutral assembly of interactions, probability of interactions depend solely on abundance (rel.abundance of A * rel.abundance of B)
    metaMatrix$intProb <- metaMatrix$intProb * 1
    }
  
  
  return(metaMatrix)
  
}



###########################
###########################

###########################
###########################

#' calculateMetrix: function to prune interaction potential, create a binary bipartite matrix and calculate network metrics for each quatiles of realized interactions 
#' @param metaMatrix An edge list of interactions, passed from the results of createMetaMatrix function 
#' @param quatile Number of quantiles to split the potential interactions 
#' @param repNull Number of matrix replicates to create the nullmodels for nestedness and modularity
#' @return a list object with "metrics" containing the network metrics calculated and  "wnet" the realized network of each replicate 

##########
calculateMetrix <- function(metaMatrix, quantile, repNull = 100){
 
  # Create the weighted meta-network 
  intData <- data.frame(reshape2::acast(
    data = metaMatrix, 
    sp1~sp2, 
    sum,
    value.var = "intProb"))
  
    print("Subsetting metaMatrix on different int.Prob threshold")
    metaMatrix <- metaMatrix %>%
      mutate(quantile = ntile(metaMatrix$intProb, quantile))
    
    res <- c()
  
    for(i in 1:quantile){ 
      
      # Create binary matrices
      print(i)
      
      # subset network based on interaction threshold 
      
      subSet <-  metaMatrix[metaMatrix$quantile >= i,]
      # Make binary matrix
    
      intData2 <- table(subSet$sp1, subSet$sp2)
      intData2[intData2 >= 1] <- 1
      attributes(intData2)$class <- "matrix"
      
      print(dim(intData2))
      
      # null models
      null1 = NullModSen(intData2, repNull)
      
      # Calculate  z-score Modularity
      zMod.net1 = (bipartite::LPA_wb_plus(intData2)$modularity - null1[[1]]["NulMod.mean"]) / null1[[1]]["NulMod.sd"]
      
      # Calculate z-score Nestedness 
      
      zNes.net1 =  (bipartite::nested(intData2, method = "NODF2") - null1[[1]]["NulNest.mean"]) / null1[[1]]["NulNest.sd"]
      
      res[[i]] <- list("metrics" = data.frame("NODF Z-score" = zNes.net1,
                             "Q Z-score" = zMod.net1 ,
                             "con" = networklevel(intData2, weighted = F, index = c("connectance")),
                             "cor" = cor(subSet$trait1, subSet$trait2, method = "spearman")),
                       "bnet" = subSet, 
                       "specialization" = bipartite::nodespec(intData2, inf.replace = Inf) )
    }
    


  res <- list("metrics" = res, "wnet" = intData)
  return(res)
}
###########

###########################
###########################

###########################
###########################


### Wrapper function to run the simulations 


RunSimul <- function(Pool, Scenario, replicates = 10, repNull, quantile, runParal = T){
  
  
  if (runParal == T){
    
    
    SIMS <- sfLapply(1:length(Scenario$EA), function(x) replicate(replicates , # number of replicates per scenario
                                                                        simulateEF(TraitA = Scenario[x,]$TraitA,
                                                                                   TraitB = Scenario[x,]$TraitB,
                                                                                   mA = Scenario[x,]$mA,
                                                                                   mB = Scenario[x,]$mB,
                                                                                   chunkA = Scenario[x,]$chunkA,
                                                                                   chunkB = Scenario[x,]$chunkB,
                                                                                   JA = Scenario[x,]$JA,
                                                                                   JB = Scenario[x,]$JB,
                                                                                   EF.A = Scenario[x,]$EA,
                                                                                   EF.B = Scenario[x,]$EB,
                                                                                   poolA = Pool$poolA,
                                                                                   poolB = Pool$poolB,
                                                                                   sigmaA = Scenario[x,]$sigmaA,
                                                                                   sigmaB = Scenario[x,]$sigmaB,
                                                                                   repNull = repNull,
                                                                                   quantile = quantile,
                                                                                   intHyp =  Scenario[x,]$intHyp)
    )
    )
    
    
  }
  
  if (runParal == F){
    
    SIMS <- lapply(1:length(Scenario$EA), function(x) replicate(replicates , # number of replicates per scenario
                                                                      simulateEF(TraitA = Scenario[x,]$TraitA,
                                                                                 TraitB = Scenario[x,]$TraitB,
                                                                                 mA = Scenario[x,]$mA,
                                                                                 mB = Scenario[x,]$mB,
                                                                                 chunkA = Scenario[x,]$chunkA,
                                                                                 chunkB = Scenario[x,]$chunkB,
                                                                                 JA = Scenario[x,]$JA,
                                                                                 JB = Scenario[x,]$JB,
                                                                                 EF.A = Scenario[x,]$EA,
                                                                                 EF.B = Scenario[x,]$EB,
                                                                                 poolA = Pool$poolA,
                                                                                 poolB = Pool$poolB,
                                                                                 sigmaA = Scenario[x,]$sigmaA,
                                                                                 sigmaB = Scenario[x,]$sigmaB,
                                                                                 repNull = repNull,
                                                                                 quantile = quantile,
                                                                                 intHyp =  Scenario[x,]$intHyp)
    )
    )
    
    
    
  }
  
  
  
  ## Give appropiate names to the resulting object 
  
  
  names(SIMS) <- paste0(Scenario$EA,"_",
                        Scenario$EB, "_",
                        Scenario$intHyp, "_",
                        Scenario$TraitA, "_", 
                        Scenario$TraitB,"_",
                        Scenario$sigmaA,"_", 
                        Scenario$sigmaB,"_",
                        Scenario$chunkA,"_",
                        Scenario$chunkB)
  
  
  
  
  return(SIMS)
  
  
}

###########################
###########################
###########################
###########################


### Wrapper function to generate networks with varying trait optima to use as test of the null model approach  


GenBiNet <- function(Pool, Scenario, replicates = 10, repNull, quantile, runParal = T){
  
  if (runParal == T){
    
    
    SIMS <- sfLapply(1:length(Scenario$EA), function(x) replicate(replicates , # number of replicates per scenario
                                                                  simulateEF2(TraitA = Scenario[x,]$TraitA,
                                                                             TraitB = Scenario[x,]$TraitB,
                                                                             mA = Scenario[x,]$mA,
                                                                             mB = Scenario[x,]$mB,
                                                                             chunkA = Scenario[x,]$chunkA,
                                                                             chunkB = Scenario[x,]$chunkB,
                                                                             JA = Scenario[x,]$JA,
                                                                             JB = Scenario[x,]$JB,
                                                                             EF.A = Scenario[x,]$EA,
                                                                             EF.B = Scenario[x,]$EB,
                                                                             poolA = Pool$poolA,
                                                                             poolB = Pool$poolB,
                                                                             sigmaA = Scenario[x,]$sigmaA,
                                                                             sigmaB = Scenario[x,]$sigmaB,
                                                                             repNull = repNull,
                                                                             quantile = quantile,
                                                                             intHyp =  Scenario[x,]$intHyp)
    )
    )
    
    
  }
  
  if (runParal == F){
    
    SIMS <- lapply(1:length(Scenario$EA), function(x) replicate(replicates , # number of replicates per scenario
                                                                simulateEF2(TraitA = Scenario[x,]$TraitA,
                                                                           TraitB = Scenario[x,]$TraitB,
                                                                           mA = Scenario[x,]$mA,
                                                                           mB = Scenario[x,]$mB,
                                                                           chunkA = Scenario[x,]$chunkA,
                                                                           chunkB = Scenario[x,]$chunkB,
                                                                           JA = Scenario[x,]$JA,
                                                                           JB = Scenario[x,]$JB,
                                                                           EF.A = Scenario[x,]$EA,
                                                                           EF.B = Scenario[x,]$EB,
                                                                           poolA = Pool$poolA,
                                                                           poolB = Pool$poolB,
                                                                           sigmaA = Scenario[x,]$sigmaA,
                                                                           sigmaB = Scenario[x,]$sigmaB,
                                                                           repNull = repNull,
                                                                           quantile = quantile,
                                                                           intHyp =  Scenario[x,]$intHyp)
    )
    )
    
    
    
  }
  
  
  
  ## Give appropiate names to the resulting object 
  
  names(SIMS) <- paste0(Scenario$EA,"_",
                        Scenario$EB, "_",
                        Scenario$intHyp, "_",
                        Scenario$TraitA, "_", 
                        Scenario$TraitB,"_",
                        Scenario$sigmaA,"_", 
                        Scenario$sigmaB,"_",
                        Scenario$chunkA,"_",
                        Scenario$chunkB)
  
  
  
  
  
  
  return(SIMS)
  
  
}


###########################
###########################

###########################
###########################



########
## End of functions for the simulation modelling approach
########


######
# Start of functions to read and interpret the output
#######


###########################
###########################

## sumarizeResults: function to summarize  the simulation object into an interpretable dataframe 


sumarizeResults <- function(Scenario, simul, quantile, replicate){
  
  ## Generate a new dataframe to plot results (based on medians and sd)
  RES1 <- Scenario
  RES1$treat <- paste0(Scenario$EA,"_", Scenario$EB)
  # Nestedness
  RES1$NODFx <- sapply(1:length(simul),
                       function(y) median(sapply(1:replicate, 
                                                 function(x)simul[[y]][,x]$metrics[[quantile]]$metrics$NODF.Z.score), 
                                          na.rm = T))
  RES1$NODFsd <-  sapply(1:length(simul), 
                         function(y) sd(sapply(1:replicate, 
                                               function(x)simul[[y]][,x]$metrics[[quantile]]$metrics$NODF.Z.score),
                                        na.rm = T))
  # Modularity
  RES1$QZscorex <-sapply(1:length(simul),
                         function(y) median(sapply(1:replicate,
                                                   function(x)simul[[y]][,x]$metrics[[quantile]]$metrics$Q.Z.score),
                                            na.rm = T))
  RES1$QZscoresd  <- sapply(1:length(simul), 
                            function(y) sd(sapply(1:replicate,
                                                  function(x)simul[[y]][,x]$metrics[[quantile]]$metrics$Q.Z.score), 
                                           na.rm = T))
  # Connectance
  RES1$conX <- sapply(1:length(simul),
                      function(y) median(sapply(1:replicate,
                                                function(x)simul[[y]][,x]$metrics[[quantile]]$metrics$con), 
                                         na.rm = T))
  RES1$conSD  <-sapply(1:length(simul), 
                       function(y) sd(sapply(1:replicate, 
                                             function(x)simul[[y]][,x]$metrics[[quantile]]$metrics$con), 
                                      na.rm = T))
  # Spearman rank correlation
  RES1$corX <- sapply(1:length(simul), 
                      function(y) median(sapply(1:replicate,
                                                function(x)simul[[y]][,x]$metrics[[quantile]]$metrics$cor), 
                                         na.rm = T))
  RES1$corSD  <- sapply(1:length(simul), 
                        function(y) sd(sapply(1:replicate, 
                                              function(x)simul[[y]][,x]$metrics[[quantile]]$metrics$cor),
                                       na.rm = T))
  
  #  trait niche breath of interacting partners 
  
  RES1$avINBA <- apply(sapply(1:length(simul), 
                              function(y) sapply(1:replicate, 
                                                 function(x)
                                                   median(aggregate(simul[[y]][,x]$metrics[[quantile]]$bnet$trait2, 
                                                                    list(simul[[y]][,x]$metrics[[quantile]]$bnet$Var1) ,
                                                                    FUN = max)$x -
                                                            aggregate(simul[[y]][,x]$metrics[[quantile]]$bnet$trait2, 
                                                                      list(simul[[y]][,x]$metrics[[quantile]]$bnet$Var1) , 
                                                                      FUN = min)$x))), MARGIN = 2, mean)
  # Average trait niche breath of interacting partners 
  
  RES1$avINBB <- apply(sapply(1:length(simul), 
                              function(y) sapply(1:replicate, 
                                                 function(x)
                                                   median(aggregate(simul[[y]][,x]$metrics[[quantile]]$bnet$trait1, 
                                                                    list(simul[[y]][,x]$metrics[[quantile]]$bnet$Var2) ,
                                                                    FUN = max)$x -
                                                            aggregate(simul[[y]][,x]$metrics[[quantile]]$bnet$trait1, 
                                                                      list(simul[[y]][,x]$metrics[[quantile]]$bnet$Var2) , 
                                                                      FUN = min)$x))), MARGIN = 2, mean)
  #  Ratio of interaction niche trait breath  
  
  RES1$INBr <- log(RES1$avINBA/RES1$avINBB)  # Non necesary anymore 
  
  # Total trait niche of interacting partners 
  RES1$tINBA <- apply(sapply(1:length(simul), 
                              function(y) sapply(1:replicate, 
                                                 function(x)
                                                   max(aggregate(simul[[y]][,x]$metrics[[quantile]]$bnet$trait2, 
                                                                    list(simul[[y]][,x]$metrics[[quantile]]$bnet$Var1) ,
                                                                    FUN = max)$x) -
                                                            min(aggregate(simul[[y]][,x]$metrics[[quantile]]$bnet$trait2, 
                                                                      list(simul[[y]][,x]$metrics[[quantile]]$bnet$Var1) , 
                                                                      FUN = min)$x))), MARGIN = 2, mean)
  
  RES1$tINBB <- apply(sapply(1:length(simul), 
                             function(y) sapply(1:replicate, 
                                                function(x)
                                                  max(aggregate(simul[[y]][,x]$metrics[[quantile]]$bnet$trait1, 
                                                                list(simul[[y]][,x]$metrics[[quantile]]$bnet$Var2) ,
                                                                FUN = max)$x) -
                                                  min(aggregate(simul[[y]][,x]$metrics[[quantile]]$bnet$trait1, 
                                                                list(simul[[y]][,x]$metrics[[quantile]]$bnet$Var2) , 
                                                                FUN = min)$x))), MARGIN = 2, mean)
  
  
  return(RES1)
  

  
}



###########################
###########################



########
## End of functions to interpret the read the output of the simulation approach
########



########
## Start of functions and helper for plotting the results
########

###########################
###########################

### Make a palette vector with a vector column data 


f <- function(x,n=10, pal, rev = F){
  if(rev == F){ 
    rev(RColorBrewer::brewer.pal(n, pal))[cut(x,n)]
  }else{
    (RColorBrewer::brewer.pal(n, pal))[cut(x,n)]
  }
}

###########################
###########################


###########################
###########################


### Plot individual grids 

makeMesh <- function(RES1, intHyp, value.var, res = T, scale = T){ 
  
  SEF.fl <- RES1[RES1$treat != "ND_ND" & RES1$intHyp == intHyp,]
  
  if(res == T){
    
    SEF.fl$tINBres <- resLoww(SEF.fl$tINBA, SEF.fl$tINBB)
    
  }
  
  
  tData <- reshape2::acast(
    data = SEF.fl, 
    sigmaA~sigmaB, 
    sum,
    value.var = value.var)
  
  
  
  if (scale == T) {
    re = (tData - mean(tData))/sd(tData)
    re = re/max(re)
    return(re)
  }
  
  else {
    
    return(tData)
    
  }
  
  
  
  
}



###########################
###########################

###########################
###########################

## Rotate a matrix function
# from https://stackoverflow.com/questions/16496210/rotate-a-matrix-in-r
rotate <- function(x) t(apply(x, 2, rev))

makeGrid <- function(scal, main , limSim = F){
  pal2 = c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac')
  
  pval =  abs(rotate(scal))
  # Get positions & plot formatted p-values
  pos <- expand.grid(1:ncol(pval), ncol(pval):1)
  pos$Var1 <- 6-pos$Var1
  pos$Var2 <- 6-pos$Var2
  
  corrplot::corrplot(rotate(scal),
                     outline = F,
                     mar = c(2,2,2,2), 
                     tl.pos = "n",
                     method = "color",
                     is.corr = F,
                     p.mat = pval, 
                     insig = "pch",
                     pch = "*",
                     pch.cex = 5,
                     pch.col = "white",
                     sig.level = 1.96,
                     col = rev(pal2),
                     cl.offset = 12,
                     tl.cex = 3,
                     cl.cex = 2,
                     cl.lim=c(-8,10)
                     # cl.lim=c(min(scal),max(scal))
                     )
  title(main, cex.main = 3)
  axis(1, c(0.5,5.5), tick = T, labels = c("low", "high"), line = -2, cex.axis = 1.5)
  axis(2, line = 2, c(0.5,5.5), tick = T, labels = c("low", "high"), cex.axis = 1.5, las = 2)
  
  if (limSim == T){ 
    axis(1, c(3), tick = F, labels = c("Strenght of limiting similarity - Consumer"), line = 2, cex.axis = 3)
    axis(2, c(3), tick = F, labels = c("Strenght of limiting similarity - Resource"), line = 2, cex.axis = 3)
    }
  
  else{
    axis(1, c(3), tick = F, labels = c("Strenght of filtering - Consumer"), line = 2, cex.axis = 3)
    axis(2, c(3), tick = F, labels = c("Strenght of filtering - Resource"), line = 2, cex.axis = 3)
  }
  
  
}



###########################
###########################

###########################
###########################

### make transposed matrix under lim sim simulations

makeMesh2 <- function(RES1, intHyp, value.var, res = T, scale = T){ 

SEF.fl <- RES1[RES1$treat != "ND_ND" & RES1$intHyp == intHyp,]

if(res == T){
  
  SEF.fl$tINBres <- resLoww(SEF.fl$tINBA, SEF.fl$tINBB)
  
}


tData <- reshape2::acast(
  data = SEF.fl, 
  chunkA~chunkB, 
  sum,
  value.var = value.var)



if (scale == T) {
  re = (tData - mean(tData))/sd(tData)
  re = re/max(re)
  return(re)
}

else {
  
  return(tData)
  
}




}


###########################
###########################



###########################
###########################
# Plot a covex hull from a list of given vectors 



Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor, lwd = 3)
}  


###########################
###########################


###########################
###########################


### Function to export bipartite network data and trait communities 

getBnet <- function(mmMat){
  subSet <-  mmMat[mmMat$quantile >= 5,]
  # Make binary matrix
  intData2 <- table(subSet$Var1, subSet$Var2)
  attributes(intData2)$class <- "matrix"
  ret = list (
    "bnet" = intData2,
    "tlist" = na.omit(subSet))
  
  return(ret)
}





###########################
###########################

makeTraitPlot <- function(traitNet){ 
  nam <- names(traitNet)
  nam <- data.frame(stringr::str_split(nam, "_", simplify = T), stringsAsFactors = F)
  names(nam) <- c("filt1", "filt2", "intFil", "trait1", "trait2", "sigma1", "sigma2", "X", "X")
  
  
  for(i in 1:length(traitNet)){ 
    
    sd <- traitNet[[i]][,sample(1:10, 1)]$com
    
    names(traitNet[[1]][,1]$com)
    plot(0, col = "white",
         xaxt = "n",
         yaxt = "n",
         xlab = paste0("trait1, sigma = ", nam$sigma1[i]),
         ylab = paste0("trait2, sigma = ", nam$sigma2[i]),
         main = paste(nam$filt1[i], nam$filt2[i], nam$intFil[i]),
         xlim = c(0,1),
         ylim = c(0,1))
    rect(0,1,1,0, density = 90, col = "gray80")
    points(c(sd$trait1),c(sd$trait2),
           pch = 16,
           cex = (sd$relAbPl/max(sd$relAbPl)),
           col = scales::alpha(f(c(sd$intProb/max(sd$intProb)),10,"Spectral"),0.1))
    
    abline(v= c(nam$trait1[i]),h = c(nam$trait2[i]))
    abline(0,1, lty = 2, col = "gray90")
    
    
  }
}




###########################
###########################

makeIntPlot <- function(traitNet){ 
  nam <- names(traitNet)
  nam <- data.frame(stringr::str_split(nam, "_", simplify = T), stringsAsFactors = F)
  names(nam) <- c("filt1", "filt2", "intFil", "trait1", "trait2", "sigma1", "sigma2", "X", "X")
  
  
  for(i in 1:length(traitNet)){ 
    
    sd <- traitNet[[i]][,sample(1:10, 1)]$com
    
    names(traitNet[[1]][,1]$com)
    
    plot(c(sd$intProb)~c(sd$trait1/sd$trait2),
         cex= log1p(sd$relAbPol/sd$relAbPl))
    
    
    
    
    
    
  }
}



###########################
###########################


## Get the residuals from of a loweess smoothed fit 

resLoww <- function(x,y){
  # first set up a lowess fit:
  lfit <- lowess(x,y)
  
  # create a functional version of the lowess fit
  lfun <- approxfun(lfit)
  fitted <- lfun(x)
  resid <- y-fitted
  return(resid)
  
}

###########################
###########################
## Function to calculate Modularity Z-scores
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

######
# function to calculate nestedness z-scores 

function(matrix, nullSim){
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







