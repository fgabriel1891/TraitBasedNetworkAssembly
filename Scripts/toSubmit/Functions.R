###################################
## Code with custom functions used to perform the simulations and analyse of the results shown in: Trait-based inference of ecological network assembly
## Sub. to: Ec. letters 
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

## Load libraries 
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

#######
# Functions to create null models 
######

sortMATRIX <- function(MATRIX,binary,sortvar) {
  
  
  rows <- dim(MATRIX)[1]
  cols <- dim(MATRIX)[2]
  
  ## STAGE 1  - remove zero rows/columns
  
  #Sorts the binary positions of matrix MAT in terms of row and column degree.
  BMAT <- 1*(MATRIX!=0);
  
  colsum <- colSums(BMAT)
  rowsum <- rowSums(BMAT)
  
  
  JJ <- which(colsum==0);#Find indexes of columns with no interactions
  KK <- which(rowsum==0);#Find indexes of rows with no interactions
  
  index_rows <- 1:rows;
  index_cols <- 1:cols;
  
  #remove these
  sortMAT <- MATRIX;
  if (sum(JJ) > 0) {
    sortMAT <- sortMAT[,-JJ]
    index_cols<-index_cols[-JJ];
  }
  
  if (length(dim(sortMAT)) < 2) {
    sortMAT <- sortMAT[-KK]
    index_rows<-index_rows[-KK]
  } else if (sum(KK) >0) {
    sortMAT <- sortMAT[-KK,]
    index_rows<-index_rows[-KK]
  }
  
  
  ## STAGE 2 -- maximal sorting
  
  if (sortvar==1) {#If want to package matrices in degree ordered way
    
    #BINARY SORTMAT.
    BMAT <- 1*(sortMAT!=0);
    
    #Sort rows and columns in decending order to find zero rows, columns.
    #This is the sort for entire matrix which may contain zero rows and
    #columns.
    if (length(dim(BMAT)) == 2) { #if not a matrix after removing rows and columns, then sortMAT is not changed 
      new<-sort(rowSums(BMAT),decreasing=TRUE,index.return=TRUE)
      new.index_rows <- new$ix
      new<-sort(colSums(BMAT),decreasing=TRUE,index.return=TRUE)
      new.index_cols <-new$ix
      
      index_rows = index_rows[new.index_rows];
      index_cols = index_cols[new.index_cols];
      
      sortMAT <- MATRIX[index_rows,index_cols];
      
      
      ## STAGE 3 - If input matrix not binary
      
      if ((sum(sum( (MATRIX==0) + (MATRIX==1)))!=rows*cols) && (binary!=1)) {
        BMAT <- 1*(sortMAT!=0)
        
        colsum <- colSums(BMAT);
        rowsum <- rowSums(BMAT);
        
        JJ <- unique(colsum);
        KK <- unique(rowsum);
        
        if (length(JJ) < dim(BMAT)[2]) { #If more than 1 column with same degree choose ordering
          
          for (aa in 1:length(JJ)) {
            if (sum(colsum==JJ[aa])>1) { #If more than one occurance
              indexes <- which(colsum==JJ[aa]);
              
              HH <- matrix(0,length(indexes),length(indexes));
              
              for (RR in 1:dim(HH)[1]) {
                for (CC in 1:dim(HH)[2]) {
                  HH[RR,CC] <- sum(sortMAT[,indexes[RR]]>sortMAT[,indexes[CC]]);
                }
              }
              
              
              select <- colSums(HH);
              unis <- unique(select);
              if (length(unis) < length(select)) {
                for (uu in 1:length(unis)) {
                  if (sum(select==unis[uu])>1) {#if more than 1 value of this
                    nextind <- which(select==unis[uu]);
                    select[nextind] <- 0.5*(colSums(sortMAT[,nextind])/(sum(sortMAT[,nextind])));
                  }
                }
                
              }
              
              new <- sort(colSums(HH),index.return=TRUE);
              ind <- new$ix
              
              
              index_cols[indexes] <- index_cols[indexes[ind]];
              
              
              sortMAT[,indexes] <- sortMAT[,indexes[ind]];
              
              
            }
          }
        }
        
        if (length(KK) < dim(BMAT)[1]) { #If more than 1 row with same degree choose ordering
          
          for (aa in 1:length(KK)) {
            if (sum(rowsum==KK[aa])>1) { #If more than one occurance
              indexes <- which(rowsum==KK[aa]);
              
              HH <- matrix(0,length(indexes),length(indexes));
              
              for (RR in 1:dim(HH)[1]) {
                for (CC in 1:dim(HH)[2]) {
                  HH[RR,CC] <- sum(sortMAT[indexes[RR],]>sortMAT[indexes[CC],]);
                }
              }
              
              select <- rowSums(HH);
              unis <- unique(select);
              if (length(unis) < length(select)) {
                for (uu in 1:length(unis)) {
                  if (sum(select==unis[uu])>1) { #if more than 1 value of this
                    nextind <- which(select==unis[uu]);
                    select[nextind] <- 0.5*(rowSums(sortMAT[nextind,])/(sum(sortMAT[nextind,])));
                  }
                }
                
              }
              
              
              new <- sort(colSums(HH),index.return=TRUE);
              ind <- new$ix
              
              index_rows[indexes] <- index_rows[indexes[ind]];
              
              
              sortMAT[indexes,] <- sortMAT[indexes[ind],];
              
              
            }
            
          }
          
          
          
        }
        
        sortMAT <- MATRIX[index_rows,index_cols];
        
        
      }
      
    }
    
    
    
    
  }
  
  
  
  #FINALLY IF WANT TO RETURN A BINARY MATRIX
  if (binary==1) {
    sortMAT <- 1*(sortMAT!=0);
  }
  
  return(list(sortMAT=sortMAT,index_rows=index_rows,index_cols=index_cols))
  
}

DD_null <- function(MATRIX) {#%DD
  # Degreeprobable - Degreeprobable null model
  #Proportionally fills matrix depending on size and degree distribution of
  #rows and columns - as described in Bascompte et al. 2003.However, we do 
  #prevent degenerate matrices being formed - i.e. those which are empty, 
  #scalar or vectors.
  
  #J Bascompte, P Jordano, CJ Melián, JM Olesen. 2003.
  #The nested assembly of plant–animal mutualistic networks.
  #PNAS 100: 9383–9387. (http://dx.doi.org/10.1073/pnas.1633576100)
  MAT <- MATRIX
  #Find matrix dimensions and row and column degrees.
   r<-dim(MAT)[1]
   c<-dim(MAT)[2]
   
   coldegreesprop<-(colSums(MAT))/r
   rowdegreesprop<-(rowSums(MAT))/c
  # 
  TEST<- 1* ( array(runif(r*c), dim=c(r,c)) < 0.5 * (array(rep(coldegreesprop,rep(r,c)), dim=c(r,c)) + array(rep(rowdegreesprop,c),dim=c(r,c)))) 
  # # #sort
    TEST<-sortMATRIX(TEST,1,T)$sortMAT
  return(TEST)
}




###### Helper function to calculate null models for modularity and nestedness ( modified From Munoz et al. 2019. JBI )
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
NullModSen =  function(matrix, rep, nulltype ){
  matrix <-matrix[as.logical(rowSums(matrix != 0)), as.logical(colSums(matrix != 0))]
  
  matrix <- as.matrix.data.frame(matrix)
  ## generate null models
  
  if(nulltype == "curveball"){
    vn <- vegan::nullmodel(matrix, "curveball")
    vns <- simulate(vn, rep)
    null.nest <- sapply(1:rep, function(x) vegan::nestednodf(vns[,,x])$statistic["NODF"])
    null.mod <- sapply(1:rep, function(x) bipartite::LPA_wb_plus(vns[,,x])$modularity)
    
  }
  if(nulltype == "DD") {
    times <- rep
    MAT <- matrix
    vns <- c()
    for(i in 1:rep){vns[[i]] <- DD_null(MAT) }
    
    null.nest <- sapply(1:rep, function(x) vegan::nestednodf(vns[[x]])$statistic["NODF"])
    null.mod <- sapply(1:rep, function(x) bipartite::LPA_wb_plus(vns[[x]])$modularity)
    
  }
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
                       chunkA = chunkA, chunkB = chunkB, nulltype = nulltype, 
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
                         repNull = repNull, 
                         nulltype = nulltype)
  
  
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
                       chunkA = chunkA, chunkB = chunkB, nulltype,
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
calculateMetrix <- function(metaMatrix, quantile, nulltype, repNull = 100){
 
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
  
    for(i in quantile){ 
      
      # Create binary matrices
      
      # subset network based on interaction threshold 
      
      subSet <-  metaMatrix[metaMatrix$quantile >= i,]
      # Make binary matrix
    
      intData2 <- table(subSet$sp1, subSet$sp2)
      intData2[intData2 >= 1] <- 1
      attributes(intData2)$class <- "matrix"
      
      # 
      # # null models
      null1 <- NullModSen(intData2, repNull, nulltype)
      # 
      # # Calculate  z-score Modularity
      zMod.net1 <- (bipartite::LPA_wb_plus(intData2)$modularity - null1[[1]]["NulMod.mean"]) / null1[[1]]["NulMod.sd"]
      # 
      # # Calculate z-score Nestedness 
      # 
       zNes.net1 <- (bipartite::nested(intData2, method = "NODF2") - null1[[1]]["NulNest.mean"]) / null1[[1]]["NulNest.sd"]
      
      # zNes.net1<-  makeNull(matrix,repNull, nulltype = nulltype )
      # zMod.net1<-  makeNullNes(matrix,repNull, nulltype = nulltype )
      # 
      res[[i]] <- list("metrics" = data.frame("NODF Z-score" = zNes.net1,
                             "Q Z-score" = zMod.net1 ,
                             "con" = networklevel(intData2, weighted = F, index = c("connectance"))),
                       "bnet" = subSet )
    }
    


  res2 <- list("metrics" = res, "wnet" = intData)
  return(res2)
}
###########

###########################
###########################

###########################
###########################


### Wrapper function to run the simulations 


RunSimul <- function(Pool, Scenario, replicates = 10, repNull, quantile,nulltype, runParal = T){
  
  
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
                                                                                   nulltype = nulltype,
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
                                                                                 nulltype =nulltype,
                                                                                 intHyp =  Scenario[x,]$intHyp)
    )
    )
    
    
    
  }
  
  
  
  ## Give appropiate names to the resulting object 
  
  
  # names(SIMS) <- paste0(Scenario$EA,"_",
  #                       Scenario$EB, "_",
  #                       Scenario$intHyp, "_",
  #                       Scenario$TraitA, "_", 
  #                       Scenario$TraitB,"_",
  #                       Scenario$sigmaA,"_", 
  #                       Scenario$sigmaB,"_",
  #                       Scenario$chunkA,"_",
  #                       Scenario$chunkB)
  
  
  
  
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
                                                                             nulltype=nulltype,
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
                                                                           nulltype = nulltype,
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

### Make binary 

makeBinar <- function(matrix){
  matrix <-matrix[as.logical(rowSums(matrix != 0)), as.logical(colSums(matrix != 0))]
  matrix <- as.matrix.data.frame(matrix)
  return(matrix)
  
  
}


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
makeNull <- function(matrix,nullSim, nulltype= c("curveball", "DD")){
  matrix <-matrix[as.logical(rowSums(matrix != 0)), as.logical(colSums(matrix != 0))]
  
  matrix <- as.matrix.data.frame(matrix)
  mod <- bipartite::computeModules(matrix)
  
  if(nulltype == "curveball"){
  vn <- vegan::nullmodel(matrix, "curveball")
  vns <- simulate(vn, nullSim)
  null <- sapply(1:nullSim, function(x) bipartite::computeModules(vns[,,x])@likelihood)
  modES <- (mod@likelihood - mean(null))/sd(null)
  }
  if(nulltype == "DD"){
    
    vns <- replicate(nullSim,DD_null(matrix))
    print(vns)
    null <- sapply(1:nullSim, function(x) bipartite::computeModules(vns[[x]])@likelihood)
    
    modES <- (mod@likelihood  - mean(null))/sd(null)
    
  }
  
    return(modES)
    
}




######
# function to calculate nestedness z-scores 
makeNullNes <- function(matrix, nullSim, nulltype = c("curveball", "DD")){
  # make into binary network first
  matrix[matrix > 0] <- 1
  # remove 0s 
  matrix <-matrix[as.logical(rowSums(matrix != 0)), as.logical(colSums(matrix != 0))]
  matrix <- as.matrix.data.frame(matrix)
  nes <- vegan::nestednodf(matrix)$statistic["NODF"]

  if(nulltype == "curveball"){
    vn <- vegan::nullmodel(matrix, "curveball")
    vns <- simulate(vn, nullSim)
    null <- sapply(1:nullSim, function(x) vegan::nestednodf(vns[,,x])$statistic["NODF"])
    modES <- (nes - mean(null))/sd(null)
    
  }
  if(nulltype == "DD"){
    
    vns <- replicate(nullSim,DD_null(matrix))
    print(vns)
    null <- sapply(1:nullSim, function(x) vegan::nestednodf(vns[[x]])$statistic["NODF"])
    modES <- (nes - mean(null))/sd(null)
    
  }
 
  return(modES)}



#######
# CUSTOM FUNCTIONS FOR SECTION 4
#######


# Custom function to geoprocess the raster into cost-distance matrices based on points 

CustomCostDist <- function(raster, centroids ){ 
  centroidPoint <- sp::SpatialPoints(raster::coordinates(centroids[,c(1,2)]))
  ## change raster into transition matrices 
  ## (cost of moving pixel = mean btw pixels, 4 directions), and geocorrect for planar geometries
  tr <- transition( 1/raster, transitionFunction = mean, directions = 4 )
  tr <- geoCorrection( tr, type="c", multpl=FALSE, scl=FALSE)
  ## Get costDistance pairwise matrices
  eDist <- costDistance( tr, centroidPoint )
  eDist <- as.matrix(eDist)
  # Fix row and col names
  rownames(eDist) <- colnames(eDist) <- centroids$site
  return(eDist)
}


### Helper functions 

## Lets define first a series of helpful functions to be able to calculate Environmental Filtering intensity and to try different combinations of methodologies as sensitivity analisis 


## make a helper function that accepts a diversity scalar and returns a list of species randomized, based on the probabilities of the probPool

SpeciesSampler <- function(divVec, ProbPool, pool, PoA = "P"){
  # select one site location
  
  probSite <- ProbPool[rownames(ProbPool) == names(divVec),]
  # make a distribution of n sites to pick a species from it 
  sitesToSample <-replicate(divVec,names(sample(probSite, 1,prob = probSite)))
  
  # iterate over this sample and pick a random species from the pool, based on the site to match 
  if(PoA == "P"){
    spSampled <- sapply(1:length(sitesToSample), function(x) sample(unique(pool$plantCode[pool$site == sitesToSample[x]]), 1))
  }
  
  if(PoA == "A"){
    spSampled<-  sapply(1:length(sitesToSample), function(x) sample(unique(pool$animalCode[pool$site == sitesToSample[x]]), 1))
  }
  
  return(spSampled)
}


## make a funtion that feeds on "SpeciesSampler function over diversity vector, change n in `replicate` to change the number of replicates


NullTraitDiv <- function(divVec,ProbPool, pool, nrep, PoA ){
  NullSites <- sapply(1:length(divVec),
                      function(x)
                        replicate(nrep,
                                  SpeciesSampler(divVec[x], ProbPool, pool, PoA )))
  
  names(NullSites) <- names(divVec)
  return(NullSites)
  
}

# make helper function that accepts NullTr and iterates based on site name and returns the distribution of range of null models 

ZdisCalR <- function(sitename, NullTr, PoA = "P"){
  
  toAg <- NullTr[NullTr$L1 == sitename,]
  if(PoA == "P"){ 
    AgRang <- aggregate(toAg$RQLp, list(toAg$Var2), range)
  }
  if(PoA == "A"){
    AgRang <- aggregate(toAg$RQLa, list(toAg$Var2), range)
  }
  Zdis <- abs(AgRang$x[,1]-AgRang$x[,2])
  return(Zdis)
}

# make helper function that accepts NullTr and iterates based on site name and returns the distribution of sd of null models 

ZdisCalSD <- function(sitename, NullTr, PoA = "P"){
  toAg <- NullTr[NullTr$L1 == sitename,]
  if(PoA == "P"){
    AgRang <- aggregate(toAg$RQLp, list(toAg$Var2), sd)
  }
  if(PoA == "A"){
    AgRang <- aggregate(toAg$RQLa, list(toAg$Var2), sd)
  }
  return(AgRang$x)
}



## Defining the main function to test for environmental filtering, feeding on all helper functions above.

# Since we have now defined the series of helper functions, now let's bind them all together into a logical framework function that allow to try various settings for the sensitivity analisis 

#####
# Test4EF: Function to test for environmental filtering based on different choices of process based or random based species pools 
##### 
test4EF <- function(pool,
                    traitName, 
                    nrep = 100,
                    side = c("P","A"),
                    RoP = c("Prob","Random"), 
                    sdOrRange = c("sd", "range"), 
                    ProbPool = ProbPool){ 
  
  
  
  ## calculate observed sd or range
  if(sdOrRange == "sd"){
    # calculate sd 
    obsDis <- aggregate(pool[[traitName]], list(pool$site), sd)$x
    names(obsDis) <- levels(pool$site)
    
    
  }
  if(sdOrRange == "range"){
    # calculate range
    obsMin <- aggregate(pool[[traitName]], list(pool$site), min)$x
    obsMax <- aggregate(pool[[traitName]], list(pool$site), max)$x
    obsDis <- abs(obsMin-obsMax)
    names(obsDis) <- levels(pool$site)
    
    
  }
  
  ## Create null models (random or process based)
  if(RoP == "Random"){
    # create randomized null models 
    if(side == "P"){
      # get diversity of plants
      divVec <- colSums(ifelse(table(pool$plantCode,pool$site) > 1, 1,0))
      if(sdOrRange == "sd"){
        # make null model 
        nullDis <- lapply(divVec, function(x) replicate(nrep,
                                                        sd(unlist(pool[traitName])[
                                                          match(sample(unique(pool$plantCode), x), 
                                                                pool$plantCode)] )))
      }
      if(sdOrRange == "range"){
        nullRang <- lapply(divVec, function(x) replicate(nrep,
                                                         range(unlist(pool[traitName])[
                                                           match(sample(unique(pool$plantCode), x), 
                                                                 pool$plantCode)] )))
        
        
        nullDis <- lapply(nullRang, function(x) diff(x))
      } 
    }
    if(side == "A"){
      # get diversity of animals
      divVec <- colSums(ifelse(table(pool$animalCode,pool$site) > 1, 1,0))
      # make null model 
      if(sdOrRange == "sd"){
        # make null model 
        nullDis <- lapply(divVec, function(x) replicate(nrep,
                                                        sd(unlist(pool[traitName])[
                                                          match(sample(unique(pool$animalCode), x), 
                                                                pool$animalCode)] )))
      }
      if(sdOrRange == "range"){
        nullRang <- lapply(divVec, function(x) replicate(nrep,
                                                         range(unlist(pool[traitName])[
                                                           match(sample(unique(pool$animalCode), x), 
                                                                 pool$animalCode)] )))
        
        
        nullDis <- lapply(nullRang, function(x) diff(x))
      } 
      
    }
    
  }
  if(RoP == "Prob"){
    if(side == "P"){
      divVecPlant <- colSums(ifelse(table(pool$plantCode,pool$site) > 1, 1,0))
      
      
      # plants
      NullTrP <- NullTraitDiv(divVecPlant,ProbPool,pool, nrep , "P")
      NullTrP <- reshape2::melt(NullTrP)
      NullTrP$RQLp <- pool$RQLp[match(NullTrP$value,pool$plantCode )]
      if(sdOrRange == "range"){
        nullDis <- lapply(1:length(unique(NullTrP$L1)), 
                          function(x)
                            ZdisCalR(unique(NullTrP$L1)[x],NullTrP, "P" ))
        
      }
      if(sdOrRange == "sd"){
        
        nullDis <- lapply(1:length(unique(NullTrP$L1)), 
                          function(x)
                            ZdisCalSD(unique(NullTrP$L1)[x],NullTrP, "P" ))
      }
      
      
    }
    if(side == "A"){
      divVecAn <- colSums(ifelse(table(pool$animalCode,pool$site) > 1, 1,0))
      # animals 
      NullTrA <- NullTraitDiv(divVecAn,ProbPool,pool, nrep, "A" )
      NullTrA <- reshape2::melt(NullTrA)
      NullTrA$RQLa <- pool$RQLa[match(NullTrA$value,pool$animalCode )]
      
      if(sdOrRange == "range"){
        nullDis <- lapply(1:length(unique(NullTrA$L1)), 
                          function(x)
                            ZdisCalR(unique(NullTrA$L1)[x],NullTrA, "A" ))
      }
      
      if(sdOrRange == "sd"){
        nullDis <- lapply(1:length(unique(NullTrA$L1)), 
                          function(x)
                            ZdisCalSD(unique(NullTrA$L1)[x],NullTrA, "A" ))
      }
      
    }
  }
  
  
  # calculate sd and mean from null dist 
  sdNull <- sapply(nullDis, function(x) sd(x))
  xNull <- sapply(nullDis, function(x) mean(x))
  
  # compute ses
  SES <- (obsDis-xNull )/sdNull
  
  return(SES)
  
}

###############

###########
# Crop of test4EF function to retrieve species names and not calculate effect sizes
############
SamplePoolComm <- function(pool,
                           traitName, 
                           nrep = 100,
                           side = c("P","A"),
                           ProbPool = ProbPool){
  
  if(side == "P"){
    divVecPlant <- colSums(ifelse(table(pool$plantCode,pool$site) > 1, 1,0))
    # plants
    NullTr <- NullTraitDiv(divVecPlant,ProbPool,pool, nrep , "P")
    NullTr <- reshape2::melt(NullTr)
    NullTr$RQLp <- pool$RQLp[match(NullTr$value,pool$plantCode )]
    
  }
  if(side == "A"){
    divVecAn <- colSums(ifelse(table(pool$animalCode,pool$site) > 1, 1,0))
    # animals 
    NullTr <- NullTraitDiv(divVecAn,ProbPool,pool, nrep, "A" )
    NullTr <- reshape2::melt(NullTr)
    NullTr$RQLa <- pool$RQLa[match(NullTr$value,pool$animalCode )]
    
  }
  
  names(NullTr) <- c("SpeciesSamID","PoolID", "SpecName", "Site",  names(NullTr)[5])
  return(NullTr)
}


######
# NullNetMod: function to generate null networks from a given pool and a site-selection probability metrics. Options are to toogle off between filtered and non-filtered communities.  
## Interactions are mantained neutral 
########


NullNetMod <- function(pool, ProbPool, repPool, NPlant = T, NAnimal = T, SiteName, netRep, nulltype){
  
  atab <- table(pool$animalCode, pool$site)
  atab[atab>1] <- 1
  
  ptab <- table(pool$plantCode, pool$site)
  ptab[ptab>1] <- 1
  
  
  AnimalRich <- colSums(atab)
  PlantRich <- colSums(ptab)
  
  
  # if the communities are neutral 
  
  if(NAnimal == T){
    AnimalSide <- SamplePoolComm(pool = pool, 
                                 nrep = repPool, 
                                 side = "A", 
                                 ProbPool = ProbPool )
    # simulate pseudo-abundances at the pool level
    animals <- table(droplevels(AnimalSide$SpecName[AnimalSide$Site == SiteName]))
    # standardize into relative abundances
    animals <- animals/sum(animals)
    # sample local communtiies
    animals <- sample(animals, AnimalRich[SiteName], prob = animals)
    # standarize to relative abundances
    animals <- animals/sum(animals)
    
  }else{ 
    # if communities are the observed 
    animals <- table((pool$animalCode[pool$site == SiteName]))
    animals <- animals/sum(animals)
  }
  
  if(NPlant == T){
    
    PlantSide <- SamplePoolComm(pool = pool, 
                                nrep = repPool, 
                                side = "P", 
                                ProbPool = ProbPool )
    
    # simulate pseudo-abundances at the pool level
    plants <- table((PlantSide$SpecName[PlantSide$Site == SiteName]))
    # standardize into relative abundances
    plants <- plants/sum(plants)
    # sample local communtiies
    plants <- sample(plants, PlantRich[SiteName], prob = plants)
    # standarize to relative abundances
    plants <- plants/sum(plants)
  }else{
    plants <- table(pool$plantCode[pool$site == SiteName])
    plants <- plants/sum(plants)
  }
  

  # recreate the metanetwork
  int <- data.frame(expand.grid(plants,animals), expand.grid(names(plants),names(animals)))
  names(int) <- c("relAbPlant", "relAbAni", "PlantID", "AnimalID")
  int$probInt <- int$relAbPlant * int$relAbAni
  
  
  
  # Get the interaction richness per site
  intDiv <- table(pool$intID, pool$site)
  intDiv[intDiv>1] <- 1
  intDiv <- colSums(intDiv)
  
  # Randomly sample the metanetwork based on the interaction richness observed at site (To simulate neutral interactions)
  netSam <- int[sample(1:length(int[,1]), intDiv[SiteName], prob = int$probInt  ),]
    
  # # Calculate Modularity Z-scores based on given number of replicates
  # print("Quantyfing modularity in simulated network")
  # modZ <- makeNull(matrix = table(netSam$PlantID, netSam$AnimalID),
  #                  nullSim = netRep, 
  #                  nulltype = nulltype )
  # cat("\014")
  # # Calculate Nestedness Z-scores based on given number of replicates
  # print("Quantyfing nestedness in simulated network")
  # nesZ <- makeNullNes(table(netSam$PlantID, netSam$AnimalID),netRep, nulltype = nulltype)
  # # erase screen
  # cat("\014")
  # 
  
  intData2 <- table(netSam$PlantID, netSam$AnimalID)
  intData2[intData2 >= 1] <- 1
  attributes(intData2)$class <- "matrix"
  
  # # null models
  null1 <- NullModSen(intData2, netRep, nulltype)
  # 
  # # Calculate  z-score Modularity
  zMod.net1 <- (bipartite::LPA_wb_plus(intData2)$modularity - null1[[1]]["NulMod.mean"]) / null1[[1]]["NulMod.sd"]
  # 
  # # Calculate z-score Nestedness 
  # 
  zNes.net1 <- (bipartite::nested(intData2, method = "NODF2") - null1[[1]]["NulNest.mean"]) / null1[[1]]["NulNest.sd"]
  
  
  
  return(c("Q"= zMod.net1, "NODF" = zNes.net1))
}







