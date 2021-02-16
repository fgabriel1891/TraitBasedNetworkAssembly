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

#####
# Custom functions 
#####


#' Create a Species Pool 
#' 
#' This function creates regional species pool of a given abundance, richness, and rank-distribution. 
#' The function returns 2 regional pools with the same parameters, one for each trophic level
#' 
#' @param Jpool Total number of individuals in the pool 
#' @param J Total number of species in the pool 
#' @param poolShape Rank-abundance curve shape, options = uniform or log-series

CreateSpPool <- function(Jpool, J, poolShape = c("uniform", "log-series")) {
  
  if (poolShape == "uniform") { 
    
    # With uniform trait values in the species pool
    poola <- data.frame(cbind(1:Jpool[1], 
                              rep(1:J[1], 
                                  Jpool[1]/J[1])))
    poola <- poola[order(poola$X2),]
    poola$X3 <- c(sapply(runif(J[1],0,1),
                         function(x) rnorm(10, x, 0.001)))
    colnames(poola) = c("ind"  ,"sp",   "trait")
    
    
    poolb <- data.frame(cbind(1:Jpool[2], 
                              rep(1:J[2],
                                  Jpool[2]/J[2])))
    poolb <- poolb[order(poolb$X2),]
    poolb$X3 <- c(sapply(runif(J[2],0,1), 
                         function(x) rnorm(10, x, 0.001)))
    colnames(poolb) = c("ind"  ,"sp",   "trait")
    
    
  }
  
  if(poolShape == "log-series"){
    
    
    poola <- coalesc(Jpool[1], 
                     m = 0.5,
                     theta = 50)$pool
    poolb <- coalesc(Jpool[2],
                     m = 0.5, theta = 50)$pool
    
  }
  
  
  
  
  return(list( "poolA" = poola, 
               "poolB" = poolb))
  
}




#' RunSimul: Wrapper function to run the simulations 
#' 
#' This function takes the simulation paramenter values and pass them to the function to create communities and simulate interactions
#' This functions returns a nested list object with branches being the different assembly scenarios and nodes in branches the independent assembly replicates
#' 
#' @param Pool Species pool oject
#' @param Scenario Object with the combinations of assembly scenarios
#' @param replicates number of independent replicates 
#' @param repNull number of reshuffles for z-scores 
#' @param quantile Quantile to prune interactions 
#' @param nulltype Specify the null model for the reshuffling step 
#' @param runParal T/F if run simulations in parallel



RunSimul <- function(Pool, 
                     Scenario, 
                     replicates = 10,
                     repNull, 
                     quantile,
                     nulltype, 
                     runParal = T){
  
  
  if (runParal == T){
    
    
    SIMS <- sfLapply(1:length(Scenario$EA), 
                     function(x) 
                       replicate(replicates , # number of replicates per scenario
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
    
    SIMS <- lapply(1:length(Scenario$EA), 
                   function(x) 
                     replicate(replicates , # number of replicates per scenario
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


#' SimulateEF 
#' 
#' Wrapper function to simulate community assembly mechanisms, simulate interactions and calculate network metrics, many parameters pass on different functions of to the `coalesc` function from the `eccolottery` package
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

simulateEF <- function(TraitA, 
                       TraitB,
                       mA, 
                       mB,
                       JA, 
                       JB, 
                       poolA, 
                       poolB, 
                       EF.A = c("ND","SEF", "DEF"), 
                       EF.B = c("ND","SEF", "DEF"), 
                       sigmaA, 
                       sigmaB,
                       repNull = 10, 
                       quantile = 10, 
                       nulltype = nulltype, 
                       intHyp = c("MM", "FL", "NL")) {
  
  
  comA <- SimulateSelCom(TraitA, 
                         mA, 
                         JA,
                         poolA,
                         EF = EF.A, 
                         sigma = sigmaA)
  comB <- SimulateSelCom(TraitB, 
                         mB, 
                         JB, 
                         poolB,
                         EF = EF.B, 
                         sigma = sigmaB)
  
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



#' SimulateSelCom
#' 
#' Helper function to call ecolottery package functions to 
#' simulate assembly processes in ecological communities
#' All parameters come from wrapper function "simulateEF" 
#' @param sigma Degree of spacing from the species pool. Modulating it would mean the "strenght" of the environmental filtering process
#' @param Trait Trait mean of trophic level
#' @param m Migration rate  
#' @param J Number of species in the regional pool
#' @param EF  assembly process  
#' @param pool Species pool object 


SimulateSelCom  <- function(Trait,
                            m,
                            J, 
                            pool, 
                            EF = c("ND","SEF", "DEF"), 
                            sigma = sigma){
  
  
  # Simulate environmental filtering - Stabilizing filtering
  
  if (EF == "SEF") { 
    
    sigma <- sigma
    
    comA <- coalesc(J, 
                    m, 
                    filt = function(x)
                      filt_gaussian(Trait, x),
                    pool = pool)
    
    
  }
  
  # Simulate environmental filtering - Directional filtering 
  
  if (EF == "DEF") {
    
    
    comA <- coalesc(J,
                    m, 
                    filt = function(x) 
                      1 - min(x,1), 
                    pool = pool)
    
    
  }
  
  # Simulate neutral dynamics 
  
  if (EF == "ND") { 
    
    comA <- coalesc(J,
                    m,
                    pool = pool) 
    
    
  }
  
  
  return(comA)
  
}


#' filt_gaussian
#' 
#'  Function to filter a distribution based on a gaussian distribution with sd=sigma
#' to be called within SimulateSelCom

filt_gaussian <- function(t,x) exp(-(x-t)^2/(2*sigma^2))

#' createMetaMatrix:
#' 
#'function calculate a meta-matrix of resource-consumer interactions from the filtered communities 
#'
#' @param Resource: Resource communtity
#' @param Consumer: Consumer communtity
#' @param intHyp: Interaction hypothesis 
#' @return a matrix of pairwise interactions in long format

createMetaMatrix <- function(Resource, 
                             Consumer, 
                             intHyp = c("MM", "FL", "NL"))
{ 
  
  ## Coercing column names
  names(Resource) <- c("id", "sp", "trait")
  names(Consumer) <- c("id", "sp", "trait")
  
  ### Simulating interactions 
  print("Simulating interactions")
  
  # All pairwise combinations
  metaMatrix <- expand.grid(unique(Resource$id),
                            unique(Consumer$id))
  
  # Add relative abundances after filtering (filtering probabilities)
  
  # Resources
  metaMatrix$relAbRes <- prop.table(table(Resource$id))[match(metaMatrix$Var1,
                                                              names(prop.table(table(Resource$id))))]
  
  metaMatrix$trait1 <- Resource$trait[match(metaMatrix$Var1, 
                                            Resource$id)]
  metaMatrix$sp1 <- Resource$sp[match(metaMatrix$Var1,
                                      Resource$id)]
  
  # Consumers 
  metaMatrix$relAbCons <- prop.table(table(Consumer$id))[match(metaMatrix$Var2, 
                                                               names(prop.table(table(Consumer$id))))]
  metaMatrix$trait2 <- Consumer$trait[match(metaMatrix$Var2, 
                                            Consumer$id)]
  metaMatrix$sp2 <- Consumer$sp[match(metaMatrix$Var2,
                                      Consumer$id)]
  
  # Erase non interacting species 
  metaMatrix <- na.omit(metaMatrix) 
  
  metaMatrix$intProb <- metaMatrix$relAbRes*metaMatrix$relAbCons ### As a probability of abundance (i.e. int. probability equals the product of both relative abundances)
  metaMatrix$intProb <- vegan::decostand(metaMatrix$intProb, "max")
  #metaMatrix$intProb <- metaMatrix$relAbPl^2/(metaMatrix$relAbPl^2+metaMatrix$relAbPol^2) # disregard this part and not uncomment. Was a tryout 
  
  # Change interaction probabilities based on interaction assembly hypothesis 
  
  if (intHyp == "FL"){
    ## Forbidden links 1-trait <= 0 --> int prob = 0 
    metaMatrix$intProb <- ifelse(metaMatrix$trait1-metaMatrix$trait2 >= 0,
                                 metaMatrix$intProb * 1,
                                 metaMatrix$intProb * 0)
  } 
  
  if (intHyp == "MM") {
    ## Morphological matching, interaction depends on the frequency of differences between traits 
    metaMatrix$intProb <- metaMatrix$intProb*(1- abs(metaMatrix$trait1-metaMatrix$trait2)
                                              /max(abs(metaMatrix$trait1-metaMatrix$trait2)
                                              ))
    metaMatrix$intProb <- vegan::decostand(metaMatrix$intProb, "max")
    
    
  }
  
  if (intHyp == "NL") {
    # Neutral assembly of interactions, probability of interactions depend solely on abundance (rel.abundance of A * rel.abundance of B)
    metaMatrix$intProb <- metaMatrix$intProb * 1
  }
  
  
  return(metaMatrix)
  
}



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
    null1 <- NullModSen(intData2, 
                        repNull, 
                        nulltype)
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
                                            "con" = networklevel(intData2,
                                                                 weighted = F, 
                                                                 index = c("connectance"))),
                     "bnet" = subSet )
  }
  
  
  
  res2 <- list("metrics" = res, "wnet" = intData)
  return(res2)
}



#' NullModSen: 
#'
#'Helper function to calculate null models for modularity and nestedness ( modified From Munoz et al. 2019. JBI )
#'
#' @param matrix A network organized as a matrix in with species A as columns  Species B as rows.
#' @param rep  integer specifying the number of replicates for the null model matrices, i.e. shuffles of the web, passes to bipartite::nullmodel
#' @return a list with the mean and sd of calculated modularity and nestedness for the null model.
#' @examples NullModSen(matrix,rep)

NullModSen =  function(matrix, rep, nulltype ){
  # remove non interacting rows and columns
  matrix <-matrix[as.logical(rowSums(matrix != 0)), 
                  as.logical(colSums(matrix != 0))]
  
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
    
    null.nest <- sapply(1:rep,
                        function(x) vegan::nestednodf(vns[[x]])$statistic["NODF"])
    null.mod <- sapply(1:rep,
                       function(x) bipartite::LPA_wb_plus(vns[[x]])$modularity)
    
  }
  null = list(c("NulMod" = c("mean" = mean(null.mod),
                             "sd" = sd(null.mod)),
                "NulNest" = c("mean" = mean(null.nest), 
                              "sd" = sd(null.nest))))
  return(null)
  
}



#' DD_null
#' 
#' 
#' Degreeprobable - Degreeprobable null model
#' Proportionally fills matrix depending on size and degree distribution of
#' rows and columns - as described in Bascompte et al. 2003.However, we do 
#' prevent degenerate matrices being formed - i.e. those which are empty, 
#' scalar or vectors.
#' J Bascompte, P Jordano, CJ Melián, JM Olesen. 2003.
#' The nested assembly of plant–animal mutualistic networks.
#' PNAS 100: 9383–9387. (http://dx.doi.org/10.1073/pnas.1633576100)
#' @param MATRIX a matrix of species interactions 

DD_null <- function(MATRIX) {
  
  MAT <- MATRIX
  #Find matrix dimensions and row and column degrees.
  r<-dim(MAT)[1]
  c<-dim(MAT)[2]
  
  coldegreesprop<-(colSums(MAT))/r
  rowdegreesprop<-(rowSums(MAT))/c
  # 
  TEST<- 1* ( array(runif(r*c), dim=c(r,c)) < 0.5 *
                (array(rep(coldegreesprop,rep(r,c)),
                       dim=c(r,c)) + 
                   array(rep(rowdegreesprop,c),dim=c(r,c)))) 
  # # #sort
  TEST<-sortMATRIX(TEST,1,T)$sortMAT
  return(TEST)
}


#' sortMATRIX
#' 
#' Degreeprobable - Degreeprobable null model
#' Proportionally fills matrix depending on size and degree distribution of
#' rows and columns - as described in Bascompte et al. 2003.However, we do 
#' prevent degenerate matrices being formed - i.e. those which are empty, 
#' scalar or vectors.
#' J Bascompte, P Jordano, CJ Melián, JM Olesen. 2003.
#' The nested assembly of plant–animal mutualistic networks.
#' PNAS 100: 9383–9387. (http://dx.doi.org/10.1073/pnas.1633576100)

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
  
  return(list(sortMAT=sortMAT,
              index_rows=index_rows,
              index_cols=index_cols))
  
}





