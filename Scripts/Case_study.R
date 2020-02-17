### Proposed implementation for the case-study

## Gabriel Muñoz 
# January 2020

source("scripts/CustomFunctions.R") # path might change 

############################################################
############################################################
## 1) define paramter objects and create combinations of assembly scenarios to simulate
############################################################
############################################################



# SD of trait values in the local community and 
# in the regional species pool (environmental filtering)

## Steps 1 + 2
# bind a pool of interactions for the whole mountain. 
pool <- data.frame(
        dataSet1$N,
        dataSet1$pT[match(dataSet1$N$plantCode, dataSet1$pT$plantCode),][,c(2:5)],
        dataSet1$aT[match(dataSet1$N$animalCode, dataSet1$aT$animalCode),][,c(2:5)],
        "RQLp" = intRQL$RLQan[match(dataSet1$N$animalCode, intRQL$animalCode)],
        "RQLa" = intRQL$RLQan[match(dataSet1$N$plantCode, intRQL$plantCode)])

# define function to test for environmental filtering 
test4EF <- function(pool, traitName, side = "P"){ 
        # make trait site elevation table
        tb <- table(unlist(pool[traitName]), 
                    pool$site,
                    pool$elevation)
        # calculate sd 
        obsSD <- rowSums(sapply(1:3,
                                function(x) apply(tb[,,x],2, sd)))
        # create null models 
        if(side == "P"){
                
                # get diversity of plants
                divVec <- colSums(ifelse(table(pool$plantCode,pool$site) > 1, 1,0))
                # make null model 
                nullDis <- lapply(divVec, function(x) replicate(100,
                                                                sd(unlist(pool[traitName])[
                                                                        match(sample(unique(pool$plantCode), x), 
                                                                              pool$plantCode)] )))
        } 
        if(side == "A"){
                # get diversity of animals
                divVec <- colSums(ifelse(table(pool$animalCode,pool$site) > 1, 1,0))
                # make null model 
                nullDis <- lapply(divVec, function(x) replicate(100,
                                                                sd(unlist(pool[traitName])[
                                                                        match(sample(unique(pool$animalCode), x), 
                                                                              pool$animalCode)] )))
                
                
        }
        else{print("revise function parameters")}
        # calculate sd and mean from null dist 
        sdNull <- sapply(nullDis, function(x) sd(x))
        xNull <- sapply(nullDis, function(x) mean(x))
        # compute ses
        SES <- (obsSD-xNull )/sdNull
        
        return(SES)
        
}

# test for EF at each plant trait
EFPlant <- sapply(names(pool)[c(6:9, 14)], function(x) test4EF(pool, traitName = x, side = "P"))
EFPlant <-  data.frame(EFPlant)
EFPlant$elev <- pool$elevation[match(rownames(EFPlant), pool$site)]
EFPlant$site <- pool$site[match(rownames(EFPlant), pool$site)]



# test for EF at each animal trait

EFAnimal <- sapply(names(pool)[c(10:13, 15)], function(x) test4EF(pool, traitName = x, side = "A"))
EFAnimal <-  data.frame(EFAnimal)
EFAnimal$elev <- pool$elevation[match(rownames(EFAnimal), pool$site)]
EFAnimal$site <- pool$site[match(rownames(EFAnimal), pool$site)]

# Examine e. filt effect sizes. 


plot(log1p(EFAnimal[,5])~EFAnimal$elev, ylim = c(0.001,8),
     ylab = "SES RLQ (log)", xlab = "Elevation", 
     col = "skyblue",
     cex = 2, pch = 16)
abline(h=0)
abline(h=c(1.96,-1.96), lty = 2)
points(log1p(EFAnimal[,5])~EFAnimal$elev,cex=2)
points(log1p(EFPlant[,5])~EFPlant$elev, 
       cex = 2, pch = 16, col = "firebrick")
points(log1p(EFPlant[,5])~EFPlant$elev,cex = 2)



## simulate a species pool with uniform distribution of traits 

simPool = CreateSpPool(c(1020, 1270), c(102,127), "uniform")

## 1) How do we assign relative abundances? do we assume more generalist are also more abundant? 
## all species have the same relative abundances 

##trait distribution in the simulated pools match the distribution of traits in the metanetwork
simPool$poolA$trait <- sapply(simPool$poolA$sp, function(x) RLQ$mR$NorS1[x])
simPool$poolB$trait <- sapply(simPool$poolB$sp, function(x) RLQ$mQ$NorS1[x])
hist(simPool$poolB$trait)



customNetw <- function(simPool,NassA = T,NassB = T,Ja, Jb, intHyp = "NL", funA, funB){
        # Coalesce communities
        
        if(NassA == T){ 
                sa <- ecolottery::coalesc(J = Ja, m = 0.5,  
                                          pool = simPool$poolA)
        }
        if(NassB == T){ 
                sb <- ecolottery::coalesc(J = Jb, m = 0.5,  
                                          pool = simPool$poolB)
        }else{ 
                sa <- ecolottery::coalesc(J = Ja, m = 0.5,  
                                          filt = function(x) 
                                                  ifelse(is.na(funA(x)), 0, funA(x)),
                                          pool = simPool$poolA)
                sb <- ecolottery::coalesc(J = Jb, m = 0.5,  
                                          filt = function(x) 
                                                  ifelse(is.na(funB(x)), 0, funB(x)),
                                          pool = simPool$poolB)     
                
                
        }
        
        ta <- table(sa$com$ind)/max(table(sa$com$ind))
        tb <- table(sb$com$ind)/max(table(sb$com$ind))
        id <- expand.grid(rownames(ta), rownames(tb))
        int <- expand.grid(ta,tb)
        
        int <- data.frame(id, int)
        int$pint <- int$Var1.1* int$Var2.1
        names(int) <-c("a", "b", "ra", "rb", "pab")
        int$traitA <- sa$com$trait[match(int$a, sa$com$ind)]
        int$traitB <- sb$com$trait[match(int$b, sb$com$ind)]
        int$spA <- sa$com$sp[match(int$a, sa$com$ind)]
        int$spB <- sb$com$sp[match(int$b, sb$com$ind)]
        if(intHyp == "FL"){
                ## Forbidden links 1-trait <= 0 --> int prob = 0 
                int$pab <- ifelse(int$traitA-int$traitB >= 0, int$pab * 1, int$pab * 0)
        } 
        if(intHyp == "MM") {
                ## Morphological matching, interaction depends on the frequency of differences between traits 
                int$pab <- 1- abs(int$traitA-int$traitB)/max(abs(int$traitA-int$traitB))
                
        }
        if(intHyp == "NL"){ 
                int$pab <- int$pab }
        return(int)
        
        
        
        
}

# calculate delta network structures
deltaNetSt <- function(pool, simulNet, site){
        # define observed and null netowkrs 
        obs <- pool[pool$site == site,]
        ln <- length(unique(obs$intID))
        print(ln)
        null  <- simulNet[sample(1:length(simulNet$a), ln, prob = simulNet$pab),]

        
        # create matrix 
        nullNet <- xtabs(pab ~ spA + spB, null)
        obsNet <- xtabs(frequency ~ plantCode + animalCode, obs)
        
        # normalize probabilities
        nullNet <- nullNet/max(nullNet)
        # calculate modularity
        nulMod <- makeNull(nullNet, 10)
        obsMod <- makeNull(obsNet, 10)  
        deltaMod <- abs(obsMod-nulMod)
        # calculate nestedness 
        nulNes <- makeNullNes(nullNet, 10)
        obsNes <- makeNullNes(obsNet, 10)  
        deltaNes <- abs(obsNes-nulNes)
        
        return(c(deltaMod, deltaNes))
        
        
}



# aproximate probability functions to sample individuals from the pool
funA <- lapply(unique(dataSet1$N$site), 
               function(x) approxfun(density(intRQL$RLQpla[intRQL$site == x])))
names(funA) <- unique(dataSet1$N$site)
funB <- lapply(unique(dataSet1$N$site), 
               function(x) approxfun(density(intRQL$RLQan[intRQL$site == x])))
names(funB) <- unique(dataSet1$N$site)


####################
#AnBoIn
####################
AnBoIn <- c()
delt_AnBoIn <-c()
for(i in 1:length(unique(dataSet1$N$site))) { 
        
        v <- unique(dataSet1$N$site)
        AnBoIn[[i]] <- customNetw(simPool, 120,124,
                                  Nass = T, NassB = F,
                                  intHyp = "NL",
                                  funA[[i]], funB[[i]])
        print(i)
        delt_AnBoIn[[i]] <- deltaNetSt(pool, AnBoIn[[i]], v[i])
        
        }


names(delt_AnBoIn) <- unique(dataSet1$N$site)
reshape2::melt(delt_AnBoIn)


res_AnBoIn <- data.frame("mod" = sapply(1:6, function(x) delt_AnBoIn[[x]][1]), 
                  "nes" = sapply(1:6, function(x) delt_AnBoIn[[x]][2]),
                  "site" = unique(dataSet1$N$site),
                  "RQLa" = EFAnimal$RQLa[match(unique(dataSet1$N$site),EFAnimal$site) ],
                  "RQLp" = EFPlant$RQLp[match(unique(dataSet1$N$site),EFPlant$site) ] )


####################
# AoBoIn
####################
AoBoIn <- c()
delt_AoBoIn <-c()
for(i in 1:length(unique(dataSet1$N$site))) { 
        
        v <- unique(dataSet1$N$site)
        AoBoIn[[i]] <- customNetw(simPool, 120,124,
                                  Nass = F, NassB = F,
                                  intHyp = "NL",
                                  funA[[i]], funB[[i]])
        print(v[i])
        delt_AoBoIn[[i]] <- deltaNetSt(pool, AoBoIn[[i]], v[i])
        
}


names(delt_AoBoIn) <- unique(dataSet1$N$site)
reshape2::melt(delt_AoBoIn)


res_AoBoIn <- data.frame("mod" = sapply(1:6, function(x) delt_AoBoIn[[x]][1]), 
                         "nes" = sapply(1:6, function(x) delt_AoBoIn[[x]][2]),
                         "site" = unique(dataSet1$N$site),
                         "RQLa" = EFAnimal$RQLa[match(unique(dataSet1$N$site),EFAnimal$site) ],
                         "RQLp" = EFPlant$RQLp[match(unique(dataSet1$N$site),EFPlant$site) ] )

####################
####################
#AnBnIn
####################
AnBnIn <- c()
delt_AnBnIn <-c()
for(i in 1:length(unique(dataSet1$N$site))) { 
        
        v <- unique(dataSet1$N$site)
        AnBnIn[[i]] <- customNetw(simPool, 120,124,
                                  Nass = T, NassB = F,
                                  intHyp = "NL",
                                  funA[[i]], funB[[i]])
        print(v[i])
        delt_AnBnIn[[i]] <- deltaNetSt(pool, AnBnIn[[i]], v[i])
        
}


names(delt_AnBnIn) <- unique(dataSet1$N$site)
reshape2::melt(delt_AnBnIn)


res_AnBnIn <- data.frame("mod" = sapply(1:6, function(x) delt_AnBnIn[[x]][1]), 
                         "nes" = sapply(1:6, function(x) delt_AnBnIn[[x]][2]),
                         "site" = unique(dataSet1$N$site),
                         "RQLa" = EFAnimal$RQLa[match(unique(dataSet1$N$site),EFAnimal$site) ],
                         "RQLp" = EFPlant$RQLp[match(unique(dataSet1$N$site),EFPlant$site) ] )

####################
####################
####################
#AnBnIn
####################
AoBnIn <- c()
delt_AoBnIn <-c()
for(i in 1:length(unique(dataSet1$N$site))) { 
        
        v <- unique(dataSet1$N$site)
        AoBnIn[[i]] <- customNetw(simPool, 120,124,
                                  Nass = T, NassB = F,
                                  intHyp = "NL",
                                  funA[[i]], funB[[i]])
        print(v[i])
        delt_AoBnIn[[i]] <- deltaNetSt(pool, AoBnIn[[i]], v[i])
        
}


names(delt_AoBnIn) <- unique(dataSet1$N$site)
reshape2::melt(delt_AoBnIn)


res_AoBnIn <- data.frame("mod" = sapply(1:6, function(x) delt_AoBnIn[[x]][1]), 
                         "nes" = sapply(1:6, function(x) delt_AoBnIn[[x]][2]),
                         "site" = unique(dataSet1$N$site),
                         "RQLa" = EFAnimal$RQLa[match(unique(dataSet1$N$site),EFAnimal$site) ],
                         "RQLp" = EFPlant$RQLp[match(unique(dataSet1$N$site),EFPlant$site) ] )

####################








par(mfrow = c(2,2))

plot(res_AnBnIn$nes~res_AnBnIn$mod,
     xlim = c(0,20), 
     ylim = c(0,4),
     main = "res_AnBnIn",
     ylab = "∆ Nestedness", xlab = "∆ Modularity",
     pch = 16, cex = log(res_AnBnIn$RQLa), col = scales::alpha("firebrick",0.5))
points(res_AnBnIn$nes~res_AnBnIn$mod,
       pch = 16, cex = log(res_AnBnIn$RQLp), col = scales::alpha("skyblue",0.5))
abline(0,1)

plot(res_AoBoIn$nes~res_AoBoIn$mod,
     xlim = c(0,20), 
     ylim = c(0,4),
     main = "res_AoBoIn",
     ylab = "∆ Nestedness", xlab = "∆ Modularity",
     pch = 16, cex = log(res$RQLa), col = scales::alpha("firebrick",0.5))
points(res_AoBoIn$nes~res_AoBoIn$mod,
       pch = 16, cex = log(res$RQLp), col = scales::alpha("skyblue",0.5))
abline(0,1)


plot(res_AnBoIn$nes~res_AnBoIn$mod,
     xlim = c(0,20), 
     ylim = c(0,4),
     main = "res_AnBoIn",
     ylab = "∆ Nestedness", xlab = "∆ Modularity",
     pch = 16, cex = log(res_AnBoIn$RQLa), col = scales::alpha("firebrick",0.5))
points(res_AnBoIn$nes~res_AnBoIn$mod,
       pch = 16, cex = log(res_AnBoIn$RQLp), col = scales::alpha("skyblue",0.5))
abline(0,1)

plot(res_AoBnIn$nes~res_AoBnIn$mod,
     xlim = c(0,20), 
     ylim = c(0,4),
     main = "res_AnBoIn",
     ylab = "∆ Nestedness", xlab = "∆ Modularity",
     pch = 16, cex = log(res_AoBnIn$RQLa), col = scales::alpha("firebrick",0.5))
points(res_AoBnIn$nes~res_AoBnIn$mod,
       pch = 16, cex = log(res_AoBnIn$RQLp), col = scales::alpha("skyblue",0.5))
abline(0,1)















