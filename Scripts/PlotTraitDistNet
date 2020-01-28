##################
#R script to create the parts for a symmetry network figure 
# Gabriel Muñoz ---> Emma L. Marjakangas
# Jan 2020
##################

### Uncomment this first if need to upload  the rds files (when working outside the cluster server)
#simul <- readRDS("EF_NL_simul_fixAug.rds")
#traitNet <- readRDS("traitNet.rds")
## key number for now
Poss
## Building figure parts as requested by Emma 

#### Four networks with indicated 
#   combinations of strengths of environmental filtering 
#   at the two trophic levels


################ FUNCTIONS ################ 
# function to set color palet
col01 <- function(dist){ 
  dt <- data.frame("col" = rev(RColorBrewer::brewer.pal(10, "Spectral")),"lim" = seq(0.1,1,0.1))
  s1 <- dt$col[1]
  colPal <-sapply(dist, function(x) dt$col[match(round(x, 1), dt$lim)])
  colPal <- ifelse(is.na(colPal), as.character(s1),as.character(colPal))
  
  return(colPal)
  
}
# function to set a pallet from a continuous set of data 
f <- function(x,n=10, pal, rev = F){
  if(rev == F){ 
    rev(RColorBrewer::brewer.pal(n, pal))[cut(x,n)]
  }else{
    (RColorBrewer::brewer.pal(n, pal))[cut(x,n)]
  }
}
# function to calculate trait degree and plot nice bipartite networks
# change plot = T to plot the network as a bipartite (choose a low number of max Iinteractions)
plotInt <- function(net, maxInt, plot = F){ 
  
  net <- net[sample(1:length(net$Var1),
                    maxInt, prob = net$intProb/max(net$intProb)),]
  bnet <- xtabs(intProb ~ Var1 + Var2, net)
  bnet[bnet > 1] <-1
  
  if(plot == T){
    bipartite::plotweb(bnet,labsize = 0.0001,
                       col.high = f(net$trait2, 3,"RdYlGn"),
                       bor.col.high = "white",
                       bor.col.low  = "white",
                       bor.col.interaction = scales::alpha("white",0),
                       col.low =  f(net$trait1, 3,"RdYlGn"),
                       col.interaction = scales::alpha(col01(net$intProb),
                                                       net$intProb/max(net$intProb)),
                       method = "cca")
    
  }

  

  
  return(net)
  
  }
# function to plot the results of plotInt function above. 
makeTraitDens <- function(d){
  plot(density(d$trait1),
       frame = F,
       col = "blue",
       main = "", 
       ylab = "Normalized interaction degree",
       xlab  = "Trait value",
       yaxt = "n",
       xlim = c(0,1),
       ylim = c(0,10))
  points(density(d$trait2), type = "l", col = "red")
}
################# END FUNCTIONS ############

### Plot networks for Morphological match ## big thing this are only considering networks assembled at the 5th quantile

netSS <- traitNet$SEF_SEF_MM_0.2_0.2_0.1_0.1__[,1]$com
netSW <- traitNet$SEF_SEF_MM_0.2_0.2_0.1_1__[,1]$com
netWS <- traitNet$SEF_SEF_MM_0.2_0.2_1_0.1__[,1]$com
netWW <- traitNet$SEF_SEF_MM_0.2_0.2_1_1__[,1]$com

par(mfrow = c(4,2), oma = c(0,0,3,0))
d <- plotInt(netWW, 20, plot = T)
makeTraitDens(d)
title("WeakA-WeakB")
d <- plotInt(netWS, 20, plot = T)
makeTraitDens(d)
title("WeakA-StrongB")
d <- plotInt(netSS, 20, plot = T)
makeTraitDens(d)
title("StrongA-StrongB")
d <- plotInt(netSW, 50, plot = T)
makeTraitDens(d)
title("StrongA-WeakB")
title("Environmental Filtering + Morphological matching, 50 most probable interactions", outer = T)


### Plot networks for Neutral Interactions

netSSnl <- traitNet$SEF_SEF_NL_0.2_0.2_0.1_0.1__[,1]$com
netSWnl <- traitNet$SEF_SEF_NL_0.2_0.2_0.1_1__[,1]$com
netWSnl <- traitNet$SEF_SEF_NL_0.2_0.2_1_0.1__[,1]$com
netWWnl <- traitNet$SEF_SEF_NL_0.2_0.2_1_1__[,1]$com

par(mfrow = c(4,2), oma = c(0,0,3,0))
d <- plotInt(netWWnl, 50, plot = T)
makeTraitDens(d)
title("WeakA-WeakB")
d <- plotInt(netWSnl, 50, plot = T)
makeTraitDens(d)
title("WeakA-StrongB")

d <- plotInt(netSSnl, 50, plot = T)
makeTraitDens(d)
title("StrongA-StrongB")
d <- plotInt(netSWnl, 50, plot = T)
makeTraitDens(d)
title("StrongA-WeakB")
title("Environmental Filtering + Neutral interactions, 50 most probable interactions", outer = T)



## key number for now
Poss
########

rr <- simul[[63]][,1]$metrics[[1]]$bnet
rr1 <- simul[[63]][,1]$metrics[[3]]$bnet
rr2 <- simul[[63]][,1]$metrics[[5]]$bnet
rrnl <- simul[[78]][,1]$metrics[[1]]$bnet

nr <- simul[[61]][,1]$metrics[[1]]$bnet
nr1 <- simul[[61]][,1]$metrics[[3]]$bnet
nr2 <- simul[[61]][,1]$metrics[[5]]$bnet
nrnl <- simul[[76]][,1]$metrics[[5]]$bnet

fr <- simul[[62]][,1]$metrics[[1]]$bnet
fr1 <- simul[[62]][,1]$metrics[[3]]$bnet
fr2 <- simul[[62]][,1]$metrics[[5]]$bnet
frnl <- simul[[77]][,1]$metrics[[5]]$bnet




### Make figure where we show the spatiotemporal effect of a cutoff in interactions 
par(mfrow = c(3,5), oma = c(0,0,0,0))
makeTraitDens(rr)
legend("topright", "Morphological matching", bty = "n", cex = 0.5)
title("1Qt")
abline(v = c(Poss$TraitA[63], Poss$TraitB[63]), col = "grey")
makeTraitDens(rr1)
legend("topright", "Morphological matching", bty = "n", cex = 0.5)
title("3Qt")
abline(v = c(Poss$TraitA[63], Poss$TraitB[63]), col = "grey")
makeTraitDens(rr2)
legend("topright", "Morphological matching", bty = "n", cex = 0.5)
title("5Qt")
abline(v = c(Poss$TraitA[63], Poss$TraitB[63]), col = "grey")
d <- plotInt(simul[[63]][,1]$metrics[[5]]$bnet, 100)
makeTraitDens(d)
legend("topright", "Morphological matching", bty = "n", cex = 0.5)
title("50 most probable")
abline(v = c(Poss$TraitA[63], Poss$TraitB[63]), col = "grey")
makeTraitDens(rrnl)
title("Neutral CA")
legend("topright", "Morphological matching", bty = "n", cex = 0.5)
abline(v = c(Poss$TraitA[78], Poss$TraitB[78]), col = "grey")



makeTraitDens(nr)
legend("top", "Neutral I.A.", bty = "n", cex = 0.5)
makeTraitDens(nr1)
legend("top", "Neutral I.A.", bty = "n", cex = 0.5)
makeTraitDens(nr2)
legend("top", "Neutral I.A.", bty = "n", cex = 0.5)
d <- plotInt(simul[[61]][,1]$metrics[[5]]$bnet, 100)
makeTraitDens(d)
legend("top", "Neutral I.A.", bty = "n", cex = 0.5)
makeTraitDens(nrnl)
legend("top", "Neutral I.A.", bty = "n", cex = 0.5)

makeTraitDens(fr)
legend("top", "Forbidden Links", bty = "n", cex = 0.5)
makeTraitDens(fr1)
legend("top", "Forbidden Links", bty = "n", cex = 0.5)
makeTraitDens(fr2)
legend("top", "Forbidden Links", bty = "n", cex = 0.5)
d <- plotInt(simul[[62]][,1]$metrics[[5]]$bnet, 100)
makeTraitDens(d)
legend("top", "Forbidden Links", bty = "n", cex = 0.5)
makeTraitDens(frnl)
legend("top", "Forbidden Links", bty = "n", cex = 0.5)





