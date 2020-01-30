##################
#R script to create the parts for a symmetry network figure 
# Gabriel Muñoz ---> Emma L. Marjakangas
# Jan 2020
##################

### Uncomment this first if need to upload  the rds files (when working outside the cluster server)
#simul <- readRDS("EF_NL_simul_fixAug.rds")
#traitNet <- readRDS("traitNet.rds")

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
  net$trait2 <- round(net$trait2,3)
  net$trait1 <- round(net$trait1,3)
  bnet <- xtabs(intProb ~ trait1 + trait2, net)
  bnet[bnet > 1] <-1
  
  # 
  # tb <- na.omit(net[match(colnames(bnet),net$trait2),])
  # tb <- tb[order(tb$trait1, decreasing = T),]
  # seq <-   list("seq.low" = tb$trait1, 
  #               "seq.high" = tb$trait2[order(tb$trait2, decreasing = T)])
  # t1 <- net$trait1[match(seq$seq.low,net$trait1)]
  # t2 <- net$trait2[match(seq$seq.high,net$trait2)]
  # 
  # bnet <- bnet[match(seq$seq.low,rownames(bnet))
  #              ,match(seq$seq.high,colnames(bnet))]
  # 
  as.numeric(rownames(bnet))
  if(plot == T){
    bipartite::plotweb(bnet,labsize = 0.0001,
                       col.high = col01( as.numeric(colnames(bnet))),
                       bor.col.high = scales::alpha("grey", 0.3),
                       bor.col.low  = scales::alpha("grey", 0.3),
                       bor.col.interaction = scales::alpha("grey",0.2),
                       col.low =  col01( as.numeric(rownames(bnet)))
                       col.interaction = scales::alpha("grey", 0.3),
                       # col.interaction = scales::alpha(col01(net$intProb),
                       #                                 net$intProb/max(net$intProb)),
                       method = "normal")
    
  }
  return(net)
  
}
# function to plot the results of plotInt function above. 
makeTraitDens <- function(trait1, trait2, add = F, col = col, twolev = F, col2 = col2){
  if(add == T){ 
    points(density(trait1), 
           type = "l",lwd = 2,
         col = col)
  }else{
    plot(density(trait1),
         frame = F,
         col = col,
         lwd = 2,
         main = "", 
         ylab = "Normalized interaction degree",
         xlab  = "Trait value",
         yaxt = "n",
         xlim = c(0,1),
         ylim = c(0,10))
    
    
  }
  if(twolev == T){ 
  points(density(trait2), type = "l", lwd = 2,col = col2)
    
    }

  
  
}
################# END FUNCTIONS ############

# define objects 
### Plot networks for Morphological match ## big thing this are only considering networks assembled at the 5th quantile

netSS <- traitNet$SEF_SEF_MM_0.2_0.2_0.1_0.1__[,1]$com
netSW <- traitNet$SEF_SEF_MM_0.2_0.2_0.1_1__[,1]$com
netWS <- traitNet$SEF_SEF_MM_0.2_0.2_1_0.1__[,1]$com
netWW <- traitNet$SEF_SEF_MM_0.2_0.2_1_1__[,1]$com

### Plot networks for Neutral Interactions

netSSnl <- traitNet$SEF_SEF_NL_0.2_0.2_0.1_0.1__[,1]$com
netSWnl <- traitNet$SEF_SEF_NL_0.2_0.2_0.1_1__[,1]$com
netWSnl <- traitNet$SEF_SEF_NL_0.2_0.2_1_0.1__[,1]$com
netWWnl <- traitNet$SEF_SEF_NL_0.2_0.2_1_1__[,1]$com



### Plot networks 

par(mfrow = c(1,4), mar = c(0,0,0,0))
d1 <- plotInt(netWW, 40, plot = T)
title("Weak-Weak")
d2 <- plotInt(netWS, 40, plot = T)
title("Weak-Strong")
d3 <- plotInt(netSW, 40, plot = T)
title("Strong-Weak")
d4 <- plotInt(netSS, 40, plot = T)
title("Strong-Strong")
title("Morphological matching", outer = T)


d5 <- plotInt(netWWnl, 40, plot = T)
title("Weak-Weak")
d6 <- plotInt(netWSnl, 40, plot = T)
title("Weak-Strong")
d7 <- plotInt(netSWnl, 40, plot = T)
title("Strong-Weak")
d8 <- plotInt(netSSnl, 40, plot = T)
title("Strong-Strong")
title("Neutral interactions", outer = T)

############################################################
# Plot trait distributions (if needed)
##############################

makeTraitDens(d1$trait1, col = "red")


makeTraitDens(d2)
makeTraitDens(d3)
makeTraitDens(d4)


makeTraitDens(d5)
makeTraitDens(d6)
makeTraitDens(d7)
makeTraitDens(d8)

############################## 



makeTraitDens(d)
d <- plotInt(simul[[63]][,1]$metrics[[5]]$bnet, 50)
makeTraitDens(d$trait1, col = scales::alpha(dd[1], 0.5))


se <- round(seq(100, 5000, length.out=10 ))
se <- rev(se)

par(mfrow = c(1,2), mar = c(5,3,3,3))
d <- plotInt(simul[[63]][,1]$metrics[[5]]$bnet, 50)
makeTraitDens(d$trait1, col = scales::alpha(dd[1], 0.5))
dd <- f(seq(1:10), 10,"Spectral", rev = T)
for(i in 1:length(se)){ 
  d <- plotInt(simul[[63]][,1]$metrics[[5]]$bnet, se[i])
  makeTraitDens(d$trait1, col = scales::alpha(dd[i], 0.2), add = T)
  }
       
legend("topright", 
       lty = 1,
       cex = 0.7,
       bty = "n",
       title = "Total interaction richness",
       legend = se,
       col = dd)
title("Consumer trophic level")
d <- plotInt(simul[[63]][,1]$metrics[[5]]$bnet, 50)
makeTraitDens(d$trait2, col = scales::alpha(dd[1], 0.5))
dd <- f(seq(1:10), 10,"Spectral", rev = T)
for(i in 1:length(se)){ 
  d <- plotInt(simul[[63]][,1]$metrics[[5]]$bnet, se[i])
  makeTraitDens(d$trait2, col = scales::alpha(dd[i], 0.2), add = T)
}

legend("topright", 
       lty = 1,
       cex = 0.7,
       bty = "n",
       title = "Total interaction richness",
       legend = se,
       col = dd)
title("Resource trophic level")
title("Environmental filtering + Morphological matching", outer = T)

############################## 
# EF + N.assembly
############################## 

se <- round(seq(100, 5000, length.out=10 ))
se <- rev(se)

par(mfrow = c(1,2), mar = c(5,3,3,3))
d <- plotInt(simul[[61]][,1]$metrics[[5]]$bnet, 50)
makeTraitDens(d$trait1, col = scales::alpha(dd[1], 0.5))
dd <- f(seq(1:10), 10,"Spectral", rev = T)
for(i in 1:length(se)){ 
  d <- plotInt(simul[[61]][,1]$metrics[[5]]$bnet, se[i])
  makeTraitDens(d$trait1, col = scales::alpha(dd[i], 0.2), add = T)
}

legend("topright", 
       lty = 1,
       cex = 0.7,
       bty = "n",
       title = "Total interaction richness",
       legend = se,
       col = dd)
title("Consumer trophic level")



d <- plotInt(simul[[61]][,1]$metrics[[5]]$bnet, 50)
makeTraitDens(d$trait2, col = scales::alpha(dd[1], 0.5))
dd <- f(seq(1:10), 10,"Spectral", rev = T)
for(i in 1:length(se)){ 
  d <- plotInt(simul[[61]][,1]$metrics[[5]]$bnet, se[i])
  makeTraitDens(d$trait2, col = scales::alpha(dd[i], 0.2), add = T)
}

legend("topright", 
       lty = 1,
       cex = 0.7,
       bty = "n",
       title = "Total interaction richness",
       legend = se,
       col = dd)
title("Resource trophic level")
title("Environmental filtering + Neutral assembly", outer = T)

############################## 
# Neutral + Mor.mat
############################## 

se <- round(seq(100, 5000, length.out=10 ))
se <- rev(se)

par(mfrow = c(1,2), mar = c(5,3,3,3))
d <- plotInt(simul[[78]][,1]$metrics[[5]]$bnet, 50)
makeTraitDens(d$trait1, col = scales::alpha(dd[1], 0.5))
dd <- f(seq(1:10), 10,"Spectral", rev = T)
for(i in 1:length(se)){ 
  d <- plotInt(simul[[78]][,1]$metrics[[5]]$bnet, se[i])
  makeTraitDens(d$trait1, col = scales::alpha(dd[i], 0.2), add = T)
}

legend("topright", 
       lty = 1,
       cex = 0.7,
       bty = "n",
       title = "Total interaction richness",
       legend = se,
       col = dd)
title("Consumer trophic level")



d <- plotInt(simul[[78]][,1]$metrics[[5]]$bnet, 50)
makeTraitDens(d$trait2, col = scales::alpha(dd[1], 0.5))
dd <- f(seq(1:10), 10,"Spectral", rev = T)
for(i in 1:length(se)){ 
  d <- plotInt(simul[[78]][,1]$metrics[[5]]$bnet, se[i])
  makeTraitDens(d$trait2, col = scales::alpha(dd[i], 0.2), add = T)
}

legend("topright", 
       lty = 1,
       cex = 0.7,
       bty = "n",
       title = "Total interaction richness",
       legend = se,
       col = dd)
title("Resource trophic level")
title("Neutral assembly + Morphological matching", outer = T)

############################## 






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





