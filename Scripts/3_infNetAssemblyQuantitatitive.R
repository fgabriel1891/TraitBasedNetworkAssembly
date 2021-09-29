###################
## Step 3. The influence of niche-based assembly processes on network structure (weigthed networks)
########################
## Start parallelization 
#####
library(snowfall)

sfInit(parallel=TRUE, cpus=10, type="SOCK")  # start cluster 
sfSource("Scripts/0_Functions.R") # Source functions to the cluster
# Export the general pool to sample species
sfExport("pool")
# Export the probabilistic matrix
sfExport("ProbPool")
## Run null networks in parallel
#####
####################
#CnRoIn
####################
q_CnRoIn <- sfLapply(1:length(unique(pool$site)), function(x) 
  replicate(10,NullNetMod(pool,
                          ProbPool, 
                          repPool = 1000,
                          NPlant =  F, 
                          NAnimal =  T, 
                          nulltype = "quasiswap_count",
                          binary = F,
                          unique(pool$site)[x], 100)))


names(q_CnRoIn) <- unique(pool$site)
####################
# CoRoIn
####################
q_CoRoIn <- sfLapply(1:length(unique(pool$site)), function(x) 
  replicate(10,NullNetMod(pool,
                          ProbPool, 
                          repPool = 1000,
                          NPlant =  F, 
                          NAnimal =  F, 
                          nulltype = "quasiswap_count",
                          binary = F,
                          unique(pool$site)[x], 100)))


names(q_CoRoIn) <- unique(pool$site)


####################
####################
#CnRnIn
####################
q_CnRnIn <- sfLapply(1:length(unique(pool$site)), function(x) 
  replicate(10,NullNetMod(pool,
                          ProbPool, 
                          repPool = 1000,
                          NPlant =  T, 
                          NAnimal =  T, 
                          nulltype = "quasiswap_count",
                          binary = F,
                          unique(pool$site)[x], 100)))
names(q_CnRnIn) <- unique(pool$site)


####################

####################
#CoRnIn
####################
q_CoRnIn <- sfLapply(1:length(unique(pool$site)), function(x) 
  replicate(10,NullNetMod(pool,
                          ProbPool, 
                          repPool = 1000,
                          NPlant =  T, 
                          NAnimal =  F, 
                          nulltype = "quasiswap_count",
                          binary = F,
                          unique(pool$site)[x], 100)))





# CoRnIn <- CoRnIn[c(1,3,5,6,8,10)] 
names(q_CoRnIn) <- unique(pool$site)


############
# Get the mean of replicates
######

q_CnRnIn_x <- sapply(1:6, function(x) apply(q_CnRnIn[[x]], 1, mean))
colnames(q_CnRnIn_x)<- names(q_CnRnIn)

q_CoRoIn_x <- sapply(1:6, function(x) apply(q_CoRoIn[[x]], 1, mean))
colnames(q_CoRoIn_x)<- names(q_CoRoIn)

q_CnRoIn_x <- sapply(1:6, function(x) apply(q_CnRoIn[[x]], 1, mean))
colnames(q_CnRoIn_x)<- names(q_CnRoIn)

q_CoRnIn_x <- sapply(1:6, function(x) apply(q_CoRnIn[[x]], 1, mean))
colnames(q_CoRnIn_x)<- names(q_CoRnIn)


############
# Get the sd of replicates
######

q_CnRnIn_sd <- sapply(1:6, function(x) apply(q_CnRnIn[[x]], 1, sd))
colnames(CnRnIn_x)<- names(q_CnRnIn)

q_CoRoIn_sd <- sapply(1:6, function(x) apply(q_CoRoIn[[x]], 1, sd))
colnames(q_CoRoIn_x)<- names(q_CoRoIn)

q_CnRoIn_sd <- sapply(1:6, function(x) apply(q_CnRoIn[[x]], 1, sd))
colnames(CnRoIn_x)<- names(q_CnRoIn)

q_CoRnIn_sd <- sapply(1:6, function(x) apply(q_CoRnIn[[x]], 1, sd))
colnames(q_CoRnIn_x)<- names(q_CoRnIn)




#####
# Calculate SESnet 
#####

# calculate observed modularity and nestedness 
netSites <- xtabs(frequency~plantCode + animalCode + site,pool)
nullDist <- apply(netSites, 3, NullModSen,100, "DD")
# make binary to calculate observed modularity
netSites <- apply(netSites, 3, makeBinar)
obsModular <- sapply(1:6, function(x) bipartite::computeModules(netSites[[x]])@likelihood)
obsNested <- sapply(1:6, function(x) vegan::nestednodf(netSites[[x]])$statistic[["NODF"]])
Zscores <- sapply(1:6, function(x) (obsModular[x] -nullDist[[x]][[1]]["NulMod.mean"]) / nullDist[[x]][[1]]["NulMod.sd"])
Zscores_nes <- sapply(1:6, function(x) (obsNested[x] -nullDist[[x]][[1]]["NulNest.mean"]) / nullDist[[x]][[1]]["NulNest.sd"])



# Calculate effect sizes from null networks 

# modularity 
q_CnRnIn_SES_mod <- (Zscores-q_CnRnIn_x[1,]) / q_CnRnIn_sd[1,]
q_CoRoIn_SES_mod <- (Zscores-q_CoRoIn_x[1,]) / q_CoRoIn_sd[1,]
q_CnRoIn_SES_mod <- (Zscores-q_CnRoIn_x[1,]) / q_CnRoIn_sd[1,]
q_CoRnIn_SES_mod <- (Zscores-q_CoRnIn_x[1,]) / q_CoRnIn_sd[1,]


## Partialying out separate effects of consumer, resources, and interactions 

# modularity 

q_INTef_SES_mod <- q_CnRnIn_SES_mod + q_CoRnIn_SES_mod - q_CnRoIn_SES_mod - q_CoRoIn_SES_mod
q_INT_SES_mod <- q_CoRoIn_SES_mod
q_R_SES_mod <- q_CoRnIn_SES_mod- q_CoRoIn_SES_mod
q_C_SES_mod <- q_CnRoIn_SES_mod- q_CoRoIn_SES_mod


## Create a new dataframe with ef sizes

# modularity 
q_SESnet_mod <- data.frame(
  q_INTef_SES_mod, 
  q_INT_SES_mod,
  q_R_SES_mod,
  q_C_SES_mod)

#SESnet_mod <- SESnet_mod/max(abs(SESnet_mod))
colnames(q_SESnet_mod) <- c("InterEf", "Interactions", "Resources", "Consumers")


# find partial sum of effects (no accounting for the interaction between effects)

# modularity 
sumAbsEf <- abs(q_SESnet_mod$Resources)+abs(q_SESnet_mod$Consumers)+abs(q_SESnet_mod$Interactions)


# Normalize data by the partial sum to find relative magnitudes

# modularity
q_SESnet_mod_norm <- q_SESnet_mod[,c(2,3,4)]/(sumAbsEf)

q_SESnet_mod_norm$elev <- dataSet1$N$elevation[match(rownames(q_SESnet_mod_norm),
                                                   dataSet1$N$site)]



q_SESnet_mod <- q_SESnet_mod[match(rownames(PSdir),rownames(q_SESnet_mod)),]



q_SESnet_mod$elev <- dataSet1$N$elevation[match(rownames(q_SESnet_mod),
                                              dataSet1$N$site)]





round(q_SESnet_mod, 2)




png("AssemblyMode.png", width = 1000, height = 1000, pointsize = 20, res = 110)
par(mfrow = c(1,1), oma = c(2,2,2,2), mar = c(5,5,2,2))

assMode <- log(abs(q_SESnet_mod$Consumers)/abs(q_SESnet_mod$Resources))


plot(elev_vec2[order]~assMode[order],
     pch = 21, 
     frame = F,
     xlab = "Network assembly mode",
     ylab = "Elevation",
     cex = 1.5,
     cex.lab = 1.5,
     ylim = c(0,4000),
     bg = "#9C413D",
     xlim = c(-2,2))
segments(0,elev_vec2[order],assMode[order],
         elev_vec2[order], 
         col = "#9C413D",
         lwd = 2)

abline(h = c(500,1500,2500,3500), v = 0,lty = 2)

mtext("Bottom-up",3, outer = F, adj = 0 , cex = 1)
mtext("Top-down",3, outer = F, adj = 1 , cex = 1)

legend("bottomright", 
       pch=c(16),
       col = "#9C413D",
       bty = "n",
       legend = c("wModularity"))
dev.off()



