###################################
## Code to replicate the analyses and figures shown in: Trait-based inference of ecological network assembly
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

## Load custom functions (important!)
source("scripts/CustomFunctions.R") # path might change 

###############
## Section 4: SIMULATING NETWORK ASSEMBLY
###############
## dependencies 
#####
library(foreach) # apply functions sequencially 
library(MASS) # stats
library(vegan) # eco stats 
library(yhat) # variance partitioning 

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
ss = CreateSpPool(c(5000, 5000), c(500,500), "uniform")

## hard save the species pool object (for later reproducibility )
saveRDS(ss, "pool.rds")

# Important note if running in a cluster:
### Parallelization starts from here. Thus, set the numbers of cpus accordingly before running this section. 
## Start parallelization (Make sure the snowfall library has been loaded)
library(snowfall)
sfInit(parallel=TRUE, cpus=80, type="SOCK")  # start cluster 
sfSource("scripts/CustomFunctions.R") # Source functions to the cluster
# Export the scenario object(s) to the cluster (i.e. objects created in step 1)
sfExport("ss")
# Fixed trait optima
sfExport("Poss")
sfExport("PossNETall")

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
saveRDS(simulx, "EF_NL_simul_fixAug.rds")
saveRDS(simul_0, "BUTD_ef.rds.")
# (End of simulations)

#### End cluster
snowfall::sfStop()
################


############################################################
## Step 3) Interpret and the raw simulation output and isolate the niche based effects.
############################################################


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
RES_DD$PSSdir <- log((RES_DD$sigmaA/RES_DD$sigmaB)^2)
RES_DD$PSSmag <- log(sqrt((RES_DD$sigmaA-RES_DD$sigmaB)^2)+1)
RES_DD$PSSdir <- RES_DD$PSSdir /max(RES_DD$PSSdir, na.rm = T)
RES_DD$PSSmag <- RES_DD$PSSmag /max(RES_DD$PSSmag, na.rm = T)

RES_DD_2$PSSdir <- log((RES_DD_2$sigmaA/RES_DD_2$sigmaB)^2)
RES_DD_2$PSSmag <- log(sqrt((RES_DD_2$sigmaA-RES_DD_2$sigmaB)^2)+1)
RES_DD_2$PSSdir <- RES_DD_2$PSSdir /max(RES_DD_2$PSSdir, na.rm = T)
RES_DD_2$PSSmag <- RES_DD_2$PSSmag /max(RES_DD_2$PSSmag, na.rm = T)

#########


##############
## Figure 3. Effects of community and interaction assembly processes on the structure of simulated resource-consumer networks
##############

########
## Figure 3a: Two dimensional space defined by modularity and nestedness axis.
#####
#######
par(mar=c(6,6,4,4), las = 1)
plot(NODFx~QZscorex,
     frame = F,
     xlim = c(-50, 200),
     ylim = c(-50,200),
     pch = ifelse( RES_DD$intHyp == "NL", 16, 
                   ifelse( RES_DD$intHyp == "FL", 15, 17)),
     col = "grey90",
     ylab = "", 
     xlab = "", 
     cex.axis = 1,
     cex = 1, 
     data = RES_DD)
Plot_ConvexHull( RES_DD$QZscorex[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL")],
                 RES_DD$NODFx[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL")], 
                 scales::alpha("gray",0.5))

Plot_ConvexHull( RES_DD$QZscorex[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM")],
                 RES_DD$NODFx[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM")], 
                 scales::alpha("gray",1))

Plot_ConvexHull( RES_DD$QZscorex[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL")],
                 RES_DD$NODFx[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL")], 
                 scales::alpha("gray",0.5))
points(NODFx~QZscorex,data = RES_DD, col = "grey60", cex=  1,
       pch = ifelse( RES_DD$intHyp == "NL", 1, 
                     ifelse( RES_DD$intHyp == "FL", 0, 2)))


#### 
# neutral arrows

segments("x0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$QZscorex
         -RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$QZscoresd,
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$QZscorex
         +RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$QZscoresd,
         "y0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$NODFx, 
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$NODFx,
         lwd = 2, lty = 1, col = scales::alpha("black",0.9))

segments("y0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$NODFx
         -RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("NL"),]$NODFsd,
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$NODFx
         +RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("NL"),]$NODFsd,
         "x0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$QZscorex, 
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$QZscorex, 
         lwd = 2, lty = 1,col = scales::alpha("black",0.9))

segments("x0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$QZscorex
         -RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$QZscoresd,
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$QZscorex
         +RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$QZscoresd,
         "y0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$NODFx, 
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$NODFx,
         lwd = 2, lty = 1, col = scales::alpha("black",0.9))

segments("y0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$NODFx
         -RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("MM"),]$NODFsd,
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$NODFx
         +RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("MM"),]$NODFsd,
         "x0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$QZscorex, 
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$QZscorex, 
         lwd = 2, lty = 1,col = scales::alpha("black",0.9))

segments("x0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$QZscorex
         -RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$QZscoresd,
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$QZscorex
         +RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$QZscoresd,
         "y0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$NODFx, 
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$NODFx,
         lwd = 2, lty = 1, col = scales::alpha("black",0.9))

segments("y0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$NODFx
         -RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("FL"),]$NODFsd,
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$NODFx
         +RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("FL"),]$NODFsd,
         "x0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$QZscorex, 
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$QZscorex, 
         lwd = 2, lty = 1,col = scales::alpha("black",0.9))
abline(h=0, lty = 2)
abline(v=0, lty = 2)

#######
# Figure 3b: Bottom-up and Top-down effects of resource and community assembly processes. 
########

## define a function to calculate effect sizes 
con_res_efz <- function(RES_butd, RES1, filt, int, ab = c("a","b")){
  
  if(ab == "a"){
    
    n_xdif <- RES_butd$NODFx[RES_butd$sigmaA == filt & RES_butd$intHyp == int]-RES1$NODFx[RES1$sigmaA == filt & RES1$intHyp == int]
    n_sd_pool <- sqrt(RES_butd$NODFsd[RES_butd$sigmaA == filt & RES_butd$intHyp == int]^2+RES1$NODFsd[RES1$sigmaA == filt & RES1$intHyp == int]^2)
    n_efz <- n_xdif/n_sd_pool
    
    q_xdif <- RES_butd$QZscorex[RES_butd$sigmaA == filt & RES_butd$intHyp == int]-RES1$QZscorex[RES1$sigmaA == filt & RES1$intHyp == int]
    q_sd_pool <- sqrt(RES_butd$QZscoresd[RES_butd$sigmaA == filt & RES_butd$intHyp == int]^2+RES1$QZscoresd[RES1$sigmaA == filt & RES1$intHyp == int]^2)
    q_efz <- q_xdif/q_sd_pool
    pss <- RES1$PSSdir[RES1$sigmaA == filt & RES1$intHyp == int]
    
    psm <- RES1$PSSmag[RES1$sigmaA == filt & RES1$intHyp == int]
    
  }
  
  
  if(ab == "b"){
    n_xdif <- RES_butd$NODFx[RES_butd$sigmaB == filt & RES_butd$intHyp == int]-RES1$NODFx[RES1$sigmaB == filt & RES1$intHyp == int]
    n_sd_pool <- sqrt(RES_butd$NODFsd[RES_butd$sigmaB == filt & RES_butd$intHyp == int]^2+RES1$NODFsd[RES1$sigmaB == filt & RES1$intHyp == int]^2)
    n_efz <- n_xdif/n_sd_pool
    
    
    q_xdif <- RES_butd$QZscorex[RES_butd$sigmaB == filt & RES_butd$intHyp == int]-RES1$QZscorex[RES1$sigmaB == filt & RES1$intHyp == int]
    q_sd_pool<- sqrt(RES_butd$QZscoresd[RES_butd$sigmaB == filt & RES_butd$intHyp == int]^2+RES1$QZscoresd[RES1$sigmaB == filt & RES1$intHyp == int]^2)
    q_efz <- q_xdif/q_sd_pool
    pss <- RES1$PSSdir[RES1$sigmaB == filt & RES1$intHyp == int]
    psm <- RES1$PSSmag[RES1$sigmaB == filt & RES1$intHyp == int]
  }
  
  
  return(data.frame(n_efz,q_efz,pss, psm))
  
  
}



cr_nl_a <- sapply(seq(0.1,1, 0.2), function(x) con_res_efz(RES_DD_2, RES_DD, x, "NL","a"))
cr_fl_a <- sapply(seq(0.1,1, 0.2), function(x) con_res_efz(RES_DD_2, RES_DD, x, "FL", "a"))
cr_mm_a <- sapply(seq(0.1,1, 0.2), function(x) con_res_efz(RES_DD_2, RES_DD, x, "MM","a"))


cr_nl_b <- sapply(seq(0.1,1, 0.2), function(x) con_res_efz(RES_DD_2, RES_DD, x, "NL","b"))
cr_fl_b <- sapply(seq(0.1,1, 0.2), function(x) con_res_efz(RES_DD_2, RES_DD, x, "FL", "b"))
cr_mm_b <- sapply(seq(0.1,1, 0.2), function(x) con_res_efz(RES_DD_2, RES_DD, x, "MM","b"))




### define a function to arrange the data.frame to plot
toplot <- function(cr_nl_a){
  toplot <- data.frame("nodf" = c(cr_nl_a[,1]$n_efz,cr_nl_a[,2]$n_efz,cr_nl_a[,3]$n_efz,cr_nl_a[,4]$n_efz,cr_nl_a[,5]$n_efz),
                       "q"= c(cr_nl_a[,1]$q_efz,cr_nl_a[,2]$q_efz,cr_nl_a[,3]$q_efz,cr_nl_a[,4]$q_efz,cr_nl_a[,5]$q_efz),
                       "size1" = c(1-sapply(seq(0.1,1,0.2),function(x)rep(x,5))),
                       "size2" = rep(1-seq(0.1,1,0.2),5),
                       "pss" = c(cr_nl_a[,1]$pss,cr_nl_a[,2]$pss,cr_nl_a[,3]$pss,cr_nl_a[,4]$pss,cr_nl_a[,5]$pss),
                       "psm" = c(cr_nl_a[,1]$psm,cr_nl_a[,2]$psm,cr_nl_a[,3]$psm,cr_nl_a[,4]$psm,cr_nl_a[,5]$psm)
  )
  return(toplot)
  
}


plot_nl_a <- toplot(cr_nl_a)
plot_fl_a <- toplot(cr_fl_a)
plot_mm_a <- toplot(cr_mm_a)

plot_nl_b <- toplot(cr_nl_b)
plot_fl_b <- toplot(cr_fl_b)
plot_mm_b <- toplot(cr_mm_b)


newplot_fl <- rbind(plot_fl_b,plot_fl_a)
newplot_mm <- rbind(plot_mm_b,plot_mm_a)
newplot_nl <- rbind(plot_nl_b,plot_nl_a)



addBin <- function(newplot){
  a <- newplot[!newplot$pss == 0,]
  b <- newplot[newplot$pss == 0,]
  
  
  a$bin <- as.numeric(cut(abs(a$pss),5))
  b$bin <- rep(0,10)
  c <- rbind(a,b)
  return(c)
}


newplot_fl <- addBin(newplot_fl)
newplot_mm <- addBin(newplot_mm)
newplot_nl <- addBin(newplot_nl)



par(mfrow = c(1,3), mar = c(4,0,1,0), oma = c(2,2,2,2))
boxplot(newplot_nl$nodf~(newplot_nl$bin), frame = F,
        ylim = c(-10,10), xaxt = "n",
        cex.axis = 2, cex = 1, 
        ylab = c(""), xlab = c(""),
        pch = 16,main = "",
        col = c("white", 
                sapply(seq(0.2,1,0.2), function(x) scales::alpha("gray20", x))))
abline( h = 0, lwd = 2)
boxplot(newplot_fl$nodf~(newplot_fl$bin), frame = F,
        ylim = c(-10,10), xaxt = "n",yaxt = "n",
        cex.axis = 2, cex = 1, 
        ylab = c(""), xlab = c(""),
        pch = 16,main = "",
        col = c("white", 
                sapply(seq(0.2,1,0.2), function(x) scales::alpha("gray20", x))))
abline( h = 0, lwd  = 2)
boxplot(newplot_mm$nodf~(newplot_mm$bin),frame = F,
        ylim = c(-10,10), xaxt = "n",yaxt = "n",
        cex.axis = 2, cex = 1, 
        ylab = c(""), xlab = c(""),
        pch = 16,main = "",
        col = c("white", 
                sapply(seq(0.2,1,0.2), function(x) scales::alpha("gray20", x))))
abline( h = 0, lwd  = 2)








par(mfrow = c(1,3), mar = c(4,0,1,0), oma = c(2,2,2,2))
boxplot(newplot_nl$q~(newplot_nl$bin), frame = F,
        ylim = c(-10,10), xaxt = "n",yaxt = "n",
        cex.axis = 2, cex = 1, 
        ylab = c(""), xlab = c(""),
        pch = 16,main = "",
        col = c("white", 
                sapply(seq(0.2,1,0.2), function(x) scales::alpha("gray20", x))))
abline( h = 0, lwd = 2)
boxplot(newplot_fl$q~(newplot_fl$bin), frame = F,
        ylim = c(-10,10), xaxt = "n",yaxt = "n",
        cex.axis = 2, cex = 1, 
        ylab = c(""), xlab = c(""),
        pch = 16,main = "",
        col = c("white", 
                sapply(seq(0.2,1,0.2), function(x) scales::alpha("gray20", x))))
abline( h = 0, lwd  = 2)
boxplot(newplot_mm$q~(newplot_mm$bin),
        ylim = c(-10,10), xaxt = "n",yaxt = "n",
        cex.axis = 2, cex = 1, frame = F, 
        ylab = c(""), xlab = c(""),
        pch = 16,main = "",
        col = c("white", 
                sapply(seq(0.2,1,0.2), function(x) scales::alpha("gray20", x))))
abline( h = 0, lwd  = 2)
###############
###############

################
################
# Section 5: INFERRING NETWORK ASSEMBLY 
################
################


############
# Step 1. The type and strength of niche-based community and interaction assembly processes
############


######
# a) Process-based species pool delineation
#####

# Slope and elevation are derived from a digital elevation model obtained with radar satellite remote sensing (NASA SRTM shuttle) (https://cgiarcsi.community/data/srtm-90m-digital-elevation-database-v4-1/)
# Temperature and precipitation correspond to the variables BIO01 and BIO12 of WORLDCLIM (https://developers.google.com/earth-engine/datasets/catalog/WORLDCLIM_V1_BIO#bands)


# Load libraries and corresponding data
library(gdistance)
library(ade4)
library(vegan)

# load raw coordinates of collection data 
locCoord <- readxl::read_xlsx("data_package/data/CoordinatesEcuador.xlsx") ## path can change
# get centroids for each of collected sites
centroids <- data.frame("lon" = aggregate(locCoord$Lon, by = list(locCoord$Plot), mean)$x,
                        "lat" = aggregate(locCoord$Lat, by = list(locCoord$Plot), mean)$x,
                        "site" = unique(locCoord$Plot)[order(unique(locCoord$Plot))])

# Load environmental resistance layers (exported from earthengine, https://code.earthengine.google.com/14fe67bd4749c1f90a253f6cfaa05ff9)

# path can change 
SlopeLoja <- raster::raster("data_package/data/LojaSlope.tif")
ClimLoja <- raster::raster("data_package/data/LojaClim.tif")
PrecLoja <- raster::raster("data_package/data/LojaPrec.tif")
ElevLoja <- raster::raster("data_package/data/LojaElevation.tif")

# Visualize sites in environmental space

par(mfrow = c(2,2), mar = c(2,2,2,2))
raster::plot(SlopeLoja, main = "Slope")
points(centroids$lon,centroids$lat, pch = "+")

raster::plot(ClimLoja, main = "Temp")
points(centroids$lon,centroids$lat, pch = "+")

raster::plot(PrecLoja, main = "Prec")
points(centroids$lon,centroids$lat, pch = "+")

raster::plot(ElevLoja, main = "Elev")
points(centroids$lon,centroids$lat, pch = "+")


# Let's apply the custom made function to calculate cost-path distances in environmental space and standardize the distance matrices based on their range maxima 

eDistSlope <- CustomCostDist(SlopeLoja, centroids)
eDistTemp <- CustomCostDist(ClimLoja, centroids)
eDistPrec <- CustomCostDist(PrecLoja, centroids)
eDistElev <- CustomCostDist(ElevLoja, centroids)

# standardize variables

slopeDist <- decostand(eDistSlope, "range")
ClimDist <- decostand(eDistTemp, "range")
eDistPrec <- decostand(eDistPrec, "range")
eDistElev <- decostand(eDistElev, "range")

# Reducing the dimensionality of the environmental distance matrices into linear components

pcaSlope <- princomp(slopeDist)
pcaTemp <- princomp(ClimDist)
pcaPrec <- princomp(eDistPrec)
pcaElev <- princomp(eDistElev)


# How many components are relevant for each variable? 

par(mfrow = c(2,2))
plot(pcaSlope)
plot(pcaTemp)
plot(pcaPrec)
plot(pcaElev)

# It seems that 2 components for each variable are sufficient and explain most of the variance of each environmental distance matrix
# Let's put together  a new data.frame with the relevant components 

EnvVar <- data.frame("Temp" =scores(pcaTemp)[,c(1:2)],
                     "Slope" =scores(pcaSlope)[,c(1:2)],
                     "Prec" =scores(pcaPrec)[,c(1:2)],
                     "Elev" =scores(pcaElev)[,c(1:2)])


# make a new PCA to reduce dimensionality
EnvPCA <- vegan::rda(EnvVar)
dev.off()
biplot(EnvPCA) # first axis = 59% of variation 
# extract scores from the first axis
siteScores <- scores(EnvPCA)$sites[,1]
# calculate pairwise euclidian distances from first axis 
# and normalize to inverse of max distance to create a probability distribution
distAx1 <- dist(siteScores)
probMatrix <- 1-(distAx1/max(distAx1))


# We have defined now our probability matrix for the environmental variables, this represents the probabilities of species contribution 
# from each of the site communities to generate the species pools for any other site. 
# Based only on environmental variables (slope, elevation, precipitation, temperature)

probMatrix # environmental probability matrix 

# Since we have defined our environmental probabilities, let's now define the geographic probabilities. 
# We will do this by calculating simple euclidean distances from one site to another (i.e. as the crow flies distances). 
# We will transform this matrix into a probabilistic one, using the same equation we used above for the environmental matrix 
# $$ 1-\frac{D}{max(D)}$$ being now D, the pairwise euclidian distances in geographical space. 

## construct now a probability matrix based on euclidean distances (as the crow flies)
EucDist <- spatstat::pairdist(centroids$lat, centroids$lon)
# name matrix appropiately
rownames(EucDist) <- colnames(EucDist) <- centroids$site
# standardize 
EucDist <- decostand(EucDist, "range")
# make into probMat
EucDist <- 1-EucDist
EucDist # probability matrix

# Since we have now our environmental and geographical matrices let's combine them into a single one. 
# We will do this by assigning equal weights to each one. Since they are both normalized into a 0-1 range we can use the sum of both matrices multiplied each by the same scalar (0.5)
# $$ PM = 0.5EM + 0.5GM $$ where $PM$ = the probabilistic matrix to sample sites when generating our species pools, $EM$ = Environmental probability matrix, $GM$ Geographical probability matrix. 

## weight the probability matrix based on environmental variables, with the one of simple euclidean distances, assign equal weights
ProbPool <- (0.5*(EucDist))+ (0.5*(as.matrix(probMatrix)))

# Let's visualize how our probabilisitic species pools are delineated among sites. 
# Visualize the delination of species pool selection 
par(oma =c(4,4,2,2))
heatmap(ProbPool)
dev.off()

##########
#####
# b) RQL analyses to get a scaled trait for plants and frugivores
#####
###############################################################################
# load packages and data
################################################################################
library(bipartite)
library(ade4)
library(FD)
library(picante)
load("data_package/data/dataSet1.RData")
################################################################################
# check data
################################################################################
str(dataSet1)

data.frame(dataSet1$N[match(dataSet1$pT$plantCode, dataSet1$N$plantCode[-c(1:4)]),][c("site", "elevation")],
           dataSet1$N[match(dataSet1$pT$plantCode, dataSet1$N$plantCode[-c(1:4)]),][c("site", "elevation")])


dataSet1$pT[c("plantCode", "fruitLength")]
# pT - plant traits
# aT - animal traits
# N - plant-animal interactions
################################################################################
# get some first impression: plot interaction network on each elevation
# and count number of links and species
################################################################################
N <- xtabs(frequency ~ plantCode + animalCode + elevation, data = dataSet1$N)
par(mfrow = c(3, 1))
apply(N, 3, plotweb)
apply(N, 3, FUN = function(x) sum(x > 0)) # number of unique links on each elevation
apply(N, 3, FUN = function(x) sum(rowSums(x) > 0)) # number of plant species on each elevation
apply(N, 3, FUN = function(x) sum(colSums(x) > 0)) # number of animal species on each elevation
################################################################################
# setup data matrices for RLQ analysis based on the pooled binary metaweb
# (as in Albrecht et al. 2018, Nat. Comm. 9:3177)
################################################################################
pTRLQ <- sqrt(dataSet1$pT[, 2:5]) # sqrt-transform traits to improve normality
aTRLQ <- sqrt(dataSet1$aT[, 2:5]) # sqrt-transform traits to improve normality
nRLQ <- xtabs(frequency ~ plantCode + animalCode, data = dataSet1$N)
nRLQ <- nRLQ[rownames(pTRLQ), rownames(aTRLQ)]
nRLQ <- decostand(nRLQ, "pa")
################################################################################
# perform correspondence and principal components analyses
################################################################################
coa <- dudi.coa(as.data.frame(nRLQ), scannf = FALSE, nf = 4)
duP <- dudi.pca(as.data.frame(pTRLQ), scannf = FALSE, center = TRUE, scale = TRUE, nf = 4, row.w = coa$lw)
duA <- dudi.pca(as.data.frame(aTRLQ), scannf = FALSE, center = TRUE, scale = TRUE, nf = 4, row.w = coa$cw)
################################################################################
# perform RLQ analysis
################################################################################
RLQ <- rlq(duP, coa, duA, scannf = FALSE, nf = 4)
RLQ
################################################################################
# perform permutation tests
################################################################################
# Global test for association between trait spaces of plants and animals
################################################################################
trRLQ <- fourthcorner2(as.data.frame(pTRLQ),
                       as.data.frame(nRLQ),
                       as.data.frame(aTRLQ),
                       p.adjust.method.G = "fdr",
                       modeltype = 6, nrepet = 9999)$trRLQ
trRLQ
################################################################################
# Test for axis-specific correlations between plant and animal trait spaces
# NOTE: only tests of 1st vs. 1st and 2nd vs. 2nd axis are meaningful
################################################################################
fcRLQ <- fourthcorner.rlq(RLQ, type = "axes",
                          p.adjust.method.G = "fdr",
                          p.adjust.method.D = "fdr",
                          p.adjust.D = "levels",
                          modeltype = 6, nrepet = 9999)
fcRLQ
################################################################################
# Test for association of animal traits with first and second axis of animal trait space
################################################################################
qAxes <-  fourthcorner.rlq(RLQ, type = "Q.axes",
                           p.adjust.method.G = "fdr",
                           p.adjust.method.D = "fdr",
                           p.adjust.D = "levels",
                           modeltype = 6, nrepet = 9999)
qAxes
################################################################################
# Test for association of plant traits with first and second axis of plant trait space
################################################################################
rAxes <- fourthcorner.rlq(RLQ, type = "R.axes",
                          p.adjust.method.G = "fdr",
                          p.adjust.method.D = "fdr",
                          p.adjust.D = "levels",
                          modeltype = 6, nrepet = 9999)
rAxes
################################################################################
################################################################################
# Figure 1. Plot of RLQ-fourth-corner analysis
################################################################################
#pdf("output/figure_01.pdf", width = 0.394 * 12, height = 0.394 * 6)
m <- matrix(rep(c(1, 2), each = 4), 1, 8, byrow = TRUE)
m <- apply(m, 2, rep, each = 4)
m[16] <- 3
layout(m)
par(cex = 0.7, mar = c(0.5, 0.5, 0.5, 0.5) + 0.1, oma = c(2, 2, 2, 0))
p <- RLQ$mR
a <- RLQ$mQ
N <- nRLQ
p <- cbind(p, (rowSums(N > 0)/max(rowSums(N > 0)) + 0.25) * 2)
a <- cbind(a, (colSums(N > 0)/max(colSums(N > 0)) + 0.25) * 2)
p <- p[order(p[, 3], decreasing = TRUE), ]
a <- a[order(a[, 3], decreasing = TRUE), ]
Ntab <- as.data.frame.table(N > 0)
Ntab$pX <- p[as.character(Ntab$plantCode), 1]
Ntab$pY <- p[as.character(Ntab$plantCode), 2]
Ntab$aX <- a[as.character(Ntab$animalCode), 1]
Ntab$aY <- a[as.character(Ntab$animalCode), 2]
Ntab <- subset(Ntab, Freq > 0)
plot.new()
plot.window(xlim = c(-1, 1) * 1.1, ylim = c(-1, 1) * 1.1, asp = 1)
abline(h = 0, v = 0, col = "gray")
arrows(0, 0, RLQ$c1[, 1], RLQ$c1[, 2], col = "firebrick1",
       code = 2, angle = 30, length = 0.05, lwd = 1, lty = c(1, 1, 5, 2))
arrows(0, 0, RLQ$l1[, 1], RLQ$l1[, 2], col = "deepskyblue1",
       code = 2, angle = 30, length = 0.05, lwd = 1, lty = c(1, 1, 5, 2))
text(RLQ$c1 * 1.1, gsub("_", " ", rownames(RLQ$c1), fixed = TRUE),
     col = "firebrick1", font = 1, cex = 1, xpd = TRUE)
text(RLQ$l1 * 1.1, gsub("_", " ", rownames(RLQ$l1), fixed = TRUE),
     col = "deepskyblue1", font = 1, cex = 1, xpd = TRUE)
mtext(side = 3, line = 0, adj = -0.1, cex = 1, letters[1], font = 2)
legend("topright", bty = "n", lty = c(1, 5, 2), title = "Trait type",
       text.col = c("black", "black", "black"),
       legend = c("Matching", "Energy", "Foraging"))
plot.new()
plot.window(xlim = range(c(p[, 1], a[, 1])), ylim = range(c(p[, 2], a[, 2])), asp = 1)
with(Ntab, arrows(pX, pY, aX, aY, code = 0, lwd = 0.5, col = gray(0.75)))
points(p[, 1:2], pch = 21, cex = p[, 5], lwd = 0.5,
       col = "black", bg = "deepskyblue1")
points(a[, 1:2], pch = 21, cex = a[, 5], lwd = 0.5,
       col = "black", bg = "firebrick1")
mtext(side = 3, line = 0, adj = 0, cex = 1, letters[2], font = 2)
legend("topleft", bty = "n",
       text.col = c("deepskyblue1", "firebrick1", "gray"),
       legend = c("Plants", "Animals", "Interactions"))
par(mar = c(0.5, 0.5, 0.5, 0.5) + 0.1)
barplot(RLQ$eig/sum(RLQ$eig), ylim = c(0, 1), axes = FALSE)
mtext(side = 3, line = 0, adj = 0.5, cex = 0.5, "Eigenvalues")
axis(2, at = seq(0, 1, 1), las = 1, mgp = c(1.5, 0.75, 0))
axis(1, labels = c(1:4),at = seq(1, 4, 1), las = 1,tick = F, mgp = c(1.5, 0.75, 0))
#dev.off()


################################################################################
################################################################################
# the scores from the RLQ analysis can be used to calculate functional diversity
# indices on each elevation and conduct null model analyses
################################################################################
# extract normed species scores of plants and animals for further analysis
# these scores represent the position of plants and animals in the
# respective trait spaces
################################################################################
# Animals - extract normed species scores in trait space & create distance matrix
################################################################################
RLQ$mQ
aDis <- dist(RLQ$mQ)
################################################################################
# Plants - extract normed species scores in trait space & create distance matrix
################################################################################
RLQ$mR
pDis <- dist(RLQ$mR)
################################################################################
## How interactions look into a trait space of RQL 
################################################################################
intRQL <- dplyr::bind_rows(RLQ$mQ,RLQ$mR)
rownames(intRQL) <- c(rownames(RLQ$mQ), rownames(RLQ$mR))
intRQL <- data.frame(dataSet1$N, "RLQan" = RLQ$mQ[match(dataSet1$N$animalCode,rownames(RLQ$mQ)),][,1])
intRQL <- data.frame(intRQL, "RLQpla" = RLQ$mQ[match(dataSet1$N$plantCode,rownames(RLQ$mR)),][,1])
################################################################################
pool <- data.frame(
  dataSet1$N,
  dataSet1$pT[match(dataSet1$N$plantCode, dataSet1$pT$plantCode),][,c(2:5)],
  dataSet1$aT[match(dataSet1$N$animalCode, dataSet1$aT$animalCode),][,c(2:5)],
  "RQLa" = intRQL$RLQan[match(dataSet1$N$animalCode, intRQL$animalCode)],
  "RQLp" = intRQL$RLQpla[match(dataSet1$N$plantCode, intRQL$plantCode)])
pool$intID <- paste0(pool$plantCode,pool$animalCode)
#########
################################################################################
### c) Make sensitivity analyses to different metrics to calculate SES 
######


# Plants
# randomized pool -sd
myEF_r_sd_p <- test4EF(pool = pool,
                       traitName = "RQLp", 
                       nrep = 100, 
                       side = "P",
                       sdOrRange = "sd",
                       RoP = "Random",
                       ProbPool = ProbPool )
myEF2_r_sd_a <- test4EF(pool = pool, 
                        traitName = "RQLa",
                        nrep = 100, 
                        side = "A",
                        sdOrRange = "sd",
                        RoP = "Random", 
                        ProbPool = ProbPool )
# randomized pool -range
myEF_r_ran_p <- test4EF(pool = pool, 
                        traitName = "RQLp", 
                        nrep = 100, 
                        side = "P",
                        sdOrRange = "range",
                        RoP = "Random", 
                        ProbPool = ProbPool )



myEF2_r_ran_a <- test4EF(pool = pool, 
                         traitName = "RQLa",
                         nrep = 100,
                         side = "A",
                         sdOrRange = "range",
                         RoP = "Random", 
                         ProbPool = ProbPool )
# probabilistic pool -sd
myEF_pro_sd_p <- test4EF(pool = pool, 
                         traitName = "RQLp",
                         nrep = 100, 
                         side = "P", 
                         sdOrRange = "sd",
                         RoP = "Prob",
                         ProbPool = ProbPool )
myEF2_pro_sd_a <- test4EF(pool = pool, 
                          traitName = "RQLa", 
                          nrep = 100, 
                          side = "A", 
                          sdOrRange = "sd",
                          RoP = "Prob", 
                          ProbPool = ProbPool )
# probabilistic pool -range
myEF_pro_ran_p <- test4EF(pool = pool, 
                          traitName = "RQLp", 
                          nrep = 100, 
                          side = "P", 
                          sdOrRange = "range",
                          RoP = "Prob", 
                          ProbPool = ProbPool )
myEF2_pro_ran_a <- test4EF(pool = pool, 
                           traitName = "RQLa",
                           nrep = 100, side = "A", 
                           sdOrRange = "range",
                           RoP = "Prob", 
                           ProbPool = ProbPool )
#####

########################
## Step 2. The influence of niche-based assembly processes on network structure
########################
## Start parallelization 
#####
sfInit(parallel=TRUE, cpus=10, type="SOCK")  # start cluster 
sfSource("scripts/CustomFunctions.R") # Source functions to the cluster
# Export the general pool to sample species
sfExport("pool")
# Export the probabilistic matrix
sfExport("ProbPool")
## Run null networks in parallel
#####
####################
#CnRoIn
####################
CnRoIn <- sfLapply(1:length(unique(pool$site)), function(x) 
  replicate(10,NullNetMod(pool,
                          ProbPool, 
                          repPool = 1000,
                          NPlant =  F, 
                          NAnimal =  T, 
                          nulltype = "DD",
                          unique(pool$site)[x], 100)))


names(CnRoIn) <- unique(pool$site)
####################
# CoRoIn
####################
CoRoIn <- sfLapply(1:length(unique(pool$site)), function(x) 
  replicate(10,NullNetMod(pool,
                          ProbPool, 
                          repPool = 1000,
                          NPlant =  F, 
                          NAnimal =  F, 
                          nulltype = "DD",
                          unique(pool$site)[x], 100)))


names(CoRoIn) <- unique(pool$site)


####################
####################
#CnRnIn
####################
CnRnIn <- sfLapply(1:length(unique(pool$site)), function(x) 
  replicate(10,NullNetMod(pool,
                          ProbPool, 
                          repPool = 1000,
                          NPlant =  T, 
                          NAnimal =  T, 
                          nulltype = "DD",
                          unique(pool$site)[x], 100)))
names(CnRnIn) <- unique(pool$site)


####################

####################
#CoRnIn
####################
CoRnIn <- sfLapply(1:length(unique(pool$site)), function(x) 
  replicate(10,NullNetMod(pool,
                          ProbPool, 
                          repPool = 1000,
                          NPlant =  T, 
                          NAnimal =  F, 
                          nulltype = "DD",
                          unique(pool$site)[x], 100)))





# CoRnIn <- CoRnIn[c(1,3,5,6,8,10)] 
names(CoRnIn) <- unique(pool$site)


############
# Get the mean of replicates
######

CnRnIn_x <- sapply(1:6, function(x) apply(CnRnIn[[x]], 1, mean))
colnames(CnRnIn_x)<- names(CnRnIn)

CoRoIn_x <- sapply(1:6, function(x) apply(CoRoIn[[x]], 1, mean))
colnames(CoRoIn_x)<- names(CoRoIn)

CnRoIn_x <- sapply(1:6, function(x) apply(CnRoIn[[x]], 1, mean))
colnames(CnRoIn_x)<- names(CnRoIn)

CoRnIn_x <- sapply(1:6, function(x) apply(CoRnIn[[x]], 1, mean))
colnames(CoRnIn_x)<- names(CoRnIn)

## Calculating confidence intervals 

sapply(1:2, function(x) Rmisc::CI(CnRnIn_x[x,], 0.95 ))

############
# Get the sd of replicates
######

CnRnIn_sd <- sapply(1:6, function(x) apply(CnRnIn[[x]], 1, sd))
colnames(CnRnIn_x)<- names(CnRnIn)

CoRoIn_sd <- sapply(1:6, function(x) apply(CoRoIn[[x]], 1, sd))
colnames(CoRoIn_x)<- names(CoRoIn)

CnRoIn_sd <- sapply(1:6, function(x) apply(CnRoIn[[x]], 1, sd))
colnames(CnRoIn_x)<- names(CnRoIn)

CoRnIn_sd <- sapply(1:6, function(x) apply(CoRnIn[[x]], 1, sd))
colnames(CoRnIn_x)<- names(CoRnIn)



## Process strength symmetry direction 


PSdir <- data.frame("PSdir_random_range" = (log((myEF2_r_ran_a/myEF_r_ran_p)^2)),
                    "PSdir_random_sd" = (log((myEF2_r_sd_a/myEF_r_sd_p)^2)),
                    "PSdir_process_range" = (log((myEF2_pro_ran_a/myEF_pro_sd_p)^2)),
                    "PSdir_process_sd" = (log((myEF2_pro_sd_a/myEF_pro_sd_p)^2)))
PSdir <- PSdir[match(colnames(CoRnIn_x), rownames(PSdir)),]

PSdir <- data.frame(apply(PSdir, 2, function(x) x/max(abs(x))))

## Process strength symmetry magnitude 

PSmag <- data.frame(
  "PSmag_random_range" = decostand(log(sqrt((-myEF2_r_ran_a-(-myEF_r_ran_p))^2)+1), "max"),
  "PSmag_random_sd" = decostand(log(sqrt((-myEF2_r_sd_a-(-myEF_r_ran_p))^2)+1), "max"),
  "PSmag_proc_range" = decostand(log(sqrt((-myEF2_pro_ran_a-(-myEF_pro_ran_p))^2)+1), "max"),
  "PSmag_proc_sd" = decostand(log(sqrt((-myEF2_pro_sd_a-(-myEF_pro_sd_p))^2)+1), "max"))

PSmag <- PSmag[match(colnames(CoRnIn_x), rownames(PSmag)),]







plot(PSdir$PSdir_process_sd,PSmag$PSmag_proc_sd)
dev.off()

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
names(Zscores) <- unique(pool$site)

# Calculate effect sizes from null networks 

# modularity 
CnRnIn_SES_mod <- (Zscores-CnRnIn_x[1,]) / CnRnIn_sd[1,]
CoRoIn_SES_mod <- (Zscores-CoRoIn_x[1,]) / CoRoIn_sd[1,]
CnRoIn_SES_mod <- (Zscores-CnRoIn_x[1,]) / CnRoIn_sd[1,]
CoRnIn_SES_mod <- (Zscores-CoRnIn_x[1,]) / CoRnIn_sd[1,]

# nestedness
CnRnIn_SES_nes <- (Zscores-CnRnIn_x[2,]) / CnRnIn_sd[2,]
CoRoIn_SES_nes <- (Zscores-CoRoIn_x[2,]) / CoRoIn_sd[2,]
CnRoIn_SES_nes <- (Zscores-CnRoIn_x[2,]) / CnRoIn_sd[2,]
CoRnIn_SES_nes <- (Zscores-CoRnIn_x[2,]) / CoRnIn_sd[2,]



## Partialying out separate effects of consumer, resources, and interactions 

# modularity 

INTef_SES_mod <- CnRnIn_SES_mod + CoRnIn_SES_mod -CnRoIn_SES_mod -  CoRoIn_SES_mod
INT_SES_mod <- CoRoIn_SES_mod
R_SES_mod <- CoRnIn_SES_mod- CoRoIn_SES_mod
C_SES_mod <- CnRoIn_SES_mod- CoRoIn_SES_mod

# nestedness 
INTef_SES_nes <- CnRnIn_SES_nes + CoRnIn_SES_nes -CnRoIn_SES_nes -  CoRoIn_SES_nes
INT_SES_nes <- CoRoIn_SES_nes
R_SES_nes <- CoRnIn_SES_nes- CoRoIn_SES_nes
C_SES_nes <- CnRoIn_SES_nes- CoRoIn_SES_nes




SESnet_mod <- data.frame(
  INTef_SES_mod, 
  INT_SES_mod,
  R_SES_mod,
  C_SES_mod)

#SESnet_mod <- SESnet_mod/max(abs(SESnet_mod))
colnames(SESnet_mod) <- c("InterEf", "Interactions", "Resources", "Consumers")


apply(SESnet_mod, 2, Rmisc::CI, 0.95)



SESnet_nes <- data.frame(
  INTef_SES_nes, 
  INT_SES_nes,
  R_SES_nes,
  C_SES_nes)

#SESnet_nes <- SESnet_nes/max(abs(SESnet_nes))
colnames(SESnet_nes) <- c("InterEf", "Interactions", "Resources", "Consumers")

apply(SESnet_nes, 2, Rmisc::CI, 0.95)















Rmisc::CI(CnRnIn_SES_mod, 0.95)
#####
# Calculate PSDir and PSMag
########
psmag <- vegan::decostand(log(sqrt((myEF2_pro_sd_a-myEF_pro_sd_p)^2)+1), "max")
psdir <- log(-myEF2_pro_sd_a^2/-myEF_pro_sd_p^2) / max(abs(log(-myEF2_pro_sd_a^2/-myEF_pro_sd_p^2)))
########
## Figure 4:Effect of process strength asymmetry on plant-frugivore network structure along an elevational gradient
#######
#### Figure 4.a) Strength of niche-based community assembly and assymetry processes for plants and frugivores across elevation
#############

par(mfrow = c(1,2), las = 1, oma = c(1,2,1,1),mar  = c(2,2,2,2))
layout(matrix(c(1,1,1,2), 1,4))
elev_vec <- dataSet1$N$elevation[match(names(myEF_pro_sd_p),dataSet1$N$site)]
elev_vec1 <- c(elev_vec-c(-150, -150, +250, +350, +250, -150))
elev_vec2 <- c(3250,1250, 2850, 750, 1850, 2250)

plot(-myEF_pro_sd_p,elev_vec1, frame= F,
     xlim = c(-3,3), pch = 16,
     ylim  = c(0,4000),xaxt = "n",xlab = "", ylab = "",
     col = scales::alpha("red", 0.7), cex = 1)
axis(1, c(-3:8), c(-3:8))
segments(0,elev_vec1,-myEF_pro_sd_p,elev_vec1, col = "red", lwd = 2)
points(-myEF_pro_sd_p,elev_vec1, cex =1)
abline(v=0)
abline(h = c(500,1500, 2500,3500), lty = 2)
points(-myEF2_pro_sd_a,c(elev_vec2),
       ylim = c(-3,8), 
       xlim  = c(0,4000),
       col = scales::alpha("skyblue", 0.7),
       pch = 16, cex = 1)
segments(0,elev_vec2,-myEF2_pro_sd_a,elev_vec2, col = "skyblue" ,lwd = 2)
points(-myEF2_pro_sd_a,c(elev_vec2),
       cex = 1, pch = 1)
##########
# Figure 4b,c: Niche-based effects of plants -b- and frugivores -c- on network structure as a function of process strength assymetry 
##########

par(mfrow = c(1,2), las = 1, oma =  c(1,2,1,1),mar  = c(2,2,2,2))
plot(SESnet_mod$Resources[match(names(psdir), rownames(SESnet_mod))]
     ~abs(psdir), frame = F,
     xlim = c(0,1),
     pch = 16, col = "orange",
     xlab = "", 
     ylim = c(-20,20),
     ylab = "",
     cex = 2)
points(SESnet_mod$Resources[match(names(psdir), rownames(SESnet_mod))]
       ~abs(psdir), cex = 2)
points(SESnet_nes$Resources[match(names(psdir), rownames(SESnet_mod))]
       ~abs(psdir),
       pch = 16, col = "purple",
       cex = 2)
points(SESnet_nes$Resources[match(names(psdir), rownames(SESnet_mod))]
       ~abs(psdir),cex =2)
abline(h = 0, lty = 1, col = "grey30", lwd = 2)

# 
# abline(lm(SESnet_mod$Resources~PSdir$PSdir_process_sd), col = "orange", lty = 1, lwd = 3)
# abline(lm(SESnet_nes$Resources~PSdir$PSdir_process_sd), col = "purple" ,lty = 2, lwd = 3)

# legend(0,0.2, bty = "n", legend = c("Modularity: R=0.39, p = 0.13"))
# legend(0,0.1, bty = "n", legend = c("Nestedness: R=0.25, p = 0.24"))

# summary(lm(SESnet_mod$Resources~PSdir$PSdir_process_sd))
# summary(lm(SESnet_nes$Interactions~PSdir$PSdir_process_sd))

plot(SESnet_mod$Consumers[match(names(psdir), rownames(SESnet_mod))]~abs(psdir), 
     # ylim = c(-1,1), 
     pch = 16, frame = F,
     col = "orange",
     ylim = c(-20,20),
     xlab = "", 
     ylab = "",
     cex = 2)
points(SESnet_mod$Consumers[match(names(psdir), rownames(SESnet_mod))]~abs(psdir), cex = 2)
points(SESnet_nes$Consumers[match(names(psdir), rownames(SESnet_mod))]~abs(psdir),
       pch = 16, col = "purple",
       cex = 2)
points(SESnet_nes$Consumers[match(names(psdir), rownames(SESnet_mod))]~abs(psdir),cex=2)
abline(h = 0, lty = 1, col = "grey30", lwd = 2)


abline(lm(SESnet_mod$Consumers~PSdir$PSdir_process_sd), col = "orange", lty = 1, lwd = 3)
abline(lm(SESnet_nes$Consumers~PSdir$PSdir_process_sd), col = "purple" ,lty = 2, lwd = 3)
# legend(0,0.2, bty = "n", legend = c("Modularity: R=0.77, p = 0.02"))
# legend(0,0.1, bt

plot(abs(psdir),elev_vec2,frame= F, xaxt = "n", yaxt = "n",
     xlim = c(0,1.4), pch = 0, cex = 1, 
     ylim  = c(0,4000),xaxt = "n",xlab = "", ylab = "")
segments(0,elev_vec2,abs(psdir),elev_vec2, col = "gray" ,lwd = 2)
points(abs(psdir),elev_vec2,cex = 1, col = ifelse(psdir <0, "firebrick1","deepskyblue" ),pch = 15)
points(abs(psdir),elev_vec2,cex = 1, pch = 0)
abline(h = c(500,1500, 2500,3500), lty = 2)
abline(v = 0)
#############





