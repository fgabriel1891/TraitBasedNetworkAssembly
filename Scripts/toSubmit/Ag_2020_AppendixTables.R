###################################
## Code to replicate the appendices shown in: Trait-based inference of ecological network assembly
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

#####
# Figure S1. Relationships of mean modularity and nestedness with varying strengths of environmental filtering asymmetry directions and magnitudes.
######
########

layout(mat = matrix(c(1,2,3,4,5,5), 3,2, T), heights = c(0.6,0.6,0.2))
par(mar = c(5,6,2,2))
plot(RES2$QZscorex~RES2$PSSdir,
     col = scales::alpha(as.numeric(RES2$intHyp) + 3, 0.7), 
     cex = 2,
     xlab = "Process strength asymmetry (Direction)",
     ylab = " Mean Q Z-score",
     pch = 16 , 
     cex.lab = 1.5)
points(RES2$QZscorex~RES2$PSSdir, cex = 2)

plot(RES2$NODFx~RES2$PSSdir,
     col = scales::alpha(as.numeric(RES2$intHyp) + 3, 0.7), 
     cex = 2,
     xlab = "Process strength asymmetry (Direction)",
     ylab = "Mean NODF Z-score",
     pch = 16 , 
     cex.lab = 1.5)
points(RES2$NODFx~RES2$PSSdir, cex = 2)

plot(RES2$QZscorex~RES2$PSSmag,
     col = scales::alpha(as.numeric(RES2$intHyp) + 3, 0.7), 
     cex = 2,
     cex.lab = 1.5,
     xlab = "Process strength asymmetry (Magnitude)",
     ylab = " Mean Q Z-score",
     pch = 16 )
points(RES2$QZscorex~RES2$PSSmag, cex = 2)

plot(RES2$NODFx~RES2$PSSmag,
     col = scales::alpha(as.numeric(RES2$intHyp) + 3, 0.7), 
     cex = 2,
     cex.lab = 1.5,
     xlab = "Process strength asymmetry (Magnitude)",
     ylab = " Mean NODF Z-score",
     pch = 16 )
points(RES2$NODFx~RES2$PSSmag, cex = 2)
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top", legend = c("Neutral interactions", "Morphological Matching", "Forbidden links"),
       col = as.numeric(RES2$intHyp)+3, pch  = c(16,16,16), horiz = T, cex=1.3)

legend("top", legend = c("Neutral interactions", "Morphological Matching", "Forbidden links"),
       pch  = c(1,1,1), bty = "n", horiz = T, cex = 1.3)

###
# Figure S2. Variability of modularity and nestedness with varying strengths of environmental filtering asymmetry
########
## Z-score variances
#####

layout(mat = matrix(c(1,2,3,4,5,5), 3,2, T), heights = c(0.6,0.6,0.2))
par(mar = c(5,6,2,2))
plot(RES2$QZscoresd~RES2$PSSdir,
     col = scales::alpha(as.numeric(RES2$intHyp) + 3, 0.7), 
     cex = 2,
     xlab = "Process strength symmetry (Direction)",
     ylab = " SD Q Z-score",
     pch = 16 , 
     cex.lab = 1.5,
     log = "y")
points(RES2$QZscoresd~RES2$PSSdir, cex = 2)

plot(RES2$QZscoresd~RES2$PSSmag,
     col = scales::alpha(as.numeric(RES2$intHyp) + 3, 0.7), 
     cex = 2,
     xlab = "Process strength symmetry (Magnitude)",
     ylab = "SD Q Z-score",
     pch = 16 , 
     cex.lab = 1.5,
     log = "y")
points(RES2$QZscoresd~RES2$PSSmag, cex = 2)

plot(RES2$NODFsd~RES2$PSSdir,
     col = scales::alpha(as.numeric(RES2$intHyp) + 3, 0.7), 
     cex = 2,
     cex.lab = 1.5,
     xlab = "Process strength symmetry (Direction)",
     ylab = " SD NODF Z-score",
     pch = 16 )
points(RES2$NODFsd~RES2$PSSdir, cex = 2)

plot(RES2$NODFsd~RES2$PSSmag,
     col = scales::alpha(as.numeric(RES2$intHyp) + 3, 0.7), 
     cex = 2,
     cex.lab = 1.5,
     xlab = "Process strength symmetry (Magnitude)",
     ylab = " SD NODF Z-score",
     pch = 16 )
points(RES2$NODFsd~RES2$PSSmag, cex = 2)
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend("top", legend = c("Neutral interactions", "Morphological Matching", "Forbidden links"),
       col = as.numeric(RES2$intHyp)+3, pch  = c(16,16,16), horiz = T, cex=1.3)

legend("top", legend = c("Neutral interactions", "Morphological Matching", "Forbidden links"),
       pch  = c(1,1,1), bty = "n", horiz = T, cex = 1.3)


#####################

#####
# Table S2a. Relative contribution of community and interaction assembly process to the variance in network structural metrics
######

## Modularity 

print(knitr::kable(
  anova(lm(QZscorex~
             + EA
           + EA * intHyp 
           + intHyp, RES_DD)
  )))

print(knitr::kable(
  summary(lm(NODFx~
               + EA
             + EA * intHyp 
             + intHyp, RES_DD)
  )))

##########
######
# Commonality coefficients 
######
apsOut=yhat::commonalityCoefficients(RES_DD,"QZscorex",
                                     list(
                                       "EA",
                                       "intHyp"))
print(apsOut)

apsOut=yhat::commonalityCoefficients(RES_DD,"NODFx",
                                     list(
                                       "EA",
                                       "intHyp"))

print(apsOut)
##########
#####
# Table S2b. Relative contributions of niche-based network assembly parameters to the variation in network metrics at distinct interaction assembly procesess 
#######


iH <- c("FL", "MM", "NL")



# Q means
for(i in 1:3){ 
  print(iH[i])
  #print(anova(lm(QZscorex~sigmaDif+sigmaA+sigmaB+sigmaDif, RES[RES$intHyp == iH[i],])))
  print(knitr::kable(
    anova(lm(QZscorex~abs(PSSdir)+PSSmag, RES2[RES2$intHyp == iH[i],]
    ))))
}
# Q means
for(i in 1:3){ 
  print(iH[i])
  # print(anova(lm(NODFx~sigmaDif+sigmaA+sigmaB, RES[RES$intHyp == iH[i],])))
  print(summary(lm(QZscorex~abs(PSSdir)+PSSmag , RES2[RES2$intHyp == iH[i],])))
}

# NODF means
for(i in 1:3){ 
  print(iH[i])
  # print(anova(lm(NODFx~sigmaDif+sigmaA+sigmaB, RES[RES$intHyp == iH[i],])))
  print(knitr::kable(anova(
    lm(NODFx~abs(PSSdir)+PSSmag , RES2[RES2$intHyp == iH[i],]
    ))))
}

# NODF means
for(i in 1:3){ 
  print(iH[i])
  # print(anova(lm(NODFx~sigmaDif+sigmaA+sigmaB, RES[RES$intHyp == iH[i],])))
  print(summary(lm(NODFx~PSSdir+PSSmag , RES2[RES2$intHyp == iH[i],])))
}


# Variance partitioning
## Commonality coefficients 
RES2$dir_abs <- abs(RES2$PSSdir)
# Q means
for(i in 1:3){ 
  print(iH[i])
  apsOut=yhat::commonalityCoefficients(RES2[RES2$intHyp == iH[i],],"QZscorex",
                                       list(
                                         "dir_abs","PSSmag"))
  print(apsOut)
}




# NODF mean
for(i in 1:3){ 
  print(iH[i])
  apsOut=yhat::commonalityCoefficients(RES2[RES2$intHyp == iH[i],],"NODFx",
                                       list(
                                         "dir_abs","PSSmag"))
  print(apsOut)
}


############




