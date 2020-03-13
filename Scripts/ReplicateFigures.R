### Replicate figures presented at:
## DOI:

#######
## Gabriel Muñoz 
# Mar 2020
#########

source("CaseStudy.R")


######
# Case study figure (Figure 3)
######

### Make the plot panels 

png("Figs/FigsTemp/CaseStudyNullperElev.png", width = 1500, height = 3000, pointsize = 60)
par(mfrow = c(3,1), las = 1)
for(i in c(1,3,5)){
  colR <- scales::alpha("darkgreen", 0.7)
  colB <- scales::alpha("black", 0.7)
  ######
  # Make the main plot
  ######
  plot(bindEFCaseStudy$AnBnIn.mod[i]~bindEFCaseStudy$AnBnIn.nes[i],
       xlim = c(0.01,6),
       ylim = c(0.01,15), 
       frame = F,
       main = paste("Elevation = ", bindEFCaseStudy$elev[i]),
       xlab = "SES Nestedness",
       ylab = "SES Modularity",
       pch = "+",
       col = colB,
       cex = 2)
  points(bindEFCaseStudy$AoBnIn.mod[i]~bindEFCaseStudy$AoBnIn.nes[i],
         pch = 15,   
         col = colB, cex = 2)
  points(bindEFCaseStudy$AnBoIn.mod[i]~bindEFCaseStudy$AnBoIn.nes[i],
         pch = 16,
         col = colB,
         cex = 2)
  points(bindEFCaseStudy$AoBoIn.mod[i]~bindEFCaseStudy$AoBoIn.nes[i],
         pch = 17,     
         col = colB,
         cex = 2)
  ######## 
  # add nice border to points
  #######
  points(bindEFCaseStudy$AoBnIn.mod[i]~bindEFCaseStudy$AoBnIn.nes[i],
         pch = 0, cex = 2)
  points(bindEFCaseStudy$AnBoIn.mod[i]~bindEFCaseStudy$AnBoIn.nes[i],
         pch = 1, cex = 2)
  points(bindEFCaseStudy$AoBoIn.mod[i]~bindEFCaseStudy$AoBoIn.nes[i],
         pch = 2, cex = 2)
  #####
  # Second network at the same elevation
  ######
  # Make the main plot
  ######
  points(bindEFCaseStudy$AnBnIn.mod[i+1]~bindEFCaseStudy$AnBnIn.nes[i+1],
         pch = "+", cex = 2, col = colR)
  points(bindEFCaseStudy$AoBnIn.mod[i+1]~bindEFCaseStudy$AoBnIn.nes[i+1],
         pch = 15, cex = 2, col = colR)
  points(bindEFCaseStudy$AnBoIn.mod[i+1]~bindEFCaseStudy$AnBoIn.nes[i+1],
         pch = 16, cex = 2, col = colR)
  points(bindEFCaseStudy$AoBoIn.mod[i+1]~bindEFCaseStudy$AoBoIn.nes[i+1],
         pch = 17, cex = 2, col = colR)
  ######## 
  # add nice border to points
  #######
  points(bindEFCaseStudy$AoBnIn.mod[i+1]~bindEFCaseStudy$AoBnIn.nes[i+1],
         pch = 0, cex = 2)
  points(bindEFCaseStudy$AnBoIn.mod[i+1]~bindEFCaseStudy$AnBoIn.nes[i+1],
         pch = 1, cex = 2)
  points(bindEFCaseStudy$AoBoIn.mod[i+1]~bindEFCaseStudy$AoBoIn.nes[i+1],
         pch = 2, cex = 2)
  ######
  # add legend
  #####
  legend("topleft", 
         fill = c("black", "darkgreen"), 
         bty = "n",
         legend = c(as.character(bindEFCaseStudy$AnBnIn.site[i]), 
                    as.character(bindEFCaseStudy$AnBnIn.site[i+1])))
  legend("topright", 
         pch = c(3, 0, 1, 2), 
         bty = "n",
         legend = c("CnRnIn", "CoRnIn", "CnRoIn", "CoRoIn"))
}
dev.off()

### load the mountain figure 
library(magick)


Mountain <- image_read("Figs/FigsTemp/Screen Shot 2020-03-13 at 11.06.31.png")
panel <- image_read("Figs/FigsTemp/CaseStudyNullperElev.png")

comb <- c(image_scale(Mountain, "2400x"), panel)
comb <- image_append(comb)
image_write(comb, path = "Figs/FigsTemp/nullMount.png", format = "png")
