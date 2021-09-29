




tiff(filename="toPress/SimulRes2.tiff",
     height=5600,width=5200,
     units="px",res=800,
     pointsize = 15,
     compression="lzw")
par(mfrow = c(1,1), 
    oma = c(1,1,1,1),
    mar = c(4,5,3,2), las = 1)

plot(NODFx~QZscorex,
     frame = F,
     xlim = c(40, 150),
     ylim = c(0,200),
     pch = ifelse( RES_DD$intHyp == "NL", 21, 
                   ifelse( RES_DD$intHyp == "FL", 22,
                           24)),
     col = scales::alpha(ifelse( RES_DD$intHyp == "NL", "#FADE43", 
                                 ifelse( RES_DD$intHyp == "FL", "#E0A738",
                                         "#9C413D")),1),
     xlab = "Modularity (Q Z-score)", 
     ylab = "Nestedness (NODF Z-score)", 
     cex.lab = 1.5,
     cex.axis = 1,
     cex = 1, 
     data = RES_DD[!RES_DD$EA == "ND",])


Plot_ConvexHull( RES_DD$QZscorex[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL")],
                 RES_DD$NODFx[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL")], 
                 scales::alpha("#FADE43",1), lwd = 3, lty = 2)

Plot_ConvexHull( RES_DD$QZscorex[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM")],
                 RES_DD$NODFx[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM")], 
                 scales::alpha("#9C413D",1), lwd = 3, lty = 2)

Plot_ConvexHull( RES_DD$QZscorex[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL")],
                 RES_DD$NODFx[!RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL")], lcolor = 
                   scales::alpha("#E0A738",1), lwd = 3, lty = 2)

#### 
# neutral arrows

segments("x0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$QZscorex
         -RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$QZscoresd,
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$QZscorex
         +RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$QZscoresd,
         "y0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("NL") ,]$NODFx, 
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$NODFx,
         lwd = 3, lty = 1, col = scales::alpha("#FADE43",1))

segments("y0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$NODFx
         -RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("NL"),]$NODFsd,
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$NODFx
         +RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("NL"),]$NODFsd,
         "x0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$QZscorex, 
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("NL") ,]$QZscorex, 
         lwd = 3, lty = 1,col = scales::alpha("#FADE43",1))

segments("x0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$QZscorex
         -RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$QZscoresd,
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$QZscorex
         +RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$QZscoresd,
         "y0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("MM") ,]$NODFx, 
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$NODFx,
         lwd = 3, lty = 1, col = scales::alpha("#9C413D",1))

segments("y0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$NODFx
         -RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("MM"),]$NODFsd,
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$NODFx
         +RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("MM"),]$NODFsd,
         "x0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$QZscorex, 
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("MM") ,]$QZscorex, 
         lwd = 3, lty = 1,col = scales::alpha("#9C413D",1))

segments("x0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$QZscorex
         -RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$QZscoresd,
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$QZscorex
         +RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$QZscoresd,
         "y0" = RES_DD[RES_DD$EA == "ND"& RES_DD$intHyp %in% c("FL") ,]$NODFx, 
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$NODFx,
         lwd = 3, lty = 1, col = scales::alpha("#E0A738",1))

segments("y0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$NODFx
         -RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("FL"),]$NODFsd,
         "y1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$NODFx
         +RES_DD[RES_DD$EA == "ND"  & RES_DD$intHyp %in% c("FL"),]$NODFsd,
         "x0" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$QZscorex, 
         "x1" = RES_DD[RES_DD$EA == "ND" & RES_DD$intHyp %in% c("FL") ,]$QZscorex, 
         lwd = 3, lty = 1,col = scales::alpha("#E0A738",1))
###################


legend("topright", 
       title = "Interaction assembly type",
       legend = c("Morphological matching", 
                  "Forbidden links",
                  "Stochastic interactions"),
       pch = c(2,0,1),
       cex = 0.7,
       bty = "n",
       col = c("#9C413D","#E0A738","#FADE43"))

legend("bottomleft", 
       title = "Community assembly",
       legend = c("Environmental filtering", 
                  "Neutral assembly"),
       lty = c(2,1), lwd  = c(1,1),
       cex = 0.7, 
       bty = "n",
       col = c("gray80","black"))

dev.off()










plot((butd_eff$q_eff)[match(RES_ass$key,butd_eff$key)]~c(assdev/max(abs(assdev))), 
     xlab = "PSA effect on modularity",
     ylab = "Network assembly mode",
     cex.lab = 1.5,
     frame =F,
     xlim = c(-1,1),
     ylim = c(-6,6),
     pch = ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "NL", 21, 
                   ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "FL", 22,
                           24)),
     bg = scales::alpha(ifelse(butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "NL", "#FADE43", 
                               ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "FL", "#E0A738",
                                       "#9C413D")),0.5),
     col = scales::alpha(ifelse(butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "NL", "#FADE43", 
                               ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "FL", "#E0A738",
                                       "#9C413D")),0.9))

abline(h=0,v=0, lty = 2)

legend("topleft", "Top-down", bty = "n")
legend("bottomleft", "Bottom-up", bty = "n")

legend("bottomright", 
       title = "Interaction assembly type",
       legend = c("Morphological matching", 
                  "Forbidden links",
                  "Stochastic interactions"),
       pch = c(2,0,1)+15,
       cex = 0.7,
       bty = "n",
       col = c("#9C413D","#E0A738","#FADE43"))

mtext("B",3, outer = T, adj = 0.35 , cex = 2)


#########

plot((butd_eff$nodf_eff)[match(RES_ass$key,butd_eff$key)]~c(assdev2/max(abs(assdev2))), 
     col = scales::alpha(ifelse(butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "NL", "#FADE43", 
                                ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "FL", "#E0A738",
                                        "#9C413D")),1),
     bg = scales::alpha(ifelse(butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "NL", "#FADE43", 
                               ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "FL", "#E0A738",
                                       "#9C413D")),0.5),
     xlab = "PSA effect on nestedness",
     ylab = "Network assembly mode",
     cex.lab = 1.5,
     frame =F,
     xlim = c(-1,1),
     ylim = c(-6,6),
     pch = ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "NL", 21, 
                   ifelse( butd_eff$intHyp[match(RES_ass$key,butd_eff$key)] == "FL", 22,
                           24)))

abline(h=0,v=0, lty = 2)

legend("topleft", "Top-down", bty = "n")
legend("bottomleft", "Bottom-up", bty = "n")

legend("bottomright", 
       title = "Interaction assembly type",
       legend = c("Morphological matching", 
                  "Forbidden links",
                  "Stochastic interactions"),
       pch = c(2,0,1)+15,
       cex = 0.7,
       bty = "n",
       col = c("#9C413D","#E0A738","#FADE43"))
mtext("C",3, outer = T, adj = 0.35*2 , cex = 2)

dev.off()



