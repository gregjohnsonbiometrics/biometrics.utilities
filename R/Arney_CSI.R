#-------------------------------------------------------------------------------
# Arney's Competitive Stress Index
#
# version 20140104
#-------------------------------------------------------------------------------

# This is a function that compute Arney's crown overlap and Competitive
# Stress Index as described in Arney (1973).

# Arney, J.D. 1973. Tables for quantifying competitive stress on individual trees.
# Pacific Forest Research Centre, Canadian Forest Service, Victoria, BC. 
# Information Report BC-X-78. 47p.

###
### function to compute maximum crown width using Arney (1973) for testing
###

# This function calculates maximum crown width (MCW) using the combined
# equation for Douglas-fir from Arney (1973). It is used for testing and in
# application MCW should use species specific equations.

# Inputs: 
#    DBH = indivicual tree diameter at breast height (inches)

# Output:
#    MCW = maximum crown width for a trees (feet)

crown.width.Arney <- function(DBH)
{
   #CW <- 1.9118 + 2.8880*DBH - 0.0707*DBH^2    # BC equation
   #CW <- 4.7071 + 2.0168*DBH - 0.0186*DBH^2    # OR equation
   CW <- 4.0223 + 2.1228*DBH - 0.0220*DBH^2    # combined
   return(CW)
}

###
### function to compute Arney's (1973) crown area overlap (AO) and competitive stress index (CSI)
###

# This function calculates corwn area overlap (AO) and Competitive Stress Index (CSI)
# using the the equation described by Arney (1973).

# Inputs:
#   tempT    data frame with a tree-list a subject tree and its competitors however
#            defined (e.g. maximum distance from subject tree). 
#               DBH = diameter at breast height (inches)
#               DIST = the distance between subject tree and competitors. The 
#                        DIST for the subject tree is 0.
#               MCW = (OPTIONAL) the species specific maximum crown width. If not
#                   provided it will be predicted using the equation forDouglas-fir
#                   from Arney (1973). This is mostly for testing and providing
#                  MCW is highly recommended.
#            In addition, the tree-list can have other variables 
#            as desired (e.g. tree number and species).
#   CSIflag  TRUE reports (DEFAULT) the CSI and FALSE just reports the area overlaps

# Outputs the input data frame with the following columns added:
#     CA = the maximum crown area for the trees (for testing)
#     Type = equation type (for testing)
#        0 = distance is 0 or trees do not intersect
#        1, 2, 3, 4 = equations 1, 2, 3, 4
#        5 = subject tree
#     AO = overlap area (ft2)
#     CSI = tree competitive stress index (if CSIflag = TRUE)

# If MCW is not provided for each tree they are computed using the equation for 
# Douglas-fir reported by Arney (1973). This is primarily for testing and it is 
# recommended that species specific MCWs should provided.

# Waring, this does not currently this does not take into account edge effects.

# to compute CSI for a tree-list (see test 6 below)
#   tempT <- data.frame(tree=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),
#                       DBH=c(4,4,3,3,3,5,5,4,5,4,3,6,4,4),
#                       DIST=c(7,13,10,4,7,12,10,12,7,13,5,9,13,0))
#   notes:
#    (1) includes DBH and DIST from subject tree and competitor.
#    (2) the distances (DIST) are the distance between the subject and competitors.
#    (3) the DIST = 0 for the subject tree
#    (4) as give the MCWs will be computed, but is is recommented that a column
#        be added to the tempT data frame with computed MCWs that are species specific. 
#   tempTx <- Arney.CSI(tempT, CSIflag=TRUE)
#   sum(tempTx$CSI) + 100           # 332.4
#   # The 100 is 100*AOsubject/AOsubject (could be added to the code if wanted)

# to compute average CSI for a variable plot tree-list (see test 7 below)
#   tempT <- data.frame(tree=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
#                       DBH=c(11,11,12,12,12,13,14,15,15,15,16,16,17,17,18),
#                       DIST=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0).
#                       EXPAN=c(30.30521,30.30521,25.46479,25.46479,25.46479,21.69781,
#                               18.70883,16.29747,16.29747,16.29747,14.32394,14.32394,
#                               12.68834,12.68834,11.31768))
#   notes
#    (1) the addition of an optional variable EXPAN which is the trees per acre
#        of each tree ( EXPAN = (BAF/(0.005454154*DBH^2))/NPTS
#    (2) all distances are 0 which returns the crown areas.
#    (3) as give the MCWs will be computed, but is is recommented that a column
#        be added to the tempT data frame with computed MCWs that are species specific. 
#   NTrees <- length(tempT$DBH)
#   TPA <- sum(tempT$EXPAN)
#   BAPA <- sum(tempT$EXPAN * 0.005454154 * tempT$DBH^2)
#   tempTx <- Arney.CSI(tempT, CSIflag=FALSE)
#   note that the CSIlag if FALSE since we are not using distances but an average here
#   CSIavg < 100 + sum(100*((tempTx$CA/43560)*tempTx$EXPAN))

Arney.CSI <- function(tempT, CSIflag=TRUE)
{
   # check for MCW
   if (!("MCW" %in% names(tempT)))
   {
      tempT$MCW <- 0
   }
   # deal with missing MCWs (this is mostly for testing)
   tempT$MCW <- ifelse (is.na(tempT$MCW) == TRUE | tempT$MCW == 0, 
                           crown.width.Arney(tempT$DBH), tempT$MCW)
    
   # compute cown area

   tempT$CA <- 0
   for (i in 1:nrow(tempT))
   {
      tempT$CA[i] <- pi*(tempT$MCW[i]/2)^2
   }

   # just get crown area for average CSI for samples

   if (length(tempT$DIST[tempT$DIST == 0]) > 1)
   {
      return(tempT)
   }

   # compute area overlap

   tempT$Type <- 0
   tempT$AO <- 0
   for (i in 1:nrow(tempT))
   {
      # arrange the largest and smallest tree for AO and CSI
      if (tempT$DBH[tempT$DIST == 0] >= tempT$DBH[i])
      {
         # largest tree is subject tree
         CW1 <- tempT$MCW[tempT$DIST == 0]
         DBH1 <- tempT$DBH[tempT$DIST == 0]
         r1 <- CW1/2   # crown radius
         # smaller tree is competitor
         CW2 <- tempT$MCW[i]
         DBH2 <- tempT$DBH[i]
         r2 <- CW2/2   # crown radius
      } else {
         # the largest tree is the competitor
         CW1 <- tempT$MCW[i]
         DBH1 <- tempT$DBH[i]
         r1 <- CW1/2   # crown radius
         # smaller tree is the subject tree
         CW2 <- tempT$MCW[tempT$DIST == 0]
         DBH2 <- tempT$DBH[tempT$DIST == 0]
         r2 <- CW2/2   # crown radius
      }
   
      # distance between subject tree to competitor
      DIST <- tempT$DIST[i]
      # compute tree crown area
      #tempT$CA[i] <- pi*(tempT$MCW[i]/2)^2

      # if this is the subject tree
      if (is.na(DIST) == TRUE | DIST <= 0)
      {
         #AO <- 0   # already set to 0
         tempT$Type[i] <- 5
         next
      }

      # check if trees actually overlap
      if (r1 + r2 < DIST)
      {
         #AO <- 0   # already set to 0
         tempT$Type[i] <- 0
         next
      }

      # r1 is largest and totally overlaps r2
      if (r1 > r2 & r1 > r2+DIST)
      {
         # arrangement 4: smaller tree is completely overlapped by the larger tree
         # (Equation 4 Appendix 1 page 5)
         tempT$AO[i] <- pi*r2^2
         tempT$Type[i] <- 4
         next
      }

      s <- (r1 + r2 + DIST)/2
      c <- (2/DIST)*sqrt(s*(s-r1)*(s-r2)*(s-DIST))
      x1 <- sqrt(r1^2 - c^2)

      DIST <- round(DIST,0)
      x1 <- round(x1,0)
      
      # arrangement 1: DIST > x1 (Equation 1 Appendix 1 page 5)
      if (DIST > x1)
      {
         tempT$AO[i] <- r1^2*asin(c/r1) + r2^2 * asin(c/r2) - DIST*c  # checks
         tempT$Type[i] <- 1
         next
      }

      # arrangement 2: DIST = x1 (Equation 2 Appendix 1 page 5)
      if (DIST == x1)
      {
         #tempT$AO[i] <- pi*r2^2 + r1^2*asin(r2/r1) - DIST*r2
         # PROBLEM -- I don't think this is correct
         #r1 <- 10; r2 <- 6
         # area of half competitor + sector
         #theda <- acos((r1^2 + r1^2 - (2*r2)^2)/(2*r1*r1))
         #tempT$AO[i] <- 0.5*pi*r2^2 + (0.5*r1^2*(theda - sin(theda)))
         # I think this is what was intended
         tempT$AO[i] <- (0.5*pi*r2^2) + (r1^2*asin(r2/r1)) - (DIST*r2)
         tempT$Type[i] <- 2
         next
      }

      if (DIST < x1)
      {
         # arrangement 3: DIST < x1 but competitor is not completely overlapped
         # (Equation 3 Appendix 1 page 5)
         x2 <- x1 - DIST
         tempT$AO[i] <- pi*r2^2 - r2^2*asin(c/r2) + x2*c + r1^2 * asin(c/r1) - x1*c
         tempT$Type[i] <- 3
         next
      }
   }

   # compute CSI -- #  CSI <- 100*aij/Aj
   
   if (CSIflag == TRUE)
   {
      tempT$CSI <- ifelse(tempT$AO == 0, 0, 100*tempT$AO/(pi*r1^2))
   }

   # return 
   return(tempT)
}

### Testing --------------------------------------------------------------------

# independent function for computing area of interections for 2 circles
#source("C:/Users/DavidFolder/WorkC/R_Packages/myGeometryR/R/circle_circle_area_intersection.R")

#R1 <- crown.width.Arney(DBH=tempT$DBH[2])/2  # larget tree
#R2 <- crown.width.Arney(DBH=tempT$DBH[1])/2  # smaller tree
#100*130/(pi*R1^2)

#circle.circle.area.intersection(Xa=0, Ya=0, Ra=R1, Xb=tempT$DIST[1], Yb=0, Rb=R2, plot.it=TRUE)  # 39.2266

#minX <- 0 - tempT$DIST[1] - R1 - R2
#maxX <- 0 + tempT$DIST[1] + R1 + R2
#minY <- 0 - tempT$DIST[1] - R1 - R2
#maxY <- 0 + tempT$DIST[1] + R1 + R2

#plot(x=c(0,tempT$DIST[1]), y=c(0,0), asp=1.0, pch=c("+"),
#    xlim=c(minX,maxX), ylim=c(minY,maxY), 
#    xlab="X Coordinates", ylab="Y Coordinates", sub="Red=Subject / Green=Competitor")
#symbols(x=0, y=0, circles=R1, add=TRUE, inches=FALSE, fg="red")
#symbols(x=tempT$DIST[1], y=0, circles=R2, add=TRUE, inches=FALSE, fg="green")

## Test 1 -- trees do not intersect (CHECK)

tempT <- data.frame(DBH=c(4,4), DIST=c(13,0))
out <- Arney.CSI(tempT, CSIflag=TRUE)
round(out$AO[1],0) == 0
round(out$CSI[1],0) == 0

#R1 <- crown.width.Arney(DBH=tempT$DBH[1])/2  # larget tree
#R2 <- crown.width.Arney(DBH=tempT$DBH[2])/2  # smaller tree
#circle.circle.area.intersection(Xa=0, Ya=0, Ra=R1, Xb=tempT$DIST[1], Yb=0, Rb=R2, plot.it=TRUE)  # 35.99891

### Test 2 -- DIST > x1 (Equation 1 Appendix 1 page 5) (CHECK)

tempT <- data.frame(DBH=c(4,4), DIST=c(7,0))
out <- Arney.CSI(tempT, CSIflag=TRUE); out
round(out$AO[1],0) == round(35.99891,0)
round(out$CSI[1],0) == 31

#R1 <- crown.width.Arney(DBH=tempT$DBH[1])/2  # larget tree
#R2 <- crown.width.Arney(DBH=tempT$DBH[2])/2  # smaller tree
#circle.circle.area.intersection(Xa=0, Ya=0, Ra=R1, Xb=tempT$DIST[1], Yb=0, Rb=R2, plot.it=TRUE)  # 35.99891

### Test 3 -- DIST = x1 (Equation 2 Appendix 1 page 5) (CHECK)

tempT <- data.frame(DBH=c(3,5.102), DIST=c(5,0))
out <- Arney.CSI(tempT, CSIflag=TRUE)
round(out$AO[1],0) == round(55.84,0)
round(out$CSI[1],0) == 35

#R1 <- crown.width.Arney(DBH=tempT$DBH[2])/2  # larget tree
#R2 <- crown.width.Arney(DBH=tempT$DBH[1])/2  # smaller tree
#circle.circle.area.intersection(Xa=0, Ya=0, Ra=R1, Xb=tempT$DIST[1], Yb=0, Rb=R2, plot.it=TRUE)  # 55.84

### Test 4 -- DIST < x1 (Equation 3 Appendix 1 page 5) (CHECK)  

tempT <- data.frame(DBH=c(2,4), DIST=c(4,0))
out <- Arney.CSI(tempT, CSIflag=TRUE)
round(out$AO[1],0) == round(39.2266,0)
round(out$CSI[1],0) == 68

#R1 <- crown.width.Arney(DBH=tempT$DBH[2])/2  # larget tree
#R2 <- crown.width.Arney(DBH=tempT$DBH[1])/2  # smaller tree
#circle.circle.area.intersection(Xa=0, Ya=0, Ra=R1, Xb=tempT$DIST[1], Yb=0, Rb=R2, plot.it=TRUE)  # 39.2266

### Test 5 --  competitor completely within large (Equation 4 Appendix 1 page 5) (CHECK)

tempT <- data.frame(DBH=c(2,7), DIST=c(4,0))
out <- Arney.CSI(tempT, CSIflag=TRUE)
round(out$AO[1],0) == round(52.55,0)
round(out$CSI[1],0) == 21

#100*52.55/(pi*R1^2)

#R1 <- crown.width.Arney(DBH=tempT$DBH[2])/2  # larget tree
#R2 <- crown.width.Arney(DBH=tempT$DBH[1])/2  # smaller tree
#circle.circle.area.intersection(Xa=0, Ya=0, Ra=R1, Xb=tempT$DIST[1], Yb=0, Rb=R2, plot.it=TRUE)  # 52.55

### Test 6 -- Table Values

# Type 3 DBHsubject > DBHcompetitor
tempT <- data.frame(DBH=c(5,10), DIST=c(7,0))
Arney.CSI(tempT, CSIflag=TRUE) # AO=130, CSI=31 Type=3 (CHECK)

# Type 3 DBHsubject < DBHcompetitor
tempT <- data.frame(DBH=c(10,5), DIST=c(7,0))
Arney.CSI(tempT, CSIflag=TRUE) # AO=130, CSI=84 Type=3 (CHECK) 

# Type 3 DBHsubject > DBHcompetitor
tempT <- data.frame(DBH=c(4,12), DIST=c(8,0))
Arney.CSI(tempT, CSIflag=TRUE) # AO=111, CSI=20 Type=3 (CHECK) 

# Type 3 DBHsubject < DBHcompetitor
tempT <- data.frame(DBH=c(12,4), DIST=c(8,0))
Arney.CSI(tempT, CSIflag=TRUE) # AO=111, CSI=96 Type=3 (CHECK) 

# Type 1 DBHsubject < DBHcompetitor
tempT <- data.frame(DBH=c(4,12), DIST=c(12,0))
Arney.CSI(tempT, CSIflag=TRUE) # AO=66, CSI=12 Type=1 (CHECK) 

# Type 1 DBHsubject < DBHcompetitor
tempT <- data.frame(DBH=c(12,4), DIST=c(12,0))
Arney.CSI(tempT, CSIflag=TRUE) # AO=66, CSI=57 Type=1 (CHECK) 

# Type 0 DBHsubject > DBHcompetitor
tempT <- data.frame(DBH=c(6,8), DIST=c(20,0))
Arney.CSI(tempT, CSIflag=TRUE) # AO=0, CSI=0 Type=0 (CHECK) 

# Type 0 DBHsubject < DBHcompetitor
tempT <- data.frame(DBH=c(8,6), DIST=c(20,0))
Arney.CSI(tempT, CSIflag=TRUE) # AO=0, CSI=0 Type=0 (CHECK) 

# tink there is a problem with equation 2 which I think is not fixed

# Type 2 DBHsubject > DBHcompetitor
# note had to add a round of DIST and x1 is code to check this for eq 2
tempT <- data.frame(DBH=c(3.91712,8.228396), DIST=c(7.999998,0))
Arney.CSI(tempT, CSIflag=TRUE) # AO=73, CSI=24 Type=2 (CHECK)

#source("C:/Users/DavidFolder/WorkC/R_Packages/myGeometryR/R/circle_circle_area_intersection.R")
#R1 <- crown.width.Arney(DBH=tempT$DBH[2])/2  # larget tree
#R2 <- crown.width.Arney(DBH=tempT$DBH[1])/2  # smaller tree
#circle.circle.area.intersection(Xa=0, Ya=0, Ra=R1, Xb=tempT$DIST[1], Yb=0, Rb=R2, plot.it=TRUE)  # 72.90

#0.5*pi*r2^2 + r1^2*asin(r2/r1) - DIST*r2 # 72.9 AO
#100*72.90/(pi*R1^2) # 23.2 CSI  

# Type 2 DBHsubject < DBHcompetitor
# note had to add a round of DIST and x1 is code to check this for eq 2
tempT <- data.frame(DBH=c(8.228396, 3.91712), DIST=c(7.999998,0))
Arney.CSI(tempT, CSIflag=TRUE) # AO=73, CSI=64 Type=2 (CSI=62 table?)
# the table value is a bit smaller but might be rounding (tables) and the DBHs I used

R1 <- crown.width.Arney(DBH=tempT$DBH[2])/2  # larget tree
R2 <- crown.width.Arney(DBH=tempT$DBH[1])/2  # smaller tree
#circle.circle.area.intersection(Xa=0, Ya=0, Ra=R1, Xb=tempT$DIST[1], Yb=0, Rb=R2, plot.it=TRUE)  # 72.90

#0.5*pi*r2^2 + r1^2*asin(r2/r1) - DIST*r2 # 72.9 AO
#100*72.90/(pi*R1^2) # 64.5 CSI

### Test 7 -- Appendix III tree-list

# tree-list (tree 14 is subject tree)
tempT <- data.frame(tree=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),
                   DBH=c(4,4,3,3,3,5,5,4,5,4,3,6,4,4),
                   DIST=c(7,13,10,4,7,12,10,12,7,13,5,9,13,0))

tempTx <- Arney.CSI(tempT, CSIflag=TRUE)

#  CSI <- 100*aij/Aj
# note the +100 is the subject tree (100*aj/Aj) = 100*1 = 100)                  # <- add this to function?
sum(trunc(tempTx$CSI)) + 100    # 328
# This is 1 unit off because tree 8 truncates to 0 and not 1
# The 100 is 100*AOsubject/AOsubject
sum(tempTx$CSI) + 100           # 332.4

### Test 8 -- Appendix II VP sample (single point) average CSI

BAF <- 20 # ft2/acre
NPTS <- 1

tempT <- data.frame(tree=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
                   DBH=c(11,11,12,12,12,13,14,15,15,15,16,16,17,17,18),
                   DIST=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))

tempT$EXPAN <- (BAF/(0.005454154*tempT$DBH^2))/NPTS

NTrees <- length(tempT$DBH); NTrees                          # 15
TPA <- sum(tempT$EXPAN); TPA                                 # 291.64 (CHECK)
BAPA <- sum(tempT$EXPAN * 0.005454154 * tempT$DBH^2); BAPA   # 300 ft2/acre (CHECK)

tempTx <- Arney.CSI(tempT, CSIflag=FALSE)

# average CSI
100 + sum(100*((tempTx$CA/43560)*tempTx$EXPAN))      # 537.5 vs 534 rounding?
# I think the first "100 + -" is the subject tree (100%)

#-------------------------------------------------------------------------------

