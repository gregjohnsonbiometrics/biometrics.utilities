#' @title biometrics.utilities
#' 
#' @description
#' A collection of biometric metrics written in C++ for computational efficiency. The metrics focus on competition measures at 
#' both the tree and stand or plot level.
#' 
#' @author Greg Johnson
#' 
#' @references
#' 
#' - Clark, P. J., & Evans, F. C. 1954. Distance to Nearest Neighbor as a Measure of Spatial Relationships in Populations. Ecology, 35(4), 445–453.
#' 
#' - Curtis, R.O. 1982. A Simple Index of Stand Density for Douglas-fir. Forest Sci., Vol. 28, No. 1, pp. 92-94.#' - Ritchie, M.W. and D.W. Hann. 1990. Equations for predicting height growth of six conifer species in southwest Oregon. Oregon State University, Forest Research Laboratory, Corvallis, Oregon. Research Paper 54. 12p.
#' 
#' - Curtis, Robert O.; Marshall, David D. 2000, Why quadratic mean diameter? Western Journal of Applied Forestry, 15 (3): 137–139.
#' 
#' - Hann, D.W. 1998. Equations for predicting the largest crown width of stand-grown trees in western Oregon. Forest Research Laboratory, Oregon State University, Corvallis. Research Contribution 17. 14p.
#'  
#' - Hann, D.W. 1999. An adjustable predictor of crown profile for stand-grown Douglas-fir trees. For. Sci. 45: 217–225.
#' 
#' - Hegyi, F. 1974. A simulation model for managing jack-pine stands. J. Fries (Ed.), Growth Models for Tree and Stand Simulation, Royal College of Forestry, Stockholm, Sweden (1974), pp. 74-90..
#' 
#' - Krajicek, J.E., K.A. Brinkman, and F.S. Gingrich. 1961. Crown competition: a measure of density. For. Sci. 7:36 – 42.
#' 
#' - Marshall D.D, G.P. Johnson, and D.W. Hann. 2003. Crown profile equations for stand-grown western hemlock trees in northwestern Oregon. Can. J. For. Res. 33: 2059–2066. 
#' 
#' - Paine, D.P., and D.W. Hann. 1982. Maximum crown width equations for southwestern Oregon tree species.  Forest Research Laboratory, Oregon State University, Corvallis. Research Bulletin 51. 9p.
#' 
#' - Reineke, L.H. 1933. Perfecting a stand density index for even-aged forests. Jour. Agric. Res. 46: 627 – 638.
#' 
#' - Wilson, F.G. 1946. Numerical expression of stocking in terms of height. Journal of Forestry, 44:758–761.
#'      
#' 
"_PACKAGE"

#' @importFrom Rcpp evalCpp
#' @useDynLib biometric.utilities
