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
#' - Ritchie, M.W. and D.W. Hann. 1990. Equations for predicting height growth of six conifer species in southwest Oregon. Oregon State University, Forest Research Laboratory, Corvallis, Oregon. Research Paper 54. 12p.
#' 
#' - Hann, D.W. 1999. An adjustable predictor of crown profile for stand-grown Douglas-fir trees. For. Sci. 45: 217–225.
#' 
#' - Krajicek, J.E., K.A. Brinkman, and F.S. Gingrich. 1961. Crown competition: a measure of density. For. Sci. 7:36 – 42.
#' 
#' - Curtis, R.O. 1982. A Simple Index of Stand Density for Douglas-fir. Forest Sci., Vol. 28, No. 1, pp. 92-94.
#' 
#' - Reineke, L.H. 1933. Perfecting a stand density index for even-aged forests. Jour. Agric. Res. 46: 627 – 638.
#' 
#' - Wilson, F.G. 1946. Numerical expression of stocking in terms of height. Journal of Forestry, 44:758–761.
#' 
#' - Curtis, Robert O.; Marshall, David D. 2000, Why quadratic mean diameter? Western Journal of Applied Forestry, 15 (3): 137–139.
#'      
#' 
"_PACKAGE"

#' @importFrom Rcpp evalCpp
#' @useDynLib biometric.utilities
