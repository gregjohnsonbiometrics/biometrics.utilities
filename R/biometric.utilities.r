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
#'  - Wilson, F.G. 1946. Numerical expression of stocking in terms of height. Journal of Forestry, 44:758–761.  
#'  
#'  - Curtis, R.O. 1982. A simple index of stand density for Douglas-fir. Forest Science, 28:92–94. 
#'  
#'  - Reineke, L. H. 1933, Perfecting a stand-density index for even-aged forest. Journal of Agricultural Research, 46: 627-638. 
#'  
#'  - Krajicek, J.E., K.A. Brinkman, and F.S. Gingrich. 1961. Crown competition: a measure of density. Forest Science, 7:36–42.   
#'      
#' 
"_PACKAGE"

#' @importFrom Rcpp evalCpp
#' @useDynLib biometric.utilities
