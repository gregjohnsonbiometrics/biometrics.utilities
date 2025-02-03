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
#' 
#' - Arney, J.D. 1973. Tables for quantifying competitive stress on individual trees. Pacific Forest Research Centre, Canadian Forest Service, Victoria, BC. Information Report BC-X-78. 47p.
#' 
#' - Clark, P.J., & Evans, F. C. 1954. Distance to Nearest Neighbor as a Measure of Spatial Relationships in Populations. Ecology, 35(4), 445–453.
#' 
#' - Curtis, R.0. 1967. Height-diameter, and height-diameter-age equations for second growth Douglas-fir. For. Sci. 365-375.
#' 
#' - Curtis, R.O. 1982. A Simple Index of Stand Density for Douglas-fir. Forest Sci., Vol. 28, No. 1, pp. 92-94.#' - Ritchie, M.W. and D.W. Hann. 1990. Equations for predicting height growth of six conifer species in southwest Oregon. Oregon State University, Forest Research Laboratory, Corvallis, Oregon. Research Paper 54. 12p.
#'
#' - Donnelly, K. 1978. Simulations to determine the variance and edge-effect of total nearest neighbour distance. In I. Hodder (ed.) Simulation studies in archaeology, Cambridge/New York: Cambridge University Press, pp 91–95.
#' 
#' - Curtis, R.O.; Marshall, D.D. 2000, Why quadratic mean diameter? Western Journal of Applied Forestry, 15 (3): 137–139.
#' 
#' - Hann, D.W. 1998. Equations for predicting the largest crown width of stand-grown trees in western Oregon. Forest Research Laboratory, Oregon State University, Corvallis. Research Contribution 17. 14p.
#'  
#' - Hann, D.W. 1999. An adjustable predictor of crown profile for stand-grown Douglas-fir trees. For. Sci. 45: 217–225.
#' 
#' - Gerrard, D.J. 1969. Competition Quotient: a new measure of the competition affecting individual forest trees. Michigan State University, Agr. Exp. Sta. Res. Bull. 20. 32pp.
#' 
#' - Glover G.R. and Hool, J.N. 1979. A basal area ratio predictor of loblolly pine plantation mortality. For. Sci. 25:275-282.#' - Hegyi, F. 1974. A simulation model for managing jack-pine stands. J. Fries (Ed.), Growth Models for Tree and Stand Simulation, Royal College of Forestry, Stockholm, Sweden (1974), pp. 74-90..
#' 
#' - Krajicek, J.E., K.A. Brinkman, and F.S. Gingrich. 1961. Crown competition: a measure of density. For. Sci. 7:36 – 42.
#' 
#' - Marshall D.D, G.P. Johnson, and D.W. Hann. 2003. Crown profile equations for stand-grown western hemlock trees in northwestern Oregon. Can. J. For. Res. 33: 2059–2066. 
#' 
#' - Paine, D.P., and D.W. Hann. 1982. Maximum crown width equations for southwestern Oregon tree species.  Forest Research Laboratory, Oregon State University, Corvallis. Research Bulletin 51. 9p.
#' 
#' - Ripley, B.D. 1977. Modelling spatial patterns (with discussion). Journal of the Royal Statistical Society, Series B, 39, 172 -- 212.
#' 
#' - Ritchie, M.W. and D.W. Hann. 1990. Equations for predicting height growth of six conifer species in southwest Oregon. Oregon State University, Forest Research Laboratory, Corvallis, Oregon. Research Paper 54. 12p.
#' 
#' - Reineke, L.H. 1933. Perfecting a stand density index for even-aged forests. Jour. Agric. Res. 46: 627 – 638.
#' 
#' - Wilson, F.G. 1946. Numerical expression of stocking in terms of height. Journal of Forestry, 44:758–761.


"_PACKAGE"

#' @importFrom Rcpp evalCpp
#' @useDynLib biometric.utilities
