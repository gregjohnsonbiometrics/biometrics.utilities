#include <Rcpp.h>
#include "utilities.hpp"

//' @title bal() compute basal area in larger trees.
//' @name bal
//'
//' @param dbh            : double | vector of diameter at breast height
//' @param expansion      : double | vector of expansion factors 
//' @param imperial_units : bool   | TRUE = imperial (default), FALSE = metric
//'
//' @description
//' Sorts the dbh and expansion vectors in decreasing order of dbh and computes cumulative basal area in larger trees. Ties in dbh are
//' handled appropriately (each tree has the same bal).
//'
//' @return
//' Returns numeric vector with the basal area in larger trees in the original tree order (unsorted):
//'
//' @examples
//' data(treelist)
//' bal( treelist$dbh, treelist$tpa )
//'
//' @export
// [[Rcpp::export]]

std::vector<double> bal( const std::vector<double> dbh, 
                         const std::vector<double> expansion, 
                         bool imperial_units = true )
{
    auto bal_vector = compute_bal( dbh, expansion, imperial_units );

    return bal_vector;
}


//' @title ccfl() compute crown competition factor in larger trees.
//' @name ccfl
//'
//' @param dbh             : double | vector of diameter at breast height
//' @param mcw             : double | vector of maximum crown width for the tree record
//' @param expansion       : double | vector of expansion factors 
//' @param imperial_units  : bool   | TRUE = imperial, FALSE = metric
//'
//' @description
//' Sorts the dbh, mcw, and expansion vectors in decreasing order of dbh and computes cumulative maximum crown area in larger trees. Ties in dbh are
//' handled appropriately (each tree has the same ccfl).
//' 
//' The algorithm computes the percentage of an acre or hectare (depending on \code{imperial_units}) covered by each tree
//' record's crown.  It uses the maximum crown width to estimate the coverage for each tree and expands it to a per
//' acre or hectare basis using the supplied expansion factor.
//'
//' @return
//' Returns numeric vector with the crown competition factor in larger trees in the original tree order (unsorted):
//'
//' @examples
//' library( dplyr )
//' data(treelist)
//' # use Hann's (1998) maximum crown width equation for Douglas-fir, western hemlock, bigleaf maple, and red alder
//' mcw.parms <- data.frame( species=c(202,263,312, 351),
//'                          c0=c(4.6198,4.5652,4.0953,8.0000),
//'                          c1=c(1.8426,1.4147,2.3849,1.5300),
//'                          c2=c(-0.011311,0.000000,-0.011630,0.00000))
//' temp <- treelist %>% inner_join( mcw.parms, by="species" ) %>%
//'             mutate( mcw=c0 + c1*dbh + c2*dbh^2)
//' ccfl( temp$dbh, temp$mcw, temp$tpa )
//'
//' @export
// [[Rcpp::export]]

std::vector<double> ccfl( const std::vector<double> dbh, 
                          const std::vector<double> mcw, 
                          const std::vector<double> expansion, 
                          bool imperial_units = true )
{
    auto ccfl_vector = compute_ccfl( dbh, mcw, expansion, imperial_units );

    return ccfl_vector;
}

//' @title cch() compute closure at tree tip without using interpolation.
//' @name cch
//' 
//' @param dbh             : double | vector of diameter at breast height
//' @param height          : double | vector of total height
//' @param crown_length    : double | vector of crown lengths
//' @param dacb            : double | vector of distance above crown base to largest width for each tree
//' @param lcw             : double | vector of largest crown width for each tree
//' @param expansion       : double | vector of expansion factors 
//' @param parameters      : double | vector of parameters (see below)
//' @param imperial_units  : bool   | TRUE = imperial (default), FALSE = metric
//'
//' @description 
//' Computes the crown area as a fraction of an acre or hectare (depending on \code{imperial_units}) that lies in a horizontal
//' plane tangential to the tip of the tree.
//' 
//' cch() depends on \code{dacb} and \code{lcw} which must be estimated for each tree independently (normally species and size dependent).
//' 
//' The \code{parameters} vector contains coefficients for an equation that adjusts \code{lcw} for the position within the crown of the form:
//' \eqn{adjustment = rp^{(\beta_0 + \beta_1 rp^{0.5} + \beta_2 height / dbh)}}
//' where \code{rp} is the relative position in the crown, and \code{adjustment} is the estimated percentage of the \code{lcw} at that point
//' in the crown. The \code{parameters} are the \eqn{\beta} values in the equation. A simple cone can be created by a \code{parameters} 
//' vector of c(1.0,0.0,0.0).
//' 
//' @return
//' A vector of \code{cch} values for each tree.
//' 
//' @examples
//' data(treelist)
//' library( dplyr )
//' # use Paine and Hann (1992) maximum crown width equations for Douglas-fir, western hemlock, bigleaf maple, and red alder
//' mcw.parms <- data.frame( species=c(202,263,312,351),
//'                          c0=c(4.6198,4.5652,4.0953,8.0000),
//'                          c1=c(1.8426,1.4147,2.3849,1.5300),
//'                          c2=c(-0.011311,0.000000,-0.011630,0.00000))
//'
//' # use Hann's (1998) largest crown width equations Douglas-fir, western hemlock, bigleaf maple, and red alder
//' lcw.parms <- data.frame( species=c(202,263,312,351), 
//'                          b0=c(0.0,         0.0, 0.0,      0.3227140 ),
//'                          b1=c(0.004363240, 0.0, 0.0,      0.0       ),
//'                          b2=c(0.6020020,   0.0, 1.470180, 0.0       ))
//'
//' # estimate distance above crown base: use Hann 1999 for Douglas-fir, Marshall et al. 2003 for western hemlock, and assuming 0.0 for 
//' # bigleaf maple and red alder
//' dacb.parms <- data.frame( species=c(202,263,312,351), 
//'                          d0=c(0.06200, 0.355270, 0.0, 0.0 ) )
//'
//' # use crown profile parameters for Douglas-fir
//' cch.parms <- c( 0.929973, -0.1352120, -0.0157579 )
//'
//' temp <- treelist %>% 
//'             inner_join( mcw.parms, by="species" ) %>%
//'             inner_join( lcw.parms, by="species" ) %>%
//'             inner_join( dacb.parms, by="species" ) %>%
//'             mutate( crown_length = height - htlc,
//'                     cr = crown_length/height,
//'                     mcw = c0 + c1*dbh + c2*dbh^2,
//'                     lcw = mcw * cr^(b0+b1*crown_length+b2*(dbh/height)),
//'                     dacb = d0 * crown_length )
//'                     
//' cch( temp$dbh, temp$height,  temp$crown_length, temp$dacb, temp$lcw, temp$tpa, cch.parms )
//'
//' @export
// [[Rcpp::export]]
std::vector<double> cch( const std::vector<double> dbh,
                         const std::vector<double> height, 
                         const std::vector<double> crown_length,
                         const std::vector<double> dacb,
                         const std::vector<double> lcw,
                         const std::vector<double> expansion,
                         const std::vector<double> parameters,
                         bool imperial_units = true )
{
    std::vector<double> cch_vector;

    for( auto &ht : height )
    {
        cch_vector.push_back( compute_cch( ht, dbh, height, crown_length, dacb, lcw, expansion, parameters, imperial_units ) );
    }

    return cch_vector;
}

//' @title dominant_height() compute the average height of trees in the dominant cohort.
//' @name dominant_height
//'
//' @param height               : double | vector of total height
//' @param dbh                  : double | vector of diameter at breast height
//' @param expansion            : double | vector of expansion factors 
//' @param dominant_cohort_size : int    | number of trees in the dominant height cohort
//' @param method               : int    | 0 (default), 1, or 2 (see below for definitions)
//'
//' @description
//' Compute the expansion factor weighted average height of trees in the defined dominant tree cohort. Each method
//' defines the cohort differently:
//' 
//' \itemize{
//'    \item 0 : average height of the dominant_cohort_size trees by decreasing dbh
//'    \item 1 : average height of the dominant_cohort_size trees by decreasing height
//'    \item 2 : Lorey height (height of the tree of average basal area)
//' }
//'
//' @return
//' Returns the specified dominant height.
//'
//' @examples
//' data(treelist)
//' # compute the height of the 40 largest trees by dbh
//' dominant_height( treelist$height, treelist$dbh, treelist$tpa, 40, 0 )
//' # compute the height of the 100 largest trees by height
//' dominant_height( treelist$height, treelist$dbh, treelist$tpa, 100, 1 )
//'
//' @export
// [[Rcpp::export]]

double dominant_height( const std::vector<double> height,
                        const std::vector<double> dbh,
                        const std::vector<double> expansion,
                        const int dominant_cohort_size,
                        const int method  = 0 )
{
    return compute_dominant_height( height, dbh, expansion, dominant_cohort_size, method );
}

//' @title qmd() compute the quadratic mean diameter of a stand or plot.
//' @name qmd
//'
//' @param dbh                  : double | vector of diameter at breast height
//' @param expansion            : double | vector of expansion factors 
//'
//' @description
//' Compute the quadratic mean diameter (basal area weighted mean diameter) of a stand or plot.
//'
//' @return
//' Returns the qmd.
//'
//' @examples
//' data(treelist)
//' # compute the quadratic mean diameter
//' qmd( treelist$dbh, treelist$tpa )
//'
//' @export
// [[Rcpp::export]]

double qmd( const std::vector<double> dbh,
            const std::vector<double> expansion )
{
    return compute_qmd( dbh, expansion );
}


//' @title relative_spacing() compute Wilson's Relative Spacing
//' @name relative_spacing
//'
//' @param expansion       : double | vector of expansion factors 
//' @param dominant_height : double | dominant height of the stand
//' @param imperial_units  : bool   | TRUE = imperial (default), FALSE = metric
//'
//' @description
//' Compute Wilson's Relative Spacing (average tree spacing relative to dominant height). The smaller
//' the relative spacing, the more crowded or overstocked the stand or plot is.
//' 
//' @return
//' Returns the relative spacing as a fraction of dominant height.
//'
//' @examples
//' data(treelist)
//' # compute the relative spacing for HT40 dominant height
//' dom.ht <- dominant_height( treelist$height, treelist$dbh, treelist$tpa, 40, 0 )
//' relative_spacing( treelist$tpa, dom.ht )
//'
//' @export
// [[Rcpp::export]]
double relative_spacing( const std::vector<double> expansion,
                         const double dominant_height,
                         bool imperial = true )
{
    return compute_relative_spacing( expansion, dominant_height, imperial );
}


//' @title curtis_rd() compute Curtis' Relative Density
//' @name curtis_rd
//'
//' @param dbh             : double | vector of diameters at breast height 
//' @param expansion       : double | vector of expansion factors 
//' @param imperial_units  : bool   | TRUE = imperial (default), FALSE = metric
//'
//' @description
//' Compute Curtis' Relative Density (Curtis 1982).
//'
//' The expression \eqn{RD = G/(Dg^{0.5})}, where G is basal area and Dg is quadratic mean
//' stand diameter, provides a simple and convenient scale of relative stand density for Douglas-fir, 
//' equivalent to other generally accepted diameter-based stand density measures. 
//' 
//' @return
//' Returns the Curtis' Relative Density.
//'
//' @examples
//' data(treelist)
//' curtis_rd( treelist$dbh, treelist$tpa )
//'
//' @export
// [[Rcpp::export]]
double curtis_rd( const std::vector<double> dbh,
                  const std::vector<double> expansion,
                  bool imperial_units = true )
{
    return compute_curtis_rd( dbh, expansion, imperial_units );
}                            


//' @title reineke_sdi() compute Reineke's Stand Density Index
//' @name reineke_sdi
//'
//' @param dbh             : double | vector of diameters at breast height 
//' @param expansion       : double | vector of expansion factors 
//' @param imperial_units  : bool   | TRUE = imperial (default), FALSE = metric
//'
//' @description
//' Compute Reineke's Stand Density Index (Reineke 1933).
//'
//' Reineke (1933) developed a stand density index (SDI) that relates the current
//' stand density to an equivalent density in a stand with a quadratic mean diameter (Dq) of 10 inches.
//' Reineke’s SDI can be expressed as:
//' \eqn{SDI = N(Dq/10)^b}
//' where SDI = Reineke’s stand-density index, N = trees per acre, Dq = quadratic mean diameter (inches),
//' b = exponent of Reineke’s equation, often reported to equal –1.605.
//'
//' The function allows for computing the equivalent index in metric units.
//' 
//' @return
//' Returns the Reineke's Stand Density Index.
//'
//' @examples
//' data(treelist)
//' reineke_sdi( treelist$dbh, treelist$tpa )
//'
//' @export
// [[Rcpp::export]]
double reineke_sdi( const std::vector<double> dbh,
                    const std::vector<double> expansion,
                    bool imperial_units = true )
{
    return compute_reineke_sdi( dbh, expansion, imperial_units );
}  


//' @title ccf() compute crown competition factor
//' @name ccf
//'
//' @param crown_width     : double | vector of open-grown crown widths 
//' @param expansion       : double | vector of expansion factors 
//' @param imperial_units  : bool   | TRUE = imperial (default), FALSE = metric
//'
//' @description
//' Compute Crown Competition Factor (Krajicek, et al. 1961). 
//' 
//' Crown competition factor is the ratio of the open-grown crown area of all trees as a percentage of an acre:
//' \eqn{CCF = 100 \sum{CA expansion} /43560}, 
//' where CA is the open-grown crown area of a given tree.
//'
//' The function allows for computing the equivalent index in metric units.
//' 
//' @return
//' Returns the crown competition factor.
//'
//' @examples
//' data(treelist)
//' library( dplyr )
//' # use Paine and Hann (1992) maximum crown width equations for Douglas-fir, western hemlock, bigleaf maple, and red alder
//' mcw.parms <- data.frame( species=c(202,263,312,351),
//'                          c0=c(4.6198,4.5652,4.0953,8.0000),
//'                          c1=c(1.8426,1.4147,2.3849,1.5300),
//'                          c2=c(-0.011311,0.000000,-0.011630,0.00000))
//'
//' # use Hann's (1998) largest crown width equations Douglas-fir, western hemlock, bigleaf maple, and red alder
//' lcw.parms <- data.frame( species=c(202,263,312,351), 
//'                          b0=c(0.0,         0.0, 0.0,      0.3227140 ),
//'                          b1=c(0.004363240, 0.0, 0.0,      0.0       ),
//'                          b2=c(0.6020020,   0.0, 1.470180, 0.0       ))
//'
//' temp <- treelist %>% 
//'             inner_join( mcw.parms, by="species" ) %>%
//'             inner_join( lcw.parms, by="species" ) %>%
//'             mutate( crown_length = height - htlc,
//'                     cr = crown_length/height,
//'                     mcw = c0 + c1*dbh + c2*dbh^2,
//'                     lcw = mcw * cr^(b0+b1*crown_length+b2*(dbh/height)) )
//'
//' ccf( temp$lcw, temp$tpa )
//'
//' @export
// [[Rcpp::export]]
double ccf( const std::vector<double> crown_width,
            const std::vector<double> expansion,
            bool imperial_units = true )
{
    return compute_ccf( crown_width, expansion, imperial_units );
}