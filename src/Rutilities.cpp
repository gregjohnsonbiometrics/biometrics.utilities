#include <Rcpp.h>
#include "utilities.hpp"
#include <unordered_map>
#include <set>

//' @title ba() compute basal area per unit area.
//' @name ba
//'
//' @param dbh            : double | vector of diameter at breast height
//' @param expansion      : double | vector of expansion factors 
//' @param imperial_units : bool   | TRUE = imperial (default), FALSE = metric
//'
//' @description
//' Computes the basal area per unit area of the trees in the tree list. 
//'
//' \eqn{ba = \sum{dbh_i^2 expansion_i k}}
//'
//' where k converts squared diameters to square feet or meters depending on \code{imperial_units}.
//'
//' @return
//' Returns the basal area per unit area or NAN if there are improper arguments.
//'
//' @examples
//' data(treelist)
//' ba( treelist$dbh, treelist$tpa )
//'
//' @export
// [[Rcpp::export]]

double ba( const std::vector<double> dbh, 
    const std::vector<double> expansion, 
    const bool imperial_units = true )
{
    auto ba = compute_ba( dbh, expansion, imperial_units );

    return ba;
}

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
//' Returns numeric vector with the basal area in larger trees in the original tree order (unsorted).
//'
//' @examples
//' data(treelist)
//' bal( treelist$dbh, treelist$tpa )
//'
//' @export
// [[Rcpp::export]]

std::vector<double> bal( const std::vector<double> dbh, 
                         const std::vector<double> expansion, 
                         const bool imperial_units = true )
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
//' ccfl( treelist$dbh, mcw( treelist$species, treelist$dbh ), treelist$tpa )
//'
//' @export
// [[Rcpp::export]]

std::vector<double> ccfl( const std::vector<double> dbh, 
                          const std::vector<double> mcw, 
                          const std::vector<double> expansion, 
                          const bool imperial_units = true )
{
    auto ccfl_vector = compute_ccfl( dbh, mcw, expansion, imperial_units );

    return ccfl_vector;
}

//' @title cch() compute closure at tree tip without using interpolation.
//' @name cch
//' 
//' @param species         : int    | vector of species codes
//' @param dbh             : double | vector of diameter at breast height
//' @param height          : double | vector of total height
//' @param crown_length    : double | vector of crown lengths
//' @param dacb            : double | vector of distance above crown base to largest width for each tree
//' @param lcw             : double | vector of largest crown width for each tree
//' @param expansion       : double | vector of expansion factors 
//' @param imperial_units  : bool   | TRUE = imperial (default), FALSE = metric
//'
//' @description 
//' Computes the crown area as a fraction of an acre or hectare (depending on \code{imperial_units}) that lies in a horizontal
//' plane tangential to the tip of the tree.
//' 
//' cch() depends on \code{dacb} and \code{lcw} which must be estimated for each tree independently (normally species and size dependent).
//'
//' @return
//' A vector of \code{cch} values for each tree.
//' 
//' @examples
//' data(treelist)
//' library( dplyr )
//'
//' # estimate distance above crown base: use Hann 1999 for Douglas-fir, Marshall et al. 2003 for western hemlock, and assuming 0.0 for 
//' # bigleaf maple and red alder
//' dacb.parms <- data.frame( species=c(202,263,312,351), 
//'                          d0=c(0.06200, 0.355270, 0.0, 0.0 ) )
//'
//' temp <- treelist %>% mutate( crown_length = height - htlc,
//'                              cr = crown_length/height, 
//'                              mcw.hat=mcw( species, dbh ), 
//'                              lcw.hat = lcw( species, dbh, cr ) ) %>%
//'             inner_join( dacb.parms, by="species" ) %>%
//'             mutate( dacb = d0 * crown_length )
//'                     
//' cch( temp$species, temp$dbh, temp$height, temp$crown_length, temp$dacb, temp$lcw.hat, temp$tpa )
//'
//' @export
// [[Rcpp::export]]
std::vector<double> cch( const std::vector<int>    &species,
                         const std::vector<double> &dbh,
                         const std::vector<double> &height, 
                         const std::vector<double> &crown_length,
                         const std::vector<double> &dacb,
                         const std::vector<double> &lcw,
                         const std::vector<double> &expansion,
                         const bool imperial_units = true )
{
    std::vector<double> cch_vector;

    try {
        for( auto &ht : height )
            cch_vector.push_back( compute_cch( species, ht, dbh, height, crown_length, dacb, lcw, expansion, imperial_units ) );
    } catch( ... ) {
        throw;
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
//' The cohort size is ignored and should be 0.0 for Lorey height (option 2).
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
//' # compute the Lorey height 
//' dominant_height( treelist$height, treelist$dbh, treelist$tpa, 0, 2 )
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

//' @title dominant_dbh() compute the average dbh of trees in the dominant cohort.
//' @name dominant_dbh
//'
//' @param height               : double | vector of total height
//' @param dbh                  : double | vector of diameter at breast height
//' @param expansion            : double | vector of expansion factors 
//' @param dominant_cohort_size : int    | number of trees in the dominant height cohort
//' @param method               : int    | 0 (default), 1, or 2 (see below for definitions)
//'
//' @description
//' Compute the expansion factor weighted average dbh of trees in the defined dominant tree cohort. Each method
//' defines the cohort differently:
//' 
//' \itemize{
//'    \item 0 : average dbh of the dominant_cohort_size trees by decreasing dbh
//'    \item 1 : average dbh of the dominant_cohort_size trees by decreasing height
//'    \item 2 : Lorey dbh (quadratic mean diameter)
//' }
//'
//' The cohort size is ignored and should be 0.0 for Lorey dbh (option 2).
//'
//' @return
//' Returns the specified dominant dbh.
//'
//' @examples
//' data(treelist)
//' # compute the dbh of the 40 largest trees by dbh
//' dominant_dbh( treelist$height, treelist$dbh, treelist$tpa, 40, 0 )
//' # compute the dbh of the 100 largest trees by height
//' dominant_dbh( treelist$height, treelist$dbh, treelist$tpa, 100, 1 )
//' # compute the Lorey dbh 
//' dominant_dbh( treelist$height, treelist$dbh, treelist$tpa, 0, 2 )
//'
//' @export
// [[Rcpp::export]]

double dominant_dbh( const std::vector<double> height,
                     const std::vector<double> dbh,
                     const std::vector<double> expansion,
                     const int dominant_cohort_size,
                     const int method  = 0 )
{
    return compute_dominant_dbh( height, dbh, expansion, dominant_cohort_size, method );
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
                         const bool imperial = true )
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
                  const bool imperial_units = true )
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
                    const bool imperial_units = true )
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
//'
//' # use Hann's (1998) largest crown width equations Douglas-fir, western hemlock, bigleaf maple, and red alder
//' lcw.parms <- data.frame( species=c(202,263,312,351), 
//'                          b0=c(0.0,         0.0, 0.0,      0.3227140 ),
//'                          b1=c(0.004363240, 0.0, 0.0,      0.0       ),
//'                          b2=c(0.6020020,   0.0, 1.470180, 0.0       ))
//'
//' temp <- treelist %>% mutate( mcw.hat=mcw( species, dbh ) ) %>%
//'             inner_join( lcw.parms, by="species" ) %>%
//'             mutate( crown_length = height - htlc,
//'                     cr = crown_length/height,
//'                     lcw = mcw.hat * cr^(b0+b1*crown_length+b2*(dbh/height)) )
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


//' @title Clark_Evans_R() compute Clark Evan's R with Donnelly Edge Correction in a Polygonal Plot
//' @name Clark_Evans_R
//'
//' @param x         : double | vector of x coordinates of trees on plot
//' @param y         : double | vector of y coordinates of trees on plot
//' @param plotarea  : double | plot area if polygon not available
//' @param poly_x    : double | vector of plot polygon x coordinates. 
//' @param poly_y    : double | vector of plot polygon y coordinates. 
//'
//' @description
//' Compute the Clark and Evans Aggregation Index (R) (1954). The aggregation index R is a measure
//' of clustering or ordering of trees on a plot. It is the ratio of the observed mean nearest neighbor
//' distance in the trees to that expected for a Poisson point process of the same intensity. A value 
//' R > 1 suggests ordering, while R < 1 suggests clustering (unequal inter-tree competition). R has been
//' proposed as a two-sided, distance-dependent tree competition metric.
//'
//' This implementation uses Donnelly's correction (Donnelly 1978) for polyonal plots if a polygon is supplied.
//' If not polygon coordinates are not supplied, the function computes `R` with no edge correction and relies on
//' the supplied `plotarea` value.
//'
//' \eqn{R = \frac{\frac{\sum{ d_i }}{N}}{(\frac{A}{N})^{0.5}/2}}
//'
//' where: \eqn{d_i} is the nearest neighbor distance for the ith tree, A is the plot area, and N is the number of 
//' trees on the plot.
//'
//' @return
//' Returns the Clark Evans R statistic.
//'
//' @examples
//' data(treelistxy)
//' min_x <- min(treelistxy$x)
//' min_y <- min(treelistxy$y)
//' max_x <- max(treelistxy$x)
//' max_y <- max(treelistxy$y)
//' poly_x <- c(min_x, max_x, max_x, min_x)
//' poly_y <- c(min_y, min_y, max_y, max_y)
//' Clark_Evans_R( treelistxy$x, treelistxy$y, 0.0, poly_x, poly_y )
//'
//' @export
// [[Rcpp::export]]

double Clark_Evans_R( const std::vector<double> &x, 
                      const std::vector<double> &y,
                      double plotarea,
                      const std::vector<double> &poly_x,
                      const std::vector<double> &poly_y )
{
    return compute_R( x, y, plotarea, poly_x, poly_y );
}


//' @title Clark_Evans_R_circle() compute Clark Evan's R with Donnelly Edge Correction in a circular plot
//' @name Clark_Evans_R_circle
//'
//' @param x              : double | vector of x coordinates of trees on plot
//' @param y              : double | vector of y coordinates of trees on plot
//' @param plotarea       : double | plot area if polygon not available
//' @param plot_center_x  : double | x coordinate of plot center. 
//' @param plot_center_y  : double | y coordinate of plot center. 
//' @param plot_radius    : double | radius of circular plot. 
//'
//' @description
//' Compute the Clark and Evans Aggregation Index (R) (1954). The aggregation index R is a measure
//' of clustering or ordering of trees on a plot. It is the ratio of the observed mean nearest neighbor
//' distance in the trees to that expected for a Poisson point process of the same intensity. A value 
//' R > 1 suggests ordering, while R < 1 suggests clustering (unequal inter-tree competition). R has been
//' proposed as a two-sided, distance-dependent tree competition metric.
//'
//' This implementation uses Donnelly's correction (Donnelly 1978) for circular plots if a circle is supplied.
//' If not circle dimensions are not supplied, the function computes `R` with no edge correction and relies on
//' the supplied `plotarea` value.
//'
//' \eqn{R = \frac{\frac{\sum{ d_i }}{N}}{(\frac{A}{N})^{0.5}/2}}
//'
//' where: \eqn{d_i} is the nearest neighbor distance for the ith tree, A is the plot area, and N is the number of 
//' trees on the plot.
//'
//' @return
//' Returns the Clark Evans R statistic.
//'
//' @examples
//' data(treelistxy)
//' # TO DO
//'
//' @export
// [[Rcpp::export]]

double Clark_Evans_R_circle( const std::vector<double> &x, 
                             const std::vector<double> &y,
                             double plotarea,
                             const double plot_center_x,
                             const double plot_center_y, 
                             const double plot_radius )
{
    return compute_R( x, y, plotarea, Point(plot_center_x, plot_center_y), plot_radius );
}

//' @title Hegyi() compute Hegyi's distance weighted size ratio.
//' @name Hegyi
//'
//' @param x              : double | vector of x coordinates of trees on plot
//' @param y              : double | vector of y coordinates of trees on plot
//' @param dbh            : double | vector of diameter at breast height 
//' @param poly_x         : double | vector of plot polygon x coordinates (use 0 if no plot boundaries available). 
//' @param poly_y         : double | vector of plot polygon y coordinates (use 0 if no plot boundaries available). 
//' @param imperial_units : bool   | TRUE = imperial, FALSE = metric (default)
//'
//' @description
//' Compute Hegyi's (1974) distance-weighted size ratio (a two-sided competition index) for trees within a 6-meter fixed radius plot.
//' The ratio for the ith tree in the 6-meter radius plot is:
//'
//' \eqn{heygi_i = \sum{\frac{ba_j/ba_i}{d_{ij}}}}
//'
//' where \eqn{ba_j} and \eqn{ba_i} are the basal areas of the jth and ith tree respectively, \eqn{d_{ij}} is the distance between tree i and j.
//'
//' Trees from a plot of arbitrary size can be used. The Hegyi ratio for each tree will be computed based on its neighbors within the 6-meter boundary. 
//' If \code{imperial_units} is TRUE, the coordinates will be converted to meters prior to calculations.
//'
//' `Hegyi` adjusts for edge effects using Ripley's (1977) edge correction if plot boundaries are supplied, otherwise edge effects are ignored.
//'
//' @return
//' Returns the Hegyi's distance weighted size ratio for each tree. A `NaN` is returned if a tree has a `dbh` of 0.0.
//'
//' @examples
//' data(treelistxy)
//' min_x <- min(treelistxy$x)
//' min_y <- min(treelistxy$y)
//' max_x <- max(treelistxy$x)
//' max_y <- max(treelistxy$y)
//' poly_x <- c(min_x, max_x, max_x, min_x)
//' poly_y <- c(min_y, min_y, max_y, max_y)
//' Hegyi( treelistxy$x, treelistxy$y, treelistxy$dbh, poly_x, poly_y, imperial_units=T )
//'
//' @export
// [[Rcpp::export]]

std::vector<double> Hegyi( const std::vector<double> &x,
                           const std::vector<double> &y,
                           const std::vector<double> dbh,
                           const std::vector<double> &poly_x,
                           const std::vector<double> &poly_y,
                           const bool imperial_units = false )
{
    std::vector<Point> plot;

    if( poly_x.size() == 4 )
    {
        for( size_t i = 0; i < 4; ++i )
            plot.emplace_back( Point( poly_x[i], poly_y[i] ));
    }

    return compute_Hegyi( x, y, dbh, plot, imperial_units );
}


//' @title Arney_CSI() compute Arney's competitive stress index.
//' @name Arney_CSI
//'
//' @param x              : double | vector of x coordinates of trees on plot (in same units as mcw)
//' @param y              : double | vector of y coordinates of trees on plot (in same units as mcw)
//' @param dbh            : double | vector of diameter at breast height 
//' @param mcw            : double | vector of maximum crown width for the tree record
//' @param imperial_units : bool   | TRUE = imperial, FALSE = metric (default)
//'
//' @description
//' Compute Arney's Competitive Stress Index as described in Arney (1973). CSI is the sum of percentage competing trees crown area overlaping
//' a subject tree to the subject tree's crown area.
//'
//' \eqn{CSI_i = 100 \sum{\frac{AO_j}{CA_i}}}
//'
//' where \eqn{AO_j} is the area of overlap of tree j on subject tree i, \eqn{CA_i} is the crown area of the subject tree i, and \eqn{CSI_i} is the
//' competitive stress index for tree i.
//'
//' This version currently does not adjust for edge effects.
//'
//' @return
//' Returns the Arney's CSI for each tree.
//'
//' @examples
//' data(treelistxy)
//' Arney_CSI( treelistxy$x, treelistxy$y, treelistxy$dbh, mcw( treelistxy$fia, treelistxy$dbh ) )
//'
//' @export
// [[Rcpp::export]]

std::vector<double> Arney_CSI( const std::vector<double> &x,
                               const std::vector<double> &y,
                               const std::vector<double> &dbh,
                               const std::vector<double> &mcw )
{
    return compute_Arney_CSI( x, y, dbh, mcw );
}


//' @title mcw() Estimate maximum crown width for a species
//' @name mcw
//'
//' @param fia            : int    | vector of FIA species code
//' @param dbh            : double | vector of diameter at breast height
//' @param imperial_units : bool   | TRUE = imperial (default), FALSE = metric
//' @param allow_donor    : bool   | TRUE (default) allows donor species fits; FALSE restricts the resolution to the 220 direct fits. 
//' @param allow_none     : bool   | TRUE (default) returns them with a warning naming them; FALSE refuses the 11 no_donor codes, whose numbers are the clade population equation carried for want of any defensible model. 
//' @param hold_domain    : bool   | TRUE (default) restricts diameters to the domain; FALSE allows estimates for diameters outside the data domain. 
//' @param apply_ceiling  : bool   | TRUE (default) restricts `mcw` estimates to define maximums; FALSE accepts `mcw` estimates as-is.
//'
//' @description
//' Estimate the maximum (open-grown) crown width for trees using new FVS equations (Weiskittel, Johnson, Marshall 2026). 
//' If an FIA species code does not have a parameter set or a donor species, the conifer or hardwood default is used.
//' 
//' The equation form used to estimate mcw is:
//' 
//' \eqn{mcw = \beta_0 dbh^{\beta_1}}
//'
//' @return
//' Returns numeric vector with maximum crown width in the original tree order (unsorted).
//'
//' @examples
//' # compute maximum crown width for a 10 inch Douglas-fir, red alder, and Atlantic White-cedar (parameters not available, will default to conifer)
//' mcw( c(202,351,43), c(10,10,10) )
//'
//' @export
// [[Rcpp::export]]

std::vector<double> mcw( const std::vector<int> &fia,
                         const std::vector<double> &dbh,
                         const bool imperial_units = true,
                         const bool allow_donor = true,
                         const bool allow_none = true,
                         const bool hold_domain = true,
                         const bool apply_ceiling = true )
{
    auto mcw_vector = compute_mcw( fia, dbh, imperial_units, allow_donor, allow_none, hold_domain, apply_ceiling );

    return mcw_vector;
}


//' @title lcw() Estimate largest crown width (`lcw`), the stand grown crown width for a species
//' @name lcw
//'
//' @param fia            : int    | vector of FIA species code
//' @param dbh            : double | vector of diameter at breast height
//' @param cr             : double | crown ratio (fraction of total height)
//' @param imperial_units : bool   | TRUE = imperial (default), FALSE = metric
//' @param allow_donor    : bool   | TRUE (default) allows donor species fits; FALSE restricts the estimates to the 220 direct species estimates. 
//' @param allow_none     : bool   | TRUE (default) returns estimates for all species with a warning naming species with borrowed parameters; FALSE refuses the species with `no_donor` codes (the softwood/hardwood population equation is used).
//' @param hold_domain    : bool   | TRUE (default) restricts diameters to the fitting dataset domain; FALSE allows estimates for diameters outside the data domain. 
//' @param apply_ceiling  : bool   | TRUE (default) restricts `mcw` estimates to defined maximums; FALSE accepts `mcw` estimates as-is.
//'
//' @description
//' Estimates the largest stand grown crown width for trees using new FVS equations (Weiskittel, Marshall, Johnson 2026). 
//' If an FIA species code does not have a parameter set or a donor species, the conifer or hardwood default is used.
//' 
//' lcw, the stand grown crown width is not fit as its own allometry. It is a bounded fraction of the maximum crown width envelope,
//'
//' \eqn{lcw = mcw * ratio}
//' 
//'  with
//'  
//'  \eqn{logit(ratio) = \beta_0 + \beta_1 cr + \beta_2 dbh }
//'  
//'  ratio is a Beta regression at \eqn{\tau = 0.50} fit against the conifer/hardwood split mcw envelope. Since the ratio lies strictly in (0, 1), 
//'  lcw is below mcw for every species at every diameter.
//'
//' @return
//' Returns numeric vector with largest crown width in the original tree order (unsorted).
//'
//' @examples
//' # compute largest crown width for a 10 inch and crown ratio of 0.5 Douglas-fir, red alder, and Atlantic White-cedar 
//' lcw( c(202,351,43), c(10,10,10), c(0.5, 0.5, 0.5) )
//'
//' @export
// [[Rcpp::export]]
 
 std::vector<double> lcw( const std::vector<int> &fia,
                          const std::vector<double> &dbh,
                          const std::vector<double> &cr,
                          const bool imperial_units = true,
                          const bool allow_donor = true,
                          const bool allow_none = true,
                          const bool hold_domain = true,
                          const bool apply_ceiling = true )
 {
     auto lcw_vector = compute_lcw( fia, dbh, cr, imperial_units, allow_donor, allow_none, hold_domain, apply_ceiling );
     
     return lcw_vector;
 }


//' @title hd_fit() Fit height-dbh curve
//' @name hd_fit
//'
//' @param fia    : int    | vector of FIA species codes
//' @param dbh    : double | vector of diameter at breast height
//' @param height : double | vector of total heights
//' @param bh     : double | height to breast height (default is 4.5 feet; use 1.37 for metric measurements)
//'
//' @description
//' Fits a height-diameter curve to each species using provided \code{height} and \code{dbh} vectors. The functional form is:
//'
//' \eqn{\widehat{height} = BH + e^{(\beta_0 + \beta_1 dbh^{-\beta_2})}}
//'
//' where BH = height to breast height (e.g., 4.5 feet for imperial measurements, 1.37 meters for metric), and \eqn{\beta} s are
//' parameters to be estimated. Fitting height-dbh curves is often difficult due to limited measurement data (either in observation
//' count or in range, or both). We are using David Marshall's technique of fitting a linearized form of the equation while iterating
//' over a range of \eqn{\beta_2} values (-0.1 to -1.0). The \eqn{\beta_2} value yielding the lowest sum of squared errors (SSE) is chosen.
//'
//' If a species has less than 3 observations, no parameters are estimated and a vector of 0.0 is returned.
//'
//' @return
//' A \code{data.frame} for use in \code{hd_predict()} with the following members:
//' \itemize{
//'    \item fia     : FIA species code.
//'    \item beta_0  : \eqn{\beta_0} parameter estimate.
//'    \item beta_1  : \eqn{\beta_1} parameter estimate.
//'    \item beta_2  : \eqn{\beta_2} parameter estimate.
//' }
//'
//'
//' @examples
//' data(treelist )
//' hd.model <- hd_fit( treelist$species, treelist$dbh, treelist$height )
//' plot( treelist$dbh, hd_predict( hd.model, treelist$species, treelist$dbh ), col="green" )
//' points( treelist$dbh, treelist$height )
//' 
//' @export
// [[Rcpp::export]]

Rcpp::DataFrame hd_fit( const std::vector<int>    &fia,
                        const std::vector<double> &dbh,
                        const std::vector<double> &height,
                        const double bh = 4.5 )
{
    std::unordered_map<int,std::vector<double>> parms;

    // build set of fia species codes
    std::set<int> fia_codes;
    for( size_t i = 0; i < fia.size(); ++i )
        fia_codes.insert( fia[i] );

    for( auto &spp : fia_codes )
    {
        std::vector<double> t_dbh;
        std::vector<double> t_ht;
        for( size_t j = 0; j < dbh.size(); ++j )
        {
            // build species-specific vectors
            if( fia[j] == spp )
            {
                t_dbh.emplace_back( dbh[j] );
                t_ht.emplace_back( height[j] );
            }

            // fit regression
            if( t_dbh.size() >= 3 )
                parms[spp] = height_dbh_fit( t_dbh, t_ht, bh );
            else
                parms[spp] = {0.0, 0.0, 0.0};
        }
    }

     std::vector<int> fia_out( fia_codes.begin(), fia_codes.end() );
     size_t n = fia_out.size();
     std::vector<double> beta_0( n, 0.0 );
     std::vector<double> beta_1( n, 0.0 );
     std::vector<double> beta_2( n, 0.0 );
     for( size_t i = 0; i < n; ++i )
     {
        beta_0[i] = parms[fia_out[i]][0]; 
        beta_1[i] = parms[fia_out[i]][1]; 
        beta_2[i] = parms[fia_out[i]][2];  
     }

     return Rcpp::DataFrame::create (
        Rcpp::Named("fia")    = fia_out,
        Rcpp::Named("beta_0") = beta_0,
        Rcpp::Named("beta_1") = beta_1,
        Rcpp::Named("beta_2") = beta_2 );
}

//' @title hd_predict() Use parameter estimates from \code{hd_fit} to predict heights for a vector of \code{dbh}.
//' @name hd_predict
//'
//' @param parms : double | data.frame returned by \code{\link{hd_fit}}
//' @param fia   : int    | vector of FIA species codes
//' @param dbh   : double | vector of diameter at breast height
//' @param bh    : double | height to breast height (default is 4.5 feet; use 1.37 for metric measurements)
//'
//' @description
//' Predicts heights for a \code{dbh} vector. See \link{hd_fit} for details on the parameter estimates used in the prediction.
//'
//' @return
//' A vector of predicted heights for each tree supplied in the original order. Note: species without sufficient observations
//' will have a height of 0.0 returned.
//'
//' @examples
//' data(treelist )
//' hd.model <- hd_fit( treelist$species, treelist$dbh, treelist$height )
//' plot( treelist$dbh, hd_predict( hd.model, treelist$species, treelist$dbh ), col="green" )
//' points( treelist$dbh, treelist$height )
//' 
//' @export
// [[Rcpp::export]]

std::vector<double> hd_predict( Rcpp::DataFrame &hd_parameters,
                                const std::vector<int>    &fia,
                                const std::vector<double> &dbh,
                                const double bh = 4.5 )
{
    std::vector<double> ht_out;
    std::vector<int> id_out;
    std::vector<int> data_order( fia.size() );
    std::iota( data_order.begin(), data_order.end(), 0 );

    std::vector<int> spp = hd_parameters[0];
    std::vector<double> beta_0 = hd_parameters[1];
    std::vector<double> beta_1 = hd_parameters[2];
    std::vector<double> beta_2 = hd_parameters[3];

    for( size_t i = 0; i < spp.size(); ++i )
    {
        std::vector<double> t_dbh;
        for( size_t j = 0; j < dbh.size(); ++j )
        {
            // build species-specific vectors
            if( fia[j] == spp[i] )
            {
                t_dbh.emplace_back( dbh[j] );
                id_out.emplace_back( data_order[j] );
            }
        }

        auto t_ht = height_dbh_predict( {beta_0[i],beta_1[i],beta_2[i]}, t_dbh, bh );

        ht_out.insert(std::end(ht_out), std::begin(t_ht), std::end(t_ht));
    }

    std::vector<double> ht_out_sorted;
    for( auto i : sort_indices( id_out, true ) )
        ht_out_sorted.emplace_back( ht_out[i] );

    return ht_out_sorted;
}

//' @title Glover_Hool() compute the Glover and Hool competition index for each tree.
//' @name Glover_Hool
//'
//' @param dbh            : double | vector of diameter at breast height
//' @param expansion      : double | vector of expansion factors 
//' @param use_arithmetic : bool   | TRUE = use arithmetic mean (default), FALSE = use quadratic mean
//' @param imperial_units : bool   | TRUE = imperial (default), FALSE = metric
//'
//' @description
//' Compute the Glover and Hool (1979) competition index. The index is interpreted as the ratio
//' of a tree's basal area to the basal area of the tree of mean diameter. Glover and Hool used
//' the arithmetic mean and a common variation is to use the quadratic mean (use \code{use_arithmetic}
//' flag to select the desired method).
//'
//' The index \eqn{G_i} is:
//'
//' \eqn{G_i = dbh_i^2 / \overline{dbh}^2}
//'
//' where: \eqn{dbh_i} is the diameter of tree \code{i} and \eqn{\overline{dbh}} is the mean diameter (either arithmetic or quadratic).
//'
//' @return
//' Returns a vector of competition indicies for each tree in their original order.
//'
//' @examples
//' data(treelist)
//' # compute the Glover and Hool index with each mean
//' Glover_Hool( treelist$dbh, treelist$tpa, use_arithmetic=T )
//' Glover_Hool( treelist$dbh, treelist$tpa, use_arithmetic=F )
//'
//' @export
// [[Rcpp::export]]

std::vector<double> Glover_Hool( const std::vector<double> &dbh,
                                 const std::vector<double> &expansion,
                                 const bool use_arithmetic = true,
                                 const bool imperial_units = true )
{
    return compute_glover_hool( dbh, expansion, use_arithmetic, imperial_units );
}                                         


//' @title APA() compute Area Potentially Available (APA)
//' @name APA
//'
//' @param x        : double | vector of x coordinates for each tree
//' @param y        : double | vector of y coordinates for each tree
//' @param dbh      : double | vector of dbh for each tree
//' @param plot_x   : double | vector of plot corners' x coordinates (bottom left, top right)
//' @param plot_y   : double | vector of plot corners' y coordinates (bottom left, top right)
//' @param weighted : bool   | TRUE = use dbh^2 as weighting factor for polygon construction, FALSE = unweighted (default)
//'
//' @description
//' Computes Area Potentially Available (APA) using Brown's (1965) method. Polygons are constructed 
//' around the subject tree by the intersection of the perpendicular bisectors of the distance between 
//' the subject tree and competitors (creating a Voronoi tesselation of the plot).
//'
//' NOTE: the \code{weighted} option currently does not work.
//'
//' @return
//' Returns a vector of APA values (in square units of measure used for the coordinates) for each tree in their original order.
//'
//' @examples
//' data(treelistxy)
//' min_x <- min(treelistxy$x)
//' min_y <- min(treelistxy$y)
//' max_x <- max(treelistxy$x)
//' max_y <- max(treelistxy$y)
//' poly_x <- c(min_x, max_x)
//' poly_y <- c(min_y, max_y)
//' APA( treelistxy$x, treelistxy$y, treelistxy$dbh, poly_x, poly_y, F )
//'
//' @export
// [[Rcpp::export]]

std::vector<double> APA( const std::vector<double> &x,
                         const std::vector<double> &y,
                         const std::vector<double> &dbh,
                         const std::vector<double> &plot_x,
                         const std::vector<double> &plot_y,
                         const bool weighted=true )
{
    std::array<Point,2> plot_corners = {{ {plot_x[0],plot_y[0]},
                                          {plot_x[1],plot_y[1]} }};

    return compute_apa( x, y, dbh, plot_corners, weighted );
}                                         

//' @title APA_Polygons() compute Area Potentially Available (APA) polygons
//' @name APA_Polygons
//'
//' @param tree_id  : int    | vector of tree identification numbers
//' @param x        : double | vector of x coordinates for each tree
//' @param y        : double | vector of y coordinates for each tree
//' @param dbh      : double | vector of dbh for each tree
//' @param plot_x   : double | vector of plot corners' x coordinates (bottom left, top right)
//' @param plot_y   : double | vector of plot corners' y coordinates (bottom left, top right)
//' @param weighted : bool   | TRUE = use dbh^2 as weighting factor for polygon construction, FALSE = unweighted (default)
//'
//' @description
//' Computes Area Potentially Available (APA) using Brown's (1965) method. Polygons are constructed 
//' around the subject tree by the intersection of the perpendicular bisectors of the distance between 
//' the subject tree and competitors (creating a Voronoi tesselation of the plot).
//'
//' NOTE: the \code{weighted} option currently does not work.
//'
//' @return
//' Returns a \code{data.frame} of APA polygons with the following members:
//' \itemize{
//'    \item tree_id
//'    \item v_x : x coordinate of the tree's APA polygon
//'    \item v_y : y coordinate of the tree's APA polygon
//' }
//'
//' @examples
//' library( ggplot2 )
//' library( dplyr )
//' data(treelistxy)
//' min_x <- min(treelistxy$x)
//' min_y <- min(treelistxy$y)
//' max_x <- max(treelistxy$x)
//' max_y <- max(treelistxy$y)
//' poly_x <- c(min_x, max_x)
//' poly_y <- c(min_y, max_y)
//' p <- APA_Polygons( treelistxy$tree, treelistxy$x, treelistxy$y, treelistxy$dbh, poly_x, poly_y, F )
//' p %>% ggplot( aes( v_x, v_y, group=tree_id )) + geom_path()
//'
//' @export
// [[Rcpp::export]]

Rcpp::DataFrame APA_Polygons( 
                         const std::vector<int> &tree_id,
                         const std::vector<double> &x,
                         const std::vector<double> &y,
                         const std::vector<double> &dbh,
                         const std::vector<double> &plot_x,
                         const std::vector<double> &plot_y,
                         const bool weighted=false )
{
    std::array<Point,2> plot_corners = {{ {plot_x[0],plot_y[0]},
                                          {plot_x[1],plot_y[1]} }};

    std::vector<std::vector<Point>> polys = get_voronoi_polygons( x, y, dbh, plot_corners, weighted );

    size_t n = tree_id.size();
    std::vector<int> id;
    std::vector<double> v_poly_x;
    std::vector<double> v_poly_y;

    for( size_t i = 0; i < n; ++i )
    {
        for( size_t j = 0; j < polys[i].size(); ++j )
        {
            id.emplace_back( tree_id[i] );
            v_poly_x.emplace_back( polys[i][j].x );
            v_poly_y.emplace_back( polys[i][j].y );
        }
    }

    return Rcpp::DataFrame::create (
        Rcpp::Named("tree_id") = id,
        Rcpp::Named("v_x")  = v_poly_x,
        Rcpp::Named("v_y")  = v_poly_y );
}


//' @title crown_profile() Estimates the crown width for each tree for a specified height above ground.
//' @name crown_profile
//'
//' @param fia            : int    | vector of FIA species code
//' @param dbh            : double | vector of diameter at breast height
//' @param height         : double | vector of total height
//' @param hcb            : double | vector of height to crown base
//' @param lcw            : double | vector of largest crown width
//' @param h              : double | height compute crown width
//' @param imperial_units : bool   | TRUE = imperial (default), FALSE = metric
//'
//' @description
//' Estimates the crown width at a specified height above ground. The function uses the equation for relative crown shape:
//'  
//' \eqn{cw_h = lcw (1 - z^x)^(1/y)}
//' 
//' with z running from 0 at the widest point in the crown scaled by lcw to 1 at the tip, x and y are species-specific parameters, 
//' and a per species rule for where the widest point sits. Any specie without its own profile parameters (currently only Douglas-fir, 
//' western hemlock, and red alder have parameters) falls back to a softwood or hardwood default.
//'
//' @return
//' Returns a numeric vector with crown width at height h for each tree in the original tree order (unsorted). Zero is returned for values of h
//' that are below a tree's height to crown base.
//'
//' @examples
//' # compute crown width at 
//' fia <- c(202,251,43)
//' dbh <- c(10,10,10)
//' height <- c(60,60,60)
//' hcb <- c(20,20,20) 
//' cr <- (1-hcb)/height
//' lcw <- lcw( fia, dbh, cr )
//' h <- 40
//' 
//' crown_profile( fia, dbh, height, hcb, lcw, h )
//'
//' @export
// [[Rcpp::export]]
std::vector<double> crown_profile( 
                      const std::vector<int>    &fia, 
                      const std::vector<double> &dbh, 
                      const std::vector<double> &height,
                      const std::vector<double> &hcb, 
                      const std::vector<double> &lcw,
                      const double &h, 
                      const bool imperial_units = true )
{
    std::vector<double> cwa( fia.size(), 0.0 );
    
    for( size_t i = 0; i < fia.size(); ++i )
    {
        cwa[i] = compute_crown_profile( fia[i], dbh[i], height[i], hcb[i], lcw[i], h, imperial_units );
    }
    
    return cwa;
}
 
 
//' @title crown_profile_v() Estimates the crown width for a tree for a vector of heights above grown
//' @name crown_profile_v
//'
//' @param fia            : int    | FIA species code
//' @param dbh            : double | diameter at breast height
//' @param height         : double | total height
//' @param hcb            : double | height to crown base
//' @param lcw            : double | largest crown width
//' @param h              : double | vector of height compute crown width
//' @param imperial_units : bool   | TRUE = imperial (default), FALSE = metric
//'
//' @description
//' Estimates the crown width for a vector of heights above ground. The function uses the equation for relative crown shape:
//'  
//' \eqn{cw_h = lcw (1 - z^x)^(1/y)}
//' 
//' with z running from 0 at the widest point in the crown scaled by lcw to 1 at the tip, x and y are species-specific parameters, 
//' and a per species rule for where the widest point sits. Any specie without its own profile parameters (currently only Douglas-fir, 
//' western hemlock, and red alder have parameters) falls back to a softwood or hardwood default.
//'
//' @return
//' Returns a numeric vector with crown width at each height h in the original order (unsorted). Zero is returned for values of h
//' that are below a tree's height to crown base.
//'
//' @examples
//' # compute crown width at 
//' fia <- 202
//' dbh <- 10
//' height <- 60
//' hcb <- 20
//' cr <- (1-hcb)/height
//' lcw <- lcw( fia, dbh, cr )
//' h <- seq( 0, height, 1)
//' 
//' crown_profile_v( fia, dbh, height, hcb, lcw, h )
//'

 //' @export
 // [[Rcpp::export]]
 std::vector<double> crown_profile_v( 
         const int    &fia, 
         const double &dbh, 
         const double &height,
         const double &hcb, 
         const double &lcw,
         const std::vector<double> &h, 
         const bool imperial_units = true )
 {
     std::vector<double> cwa( h.size(), 0.0 );
     
     for( size_t i = 0; i < h.size(); ++i )
     {
         cwa[i] = compute_crown_profile( fia, dbh, height, hcb, lcw, h[i], imperial_units );
     }
     
     return cwa;
 }

