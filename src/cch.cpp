

#include "utilities.hpp"
#include <cmath>

// Predict crown width at a relative position up the tree using the functional form:
//    adjustment = rp^(b0 + b1*rp^0.5 + b2*height/dbh)
// where b0,b1,b2 are species-specific parameters and rp is relative position in the crown.
// The adjustment is applied to the largest crown width of the tree
double _cwa( const double rp, 
             const double lcw, 
             const double dbh, 
             const double height, 
             const std::vector<double> &parameters ) 
{
    double cwa = 0.0;

    double alpha = (parameters[0] + parameters[1] * std::sqrt(rp) + parameters[2] * (height / std::fmax( dbh, 0.1 )) ) ;

    // Compute adjusted crown width. Check alpha to prevent overflows, adjustment can be > 1.0 for small trees 
    cwa =  alpha > 0.0  ? lcw * std::pow( rp, alpha ) : 1.0;

    return cwa;
}


// compute crown closure at tree tip (cch) for a tree
double compute_cch( const std::vector<int>    &species,
                    const double ht, // height of tree to compute cch for
                    // attributes of trees in the stand
                    const std::vector<double> &dbh,
                    const std::vector<double> &height,
                    const std::vector<double> &crown_length, 
                    const std::vector<double> &dacb,         // distance above crown base to largest crown width (lcw)
                    const std::vector<double> &lcw,          // largest crown width of each tree
                    const std::vector<double> &expansion,
                    const std::unordered_map<int,std::vector<double>> &parameters,   // three parameters of equation describing crown shape from base to tip by species
                    const bool imperial )                    // (see _cwa() function for details.)
{
    size_t n = species.size();

    if( n == 0 || dbh.size() != n || height.size() != n || crown_length.size() != n || dacb.size() != n || lcw.size() != n ||
        expansion.size() != n || parameters.size() == 0 )
        return NAN;

    double cch = 0.0;
    double cw;
    double rp;
    double last_cw = -999.0 ;
    double last_dbh = -999.0;
    double last_height = -999.0;
    double last_crown_length = -999.0;

    std::vector<double> default_parameters = {1.0,0.0,0.0};

    // crown area to percentage of acre constant
    const double area_conversion = 0.25 * PI / ( imperial ? 43560.0 : 10000.0 );

    // Create ordered list decreasing by height for tree attributes
    auto idx = sort_indices( height );

    // Iterate through trees in decreasing height order
    for( auto i : idx )
    {
        // if the current tree's height is less than or equal to the search height (ht), return the accumulated cch
        if( ht >= height[i] )
            return cch;

        // if the current tree's height is the same as the tree before it in height, dbh, and crown length, just accumulate cch
        if( height[i] == last_height && dbh[i] == last_dbh && crown_length[i] == last_crown_length )
        {
            cch += (last_cw * last_cw) * (area_conversion * expansion[i]);
            continue;
        }

        // remember the current tree
        last_dbh = dbh[i];
        last_height = height[i];
        last_crown_length = crown_length[i];

        // compute height to largest crown width (hlcw) of current tree
        double htlc = height[i] - crown_length[i];
        double hlcw = htlc + dacb[i];

        // compute crown width at the search height (ht)
        if( ht <= hlcw )
        {
            // search height is below or equal to hlcw so use the lcw
            cw = lcw[i];
        }
        else if( ht > hlcw && ht < height[i] )
        {
            // search height is above hlcw. Compute crown width above lcw

            // compute relative position in crown
            rp = (height[i] - ht) / (height[i] - hlcw);

            // compute adjusted crown width (cwa is a function to compute crown width at a relative position in the crown)
            try {
                cw = _cwa( rp, lcw[i], dbh[i], height[i], parameters.at(species[i]) );
            } catch( ... ) {
                cw = _cwa( rp, lcw[i], dbh[i], height[i], default_parameters );
            }
        }
        else
        {
            // at the tip of the tree, use zero crown width
            cw = 0.0;
        }

        // accumulate compute crown area at the search position
        cch += (cw * cw) * (area_conversion * expansion[i]);
        last_cw = cw;
    }

    return cch > 0.0 ? cch : 0.0;
}
